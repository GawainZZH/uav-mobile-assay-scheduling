"""Original network certificate; at most one stricter re-solve on bad error."""
import time
import highspy
from types import SimpleNamespace
import numpy as np
from scipy.sparse import coo_matrix,vstack
from fleet_sample_assay import master
from fleet_sample_assay.guarded import _dump,TOL
from fleet_sample_assay_certified.master import highs,pass_lp,native_run,check

RETRIES=[]
NATIVE_LP_CALLS=0
LP_ATTEMPTS=0

def reset():
    global NATIVE_LP_CALLS,LP_ATTEMPTS
    RETRIES.clear();NATIVE_LP_CALLS=0;LP_ATTEMPTS=0

def _linprog(c,Aeq,beq,Aub,bub,bounds,seconds,strict=False,physical_mass=1):
    global NATIVE_LP_CALLS,LP_ATTEMPTS
    LP_ATTEMPTS+=1
    if not strict:
        NATIVE_LP_CALLS+=1
        return master._linprog(c,Aeq,beq,Aub,bub,bounds,seconds)
    deadline=time.perf_counter()+seconds;h=highs()
    check(h.setOptionValue('primal_feasibility_tolerance',1e-10))
    check(h.setOptionValue('dual_feasibility_tolerance',1e-10))
    if bounds!=(0,None):raise ValueError('network nonnegative variable bounds required')
    if not np.isfinite(physical_mass) or physical_mass<=0:raise ValueError('strict network positive finite mass required')
    rho=physical_mass/min(physical_mass,1)
    A=vstack((Aeq,Aub),format='csc');lower=np.r_[np.array(beq)*rho,np.full(len(bub),-np.inf)];upper=np.r_[np.array(beq)*rho,np.array(bub)*rho]
    pass_lp(h,A,np.asarray(c),np.full(len(c),np.inf),lower,upper)
    left=deadline-time.perf_counter()
    if left<=0:raise TimeoutError('strict network build exhausted budget')
    check(h.setOptionValue('time_limit',h.getRunTime()+left))
    NATIVE_LP_CALLS+=1;check(h.run())
    native_status=h.getModelStatus()
    status=0 if native_status==highspy.HighsModelStatus.kOptimal else (2 if native_status==highspy.HighsModelStatus.kInfeasible else 1)
    sol=h.getSolution();info=h.getInfo()
    return SimpleNamespace(status=status,x=np.array(sol.col_value)/rho if sol.value_valid else None,
        fun=info.objective_function_value/rho if sol.value_valid else None,
        ineqlin=SimpleNamespace(marginals=np.array(sol.row_dual)[len(beq):] if sol.dual_valid else None))

def _network(inst,routes,values,q,deadline=None,strict=False):
    start=time.perf_counter();supply={}
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(q)):raise ValueError('nonfinite network input')
    if np.min(q,initial=0)<-TOL:raise ValueError('negative slot capacity')
    for r,x in zip(routes,values):
        if x<-TOL:raise ValueError('negative supply')
        if x<=0:continue
        for j,s in r.visits:
            key=r.g,j,s,r.recovery;supply[key]=supply.get(key,0)+float(x)
    if not supply:return dict(status=0,cost=0.0,beta=master.zero_cut(inst),assignments=[],error_bound=0.0,seconds=time.perf_counter()-start)
    variants=sorted(supply);U=sum(supply.values());edges=[];rows=[];cols=[];caprows=[];c=[]
    for i,(g,j,s,e) in enumerate(variants):
        d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
        if d<e:raise ValueError('empty network window')
        for k in range(e,d+1):
            idx=len(edges);edges.append((g,j,s,e,k));rows.append(i);cols.append(idx);caprows.append(g*inst.H+k);c.append(inst.jobs[j].weight*(k+1))
    Aeq=coo_matrix((np.ones(len(edges)),(rows,cols)),shape=(len(variants),len(edges))).tocsr()
    Aub=coo_matrix((np.ones(len(edges)),(caprows,cols)),shape=(inst.groups*inst.H,len(edges))).tocsr()
    b=np.array([supply[v]/U for v in variants])
    # Scaling by total supply protects all-small inputs. Capacity above twice
    # total demand is redundant and strictly slack, so its dual remains zero.
    cap=np.minimum(np.maximum(0,np.asarray(q,dtype=float))/U,2).ravel()
    left=60 if deadline is None else deadline-time.perf_counter()
    if left<=0:raise TimeoutError('guarded network budget')
    result=_linprog(c,Aeq,b,Aub,cap,(0,None),left,strict=strict,physical_mass=U)
    _dump('network',dict(instance=inst,variants=variants,supply=[supply[v] for v in variants],total_supply=U,
        edges=edges,costs=c,normalized_capacity=cap,original_capacity=q,native_x=result.x,native_objective=result.fun,
        status=int(result.status),native_capacity_dual=result.ineqlin.marginals))
    rec=dict(status=int(result.status),cost=None,beta=None,error_bound=None,seconds=time.perf_counter()-start)
    if result.status!=0:return rec
    if not np.all(np.isfinite(result.x)) or not np.isfinite(result.fun) or not np.all(np.isfinite(result.ineqlin.marginals)):
        raise ValueError('nonfinite network native point/dual')
    beta=np.minimum(0,np.array(result.ineqlin.marginals)).reshape(inst.groups,inst.H)
    if np.max(result.ineqlin.marginals)>TOL:raise ValueError('network capacity dual sign')
    primal_error=float(U*(np.sum(np.abs(Aeq@result.x-b))+np.sum(np.maximum(0,Aub@result.x-cap))+np.sum(np.maximum(0,-result.x))))
    if np.min(result.x,initial=0)<-TOL:raise ValueError('negative network native assignment')
    alpha=[master.alpha(inst,beta,*v) for v in variants]
    dual_value=sum(a*supply[v] for a,v in zip(alpha,variants))+float(np.sum(beta*q))
    cost=float(np.dot(c,result.x)*U)
    objective_error=abs(cost-result.fun*U)
    error_bound=max(abs(cost-dual_value),max(c)*primal_error,objective_error)
    rec.update(cost=cost,beta=tuple(tuple(float(v) for v in g) for g in beta),error_bound=error_bound,
        total_supply=U,smallest_supply=min(supply.values()),primal_l1_error=primal_error,
        assignments=[(edge,float(x*U)) for edge,x in zip(edges,result.x) if x>0])
    return rec


def assay_network(inst,routes,values,q,deadline=None):
    began=time.perf_counter();first=_network(inst,routes,values,q,deadline,strict=False)
    if first['status']!=0 or first['error_bound']<=TOL:return first
    event=dict(first_error=first['error_bound'],first_seconds=first['seconds'],retry_status=None,retry_error=None,
               retry_seconds=0.0,success=False)
    RETRIES.append(event);start=time.perf_counter()
    try:
        second=_network(inst,routes,values,q,deadline,strict=True)
        event.update(retry_status=second['status'],retry_error=second['error_bound'],
                     success=second['status']==0 and second['error_bound']<=TOL)
        second['seconds']=time.perf_counter()-began
        return second
    finally:event['retry_seconds']=time.perf_counter()-start
