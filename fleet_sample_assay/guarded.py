"""Numerical guards for the frozen first physical prototype; no source rewrite."""
import time
import json
from dataclasses import asdict,is_dataclass
from itertools import count
import numpy as np
from scipy.sparse import coo_matrix
from . import master

TOL=1e-7
RAW_DIR=None
_IDS=count()


def _serial(x):
    if is_dataclass(x):return asdict(x)
    if isinstance(x,dict):return {str(k):_serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [_serial(v) for v in x]
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    return x


def _dump(kind,obj):
    if RAW_DIR is not None:
        RAW_DIR.mkdir(parents=True,exist_ok=True)
        (RAW_DIR/(f'{next(_IDS):04d}_{kind}.json')).write_text(json.dumps(_serial(obj),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def validate_master(inst,routes,r,cuts,phase,explicit=False):
    if r['status']!=0:return r
    arrays=[r[k] for k in ('route_values','z_values','q','artificial')]
    scalar=[r['objective'],r['theta']]+list(r.get('dual',{}).values())
    if any(not np.all(np.isfinite(x)) for x in arrays+[scalar]):
        raise ValueError('nonfinite native point or dual')
    if len(r['route_values'])!=len(routes) or len(r['z_values'])!=len(r['arcs']):raise ValueError('point dimensions')
    if any(np.min(np.asarray(x),initial=0)<-TOL for x in arrays) or r['theta']<-TOL:raise ValueError('negative native primal')
    if np.max(r['z_values'],initial=0)>1+TOL:raise ValueError('vehicle upper bound')
    rows={};flow={};q=np.zeros((inst.groups,inst.H));cost=0.0
    for x,(g,p,t,qp,u,c) in zip(r['z_values'],r['arcs']):
        flow[g,p,t]=flow.get((g,p,t),0)+x;flow[g,qp,u]=flow.get((g,qp,u),0)-x
        rows['launch',g,p,t]=rows.get(('launch',g,p,t),0)-inst.drones*x
        rows['recover',g,p,t]=rows.get(('recover',g,p,t),0)-inst.drones*x
        if p==qp:q[g,t]+=inst.stations*x
        if phase==2:cost+=c*x
    for route,x in zip(routes,r['route_values']):
        for key,v in master.route_rows(inst,route,cuts,explicit).items():rows[key]=rows.get(key,0)+v*x
        if phase==2:cost+=route.cost(inst)*x
    for g in range(inst.groups):
        for p in range(len(inst.pads)):
            for t in range(inst.H+1):
                rhs=int(p==0 and t==0)-int(p==0 and t==inst.H)
                if abs(flow.get((g,p,t),0)-rhs)>TOL:raise ValueError('vehicle flow residual')
        for k in range(inst.H):
            if rows.get(('drone',g,k),0)>inst.drones+TOL or rows.get(('battery',g,k),0)>inst.batteries+TOL:
                raise ValueError('occupancy residual')
        if not explicit:
            for a in range(inst.H):
                for b in range(a,inst.H):
                    if rows.get(('hall',g,a,b),0)>sum(q[g,a:b+1])+TOL:raise ValueError('Hall residual')
    if any(v>TOL for key,v in rows.items() if key[0] in ('launch','recover')):raise ValueError('endpoint residual')
    for j,a in enumerate(r['artificial']):
        if abs(rows.get(('cover',j),0)+a-1)>TOL:raise ValueError('coverage residual')
        if phase==2 and abs(a)>TOL:raise ValueError('phase-two artificial')
    if not explicit:
        for n,beta in enumerate(cuts):
            if rows.get(('cost',n),0)+np.sum(np.array(beta)*q)-r['theta']>TOL:raise ValueError('cost residual')
    if np.max(np.abs(q-r['q']),initial=0)>TOL:raise ValueError('capacity reconstruction')
    expected=cost+r['theta'] if phase==2 else sum(r['artificial'])
    if abs(expected-r['objective'])>TOL:raise ValueError('objective residual')
    for key,value in r.get('dual',{}).items():
        if key[0] in ('launch','recover','drone','battery','hall','cost','slot') and value>TOL:
            raise ValueError('inequality dual sign')
    if 'dual' in r:
        inequality=('launch','recover','drone','battery','hall','cost','slot')
        dual={key:min(0,value) if key[0] in inequality else value for key,value in r['dual'].items()}
        r['dual']=dual;rcz=[]
        for g,p,t,qp,u,c in r['arcs']:
            col={('flow',g,p,t):1,('flow',g,qp,u):-1,('launch',g,p,t):-inst.drones,('recover',g,p,t):-inst.drones}
            if p==qp:
                if explicit:col[('slot',g,t)]=-inst.stations
                else:
                    for a in range(t+1):
                        for b in range(t,inst.H):col[('hall',g,a,b)]=-inst.stations
                    for n,beta in enumerate(cuts):col[('cost',n)]=inst.stations*beta[g][t]
            rcz.append((c if phase==2 else 0)-sum(dual[k]*v for k,v in col.items()))
        free_rc=[]
        for route in routes:
            free_rc.append((route.cost(inst) if phase==2 else 0)-sum(dual[k]*v for k,v in master.route_rows(inst,route,cuts,explicit).items()))
        theta_rc=(1 if phase==2 else 0)-dual[('assay_cost',)] if explicit else (1 if phase==2 else 0)+sum(dual[('cost',n)] for n in range(len(cuts)))
        route_rc=list(free_rc)
        free_rc.append(theta_rc)
        if phase==1:free_rc.extend(1-dual[('cover',j)] for j in range(len(inst.jobs)))
        if explicit:
            variants={(route.g,j,s,route.recovery) for route in routes for j,s in route.visits}
            for g,j,s,e in variants:
                d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
                free_rc.extend(-dual[('supply',g,j,s,e)]-dual[('slot',g,k)]+inst.jobs[j].weight*(k+1)*dual[('assay_cost',)] for k in range(e,d+1))
        if min(free_rc,default=0)<-TOL:raise ValueError('native dual reduced-cost residual')
        lower=sum(dual[('flow',g,0,0)]-dual[('flow',g,0,inst.H)] for g in range(inst.groups))
        lower+=sum(dual[('cover',j)] for j in range(len(inst.jobs)))
        lower+=sum(inst.drones*dual[('drone',g,k)]+inst.batteries*dual[('battery',g,k)] for g in range(inst.groups) for k in range(inst.H))
        lower+=sum(min(0,x) for x in rcz)
        # Finite-box correction, rather than silently treating tiny negative
        # RCs of nominally unbounded columns as dual feasible. Exact coverage
        # implies lambda<=1, each assay y<=1, artificial<=1. Theta<=H*sum(w)
        # holds in E and is an optimality-preserving bound in projected P:
        # Hall feasibility gives a matching whose cost is at most this value.
        theta_bound=inst.H*sum(j.weight for j in inst.jobs)
        lower+=sum(min(0,x) for x in route_rc)
        lower+=theta_bound*min(0,theta_rc)
        rest=free_rc[len(route_rc)+1:]
        lower+=sum(min(0,x) for x in rest)
        if abs(lower-r['objective'])>TOL:raise ValueError('full LP dual objective residual')
        r['verified_dual_bound']=lower
        r['finite_box_theta_bound']=theta_bound
    return r


def solve_master(inst,routes,**kwargs):
    start=time.perf_counter();r=master.solve_master(inst,routes,**kwargs)
    _dump('master',dict(instance=inst,routes=routes,arguments=kwargs,result=r))
    validate_master(inst,routes,r,kwargs.get('cuts',()),kwargs.get('phase',2),kwargs.get('explicit',False))
    if r['status']==0 and kwargs.get('explicit',False):
        net=assay_network(inst,routes,r['route_values'],r['q'],kwargs.get('deadline'))
        if net['status']!=0 or net['error_bound']>TOL or abs(net['cost']-r['theta'])>TOL:
            raise ValueError('explicit LP has no matching-cost substitute assay witness')
        r['assay_witness']=net
    r['seconds_with_validation']=time.perf_counter()-start
    return r


def assay_network(inst,routes,values,q,deadline=None):
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
    result=master._linprog(c,Aeq,b,Aub,cap,(0,None),left)
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
