"""Independent group-indexed complete binary calendar formulation."""
import math
import time
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix, vstack
from fleet_sample_assay.master import Matrix
from fleet_sample_assay_certified.master import highs, pass_lp, native_run


def clock(deadline):
    if deadline is not None and time.perf_counter()>=deadline:
        raise TimeoutError('Review88 complete workflow deadline')


def independent_routes(inst,deadline=None):
    from .types import Route
    result=[];P=len(inst.pads);q=inst.assay_duration;ds=inst.sampling_duration
    for j,job in enumerate(inst.jobs):
        for s in range(job.lo,job.hi+1):
            clock(deadline)
            latest=min(inst.H-q,s+job.shelf,job.deadline-q)-inst.handoff_delay
            for p in range(P):
                for b in range(P):
                    for launch in range(max(0,s-inst.endurance),s-inst.fly(p,P+j)+1):
                        earliest=s+ds+inst.fly(P+j,b)
                        last=min(launch+inst.endurance,latest,inst.H-inst.ground_duration)
                        for end in range(earliest,last+1):result.append(Route(p,launch,b,end,j,s))
    return tuple(result)


def build(inst,policy='J',deadline=None):
    if policy not in ('J','TP','QC'):raise ValueError('unknown policy')
    M=Matrix();H=inst.H;P=len(inst.pads);q=inst.assay_duration
    routes=independent_routes(inst,deadline)
    arcs=tuple((p,t,b,t+dur) for p in range(P) for t in range(H)
               for b in range(P) if t+(dur:=1 if b==p else inst.road(p,b))<=H)
    variants=sorted({(r.job,r.sample,r.recovery) for r in routes})
    for g in range(inst.groups):
        for p in range(P):
            for t in range(H+1):M.row(('flow',g,p,t),int(p==0 and t==0)-int(p==0 and t==H),True)
            for t in range(H):M.row(('launch',g,p,t));M.row(('return',g,p,t))
        for t in range(H):M.row(('drone',g,t),inst.drones);M.row(('slot',g,t))
        M.row(('pack',g),inst.batteries)
    for j in range(len(inst.jobs)):M.row(('cover',j),1,True)
    for g in range(inst.groups):
        for v in variants:M.row(('supply',g)+v,0,True)
        for i,r in enumerate(routes):
            for t in range(r.recovery,r.recovery+inst.ground_duration):M.row(('ground',g,i,t))
        if policy=='QC':
            for j in range(len(inst.jobs)):
                for t in range(H):M.row(('clear',g,j,t),1)
    ground={}
    for i,r in enumerate(routes):
        for t in range(r.recovery,r.recovery+inst.ground_duration):ground.setdefault((r.q,t),[]).append(i)
    zids={};rids={};assignments={}
    for g in range(inst.groups):
        clock(deadline)
        for ai,(p,t,b,u) in enumerate(arcs):
            col={('flow',g,p,t):1,('flow',g,b,u):-1,('launch',g,p,t):-inst.drones,('return',g,p,t):-inst.drones}
            if p==b:
                col['slot',g,t]=-inst.stations
                for i in ground.get((p,t),()):col['ground',g,i,t]=-1
            elif policy=='QC':
                for j in range(len(inst.jobs)):col['clear',g,j,t]=1
            zids[g,ai]=M.column(col,0 if p==b else inst.road_cost*(u-t),upper=1,integer=True)
        for i,r in enumerate(routes):
            if i%64==0:clock(deadline)
            col={('cover',r.job):1,('pack',g):1,('supply',g,r.job,r.sample,r.recovery):-1,
                 ('launch',g,r.p,r.launch):1,('return',g,r.q,r.recovery):1}
            for t in range(r.launch,r.recovery+inst.ground_duration):col['drone',g,t]=1
            for t in range(r.recovery,r.recovery+inst.ground_duration):col['ground',g,i,t]=1
            if policy=='QC':
                for t in range(r.recovery,H):col['clear',g,r.job,t]=1
            rids[g,i]=M.column(col,inst.flight_fixed+inst.air_cost*(r.recovery-r.launch),upper=1,integer=True)
        for j,s,e in variants:
            clock(deadline);job=inst.jobs[j]
            latest=min(H-q,s+job.shelf,job.deadline-q)
            for k in range(e+inst.handoff_delay,latest+1):
                col={('supply',g,j,s,e):1}
                for t in range(k,k+q):col['slot',g,t]=1
                if policy=='QC':
                    for t in range(k+q,H):col['clear',g,j,t]=-1
                assignments[g,j,s,e,k]=M.column(col,job.weight*(k+q-s),upper=1,integer=True)
    A=M.sparse().tocsc();transport=np.array(M.cost,dtype=float)
    for c in assignments.values():transport[c]=0
    clock(deadline)
    return SimpleNamespace(M=M,A=A,routes=routes,arcs=arcs,rids=rids,zids=zids,
                           assignments=assignments,transport=transport,inst=inst,policy=policy)


def restore(m,x):
    from .physical import validate_fleet
    plans=[]
    for g in range(m.inst.groups):
        selected=[r for i,r in enumerate(m.routes) if x[m.rids[g,i]]>.5]
        arcs=[a for i,a in enumerate(m.arcs) if x[m.zids[g,i]]>.5]
        trail=[];node=(0,0)
        while node!=(0,m.inst.H):
            nxt=[a for a in arcs if a[:2]==node]
            if len(nxt)!=1:raise ValueError('broken complete group truck path')
            a=nxt[0];trail.append(a);node=a[2:]
            if len(trail)>len(arcs):raise ValueError('truck cycle')
        if len(trail)!=len(arcs):raise ValueError('extra truck component')
        assays=[(j,s,e,k) for (gg,j,s,e,k),c in m.assignments.items() if gg==g and x[c]>.5]
        plan=dict(routes=selected,trail=trail,assays=assays)
        if m.policy=='QC':
            for p,t,b,u in trail:
                if p!=b and any(e<=t<k+m.inst.assay_duration for j,s,e,k in assays):
                    raise ValueError('queue-clear counterfactual violated')
        plans.append(plan)
    check=validate_fleet(m.inst,plans)
    if not check['valid']:raise ValueError(('independent fleet replay',check))
    return dict(group_plans=plans,check=check),check


def dump_matrix(m,c,x):
    a=m.A.tocoo()
    return dict(shape=list(a.shape),row=a.row.tolist(),col=a.col.tolist(),data=a.data.tolist(),
        cost=np.asarray(c,dtype=float).tolist(),rhs=list(m.M.rhs),equalities=list(m.M.eq),
        upper=list(m.M.upper),integer=list(m.M.integer),x=None if x is None else np.asarray(x).tolist())


def solve(inst,policy='J',seconds=30,checkpoint=None):
    start=time.perf_counter();deadline=start+seconds;reserve=min(1.,.1*seconds)
    if not math.isfinite(seconds) or seconds<0:raise ValueError('finite nonnegative budget')
    out=dict(method='PQ',policy=policy,status='unknown',stages=[],native_calls=0,fleet=None,
             lower_bound=None,objective=None,optimal=False,time_budget=seconds,proof_scope=None)
    try:
        m=build(inst,policy,deadline-reserve)
        out.update(route_count=len(m.routes),rows=m.A.shape[0],columns=m.A.shape[1],nonzeros=m.A.nnz,
                   construction_seconds=time.perf_counter()-start,
                   support=dict(routes=m.routes,arcs=m.arcs,zids=[(g,i,c) for (g,i),c in m.zids.items()],
                       rids=[(g,i,c) for (g,i),c in m.rids.items()],assignments=list(m.assignments.items())))
        for name in (('transport','joint') if policy=='TP' else ('joint',)):
            clock(deadline-reserve)
            c=m.transport if name=='transport' else np.array(m.M.cost,dtype=float)
            h=highs();rhs=np.array(m.M.rhs);eq=np.array(m.M.eq,dtype=bool)
            pass_lp(h,m.A,c,np.array(m.M.upper),np.where(eq,rhs,-np.inf),rhs,np.array(m.M.integer))
            rec=dict(stage=name,matrix=dump_matrix(m,c,None));out['stages'].append(rec)
            if checkpoint:checkpoint(name+'_prepared',out)
            out['native_calls']+=1;t=time.perf_counter()
            stop=min(deadline-reserve,start+min(12.,.4*seconds)) if name=='transport' else deadline-reserve
            status,sol,info,msg=native_run(h,stop);bound=float(info.mip_dual_bound)
            rec.update(status=status,message=msg,seconds=time.perf_counter()-t,nodes=int(info.mip_node_count),
                       lower_bound=bound if math.isfinite(bound) else None,optimal=False)
            if not sol.value_valid:
                out['status']='infeasible' if status==2 else 'time_limit';break
            x=np.asarray(sol.col_value);activity=m.A@x
            error=max(np.max(abs(activity[eq]-rhs[eq]),initial=0),np.max(activity[~eq]-rhs[~eq],initial=0),
                      np.max(-x,initial=0),np.max(x-np.asarray(m.M.upper),initial=0),
                      np.max(abs(x-np.rint(x)),initial=0),abs(float(c@x)-float(info.objective_function_value)))
            if not np.all(np.isfinite(x)) or error>1e-6:raise ValueError(('matrix residual',error))
            fleet,check=restore(m,x)
            value=float(check['road_cost']+check['flight_cost'] if name=='transport' else check['objective'])
            if abs(float(c@x)-value)>1e-6:raise ValueError('physical cost residual')
            closed=status==0 and math.isfinite(bound) and abs(bound-value)<=1e-6
            rec.update(matrix=dump_matrix(m,c,x),objective=value,optimal=closed,fleet=fleet,physical=check)
            if checkpoint:checkpoint(name+'_restored',out)
            if name=='transport':
                out['transport_optimum']=value if closed else None
                if not closed:out['status']='transport_unresolved';break
                m.M.row(('transport_face',),value,True)
                m.A=vstack((m.A,csr_matrix(c.reshape(1,-1))),format='csc')
            else:
                out.update(fleet=fleet,physical=check,objective=value,lower_bound=rec['lower_bound'],
                           optimal=closed,status='optimal' if closed else 'feasible',
                           proof_scope='native_numerical_integer_optimum' if closed else None)
        clock(deadline)
    except TimeoutError as exc:out.update(status='workflow_time_limit',reason=str(exc),optimal=False)
    except Exception as exc:out.update(status='correctness_error',reason=repr(exc),optimal=False)
    out['seconds']=time.perf_counter()-start
    if out['seconds']>seconds:
        if out['status']!='correctness_error':out['status']='workflow_time_limit'
        out.update(optimal=False,proof_scope=None,objective=None,fleet=None)
    return out
