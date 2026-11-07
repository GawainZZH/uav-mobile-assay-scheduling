"""Full physical route MILP: no restricted pool may certify the full problem."""
import math,time
from dataclasses import dataclass
from itertools import combinations
import numpy as np
from fleet_sample_assay.physical import enumerate_routes,vehicle_arcs
from fleet_sample_assay.master import Matrix
from fleet_sample_assay_certified.master import highs,pass_lp,native_run
from fleet_sample_assay_certified.network import assay_network
from fleet_sample_assay.check import integer_schedule


@dataclass
class PatternModel:
    M: object
    A: object
    routes: tuple
    arcs: tuple
    zids: list
    rids: list
    theta: int
    assignments: dict


def check_clock(deadline):
    if time.perf_counter()>=deadline:raise TimeoutError('complete pattern total budget')


def unstrengthened_build(inst,routes,deadline):
    """Build an aggregated supply network; caller must establish completeness."""
    check_clock(deadline);M=Matrix();arcs=vehicle_arcs(inst)
    for g in range(inst.groups):
        check_clock(deadline)
        for p in range(len(inst.pads)):
            for t in range(inst.H+1):
                M.row(('flow',g,p,t),int(p==0 and t==0)-int(p==0 and t==inst.H),True)
            for t in range(inst.H):M.row(('launch',g,p,t));M.row(('recover',g,p,t))
        for k in range(inst.H):
            M.row(('drone',g,k),inst.drones);M.row(('battery',g,k),inst.batteries)
            M.row(('slot',g,k))
    for j in range(len(inst.jobs)):M.row(('cover',j),1,True)
    for S in combinations(range(len(inst.jobs)),3):M.row(('sri',)+S,1)
    variants=sorted({(r.g,j,s,r.recovery) for r in routes for j,s in r.visits})
    for v in variants:M.row(('supply',)+v,0,True)
    M.row(('assay_cost',),0,True)
    zids=[];rids=[];assignments={}
    for g,p,t,q,u,cost in arcs:
        col={('flow',g,p,t):1,('flow',g,q,u):-1,
             ('launch',g,p,t):-inst.drones,('recover',g,p,t):-inst.drones}
        if p==q:col['slot',g,t]=-inst.stations
        zids.append(M.column(col,cost,upper=1,integer=True))
    for index,r in enumerate(routes):
        if index%128==0:check_clock(deadline)
        col={('cover',j):1 for j,s in r.visits}
        col['launch',r.g,r.p,r.launch]=1;col['recover',r.g,r.q,r.recovery]=1
        for k in r.drone_slots():col['drone',r.g,k]=1
        for k in r.battery_slots(inst):col['battery',r.g,k]=1
        for j,s in r.visits:col['supply',r.g,j,s,r.recovery]=-1
        visited=[j for j,s in r.visits];triples=set()
        for pair in combinations(visited,2):
            for j in range(len(inst.jobs)):
                if j not in pair:triples.add(tuple(sorted(pair+(j,))))
        for S in triples:col[('sri',)+S]=1
        rids.append(M.column(col,r.cost(inst),upper=1,integer=True))
    theta=M.column({('assay_cost',):1},1,upper=inst.H*sum(j.weight for j in inst.jobs))
    for g,j,s,e in variants:
        check_clock(deadline)
        end=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
        for k in range(e,end+1):
            assignments[g,j,s,e,k]=M.column({('supply',g,j,s,e):1,('slot',g,k):1,
                                           ('assay_cost',):-inst.jobs[j].weight*(k+1)},upper=1)
    check_clock(deadline);A=M.sparse().tocsc();check_clock(deadline)
    # Indexed sparse matrix is sufficient for validation; release dictionary copies.
    M.columns.clear()
    return PatternModel(M,A,tuple(routes),arcs,zids,rids,theta,assignments)


def build(inst,routes,deadline):
    from fleet_sample_assay_queue.bound import strengthen_pattern
    return strengthen_pattern(inst,unstrengthened_build(inst,routes,deadline))


def solve(inst,seconds=120,relax=False):
    start=time.perf_counter();deadline=start+seconds
    out=dict(status='enumeration_time_limit',upper_bound=None,lower_bound=0,objective=None,
             enumeration_complete=False,route_enumeration_calls=0,route_count=None,
             enumeration_seconds=0.0,build_seconds=0.0,solve_seconds=0.0,restore_seconds=0.0,
             schedule_check=None,seconds=0.0,native_lower_bound=None,
             candidate_completion_seconds=None,candidate_within_budget=None,late_validation=False)
    stage='enumeration';began=start
    try:
        check_clock(deadline);out['route_enumeration_calls']=1
        routes=enumerate_routes(inst,deadline)
        out.update(enumeration_complete=True,route_count=len(routes),enumeration_seconds=time.perf_counter()-began)
        stage='build';began=time.perf_counter();model=build(inst,routes,deadline)
        M=model.M;A=model.A
        c=np.array(M.cost);upper=np.array(M.upper);rhs=np.array(M.rhs);eq=np.array(M.eq,dtype=bool)
        if not all(np.all(np.isfinite(x)) for x in (A.data,c,upper,rhs)):
            raise ValueError('nonfinite complete pattern data')
        h=highs();pass_lp(h,A,c,upper,np.where(eq,rhs,-np.inf),rhs,
                          None if relax else np.array(M.integer))
        out.update(build_seconds=time.perf_counter()-began,rows=A.shape[0],columns=A.shape[1],
                   nonzeros=A.nnz,assay_variables=len(model.assignments),
                   binary_variables=sum(M.integer),native_backend=h.version())
        stage='solve';began=time.perf_counter()
        reserve=0 if relax else min(2.0,.1*seconds)
        status,sol,info,message=native_run(h,deadline-reserve)
        out.update(solve_seconds=time.perf_counter()-began,native_status=status,native_message=message,
                   nodes=int(info.mip_node_count) if not relax else None)
        if not relax:
            raw=float(info.mip_dual_bound)
            out['native_lower_bound']=raw if math.isfinite(raw) else None
            out['native_nonfinite_lower']=None if math.isfinite(raw) else str(raw)
            if math.isfinite(raw):out['lower_bound']=max(0,raw)
        if status==2:
            out['status']='lp_infeasible' if relax else 'infeasible'
            return out
        out['status']='lp_time_limit' if relax else 'time_limit'
        if not sol.value_valid:return out
        x=np.array(sol.col_value);activity=A@x;objective=float(c@x)
        error=max(np.max(np.abs(activity[eq]-rhs[eq]),initial=0),
                  np.max(np.maximum(0,activity[~eq]-rhs[~eq]),initial=0),
                  np.max(np.maximum(0,-x),initial=0),np.max(np.maximum(0,x-upper),initial=0),
                  abs(objective-info.objective_function_value))
        if not np.all(np.isfinite(x)) or not math.isfinite(objective) or error>(1e-7 if relax else 1e-6):
            raise ValueError(('pattern primal/objective residual',error))
        out.update(native_primal_error=float(error),native_objective=objective)
        if relax:
            if status==0:
                if not sol.dual_valid:raise ValueError('complete LP missing dual')
                pi=np.array(sol.row_dual);pi[~eq]=np.minimum(0,pi[~eq])
                if not np.all(np.isfinite(pi)):raise ValueError('complete LP nonfinite dual')
                bound=float(rhs@pi+np.minimum(0,c-A.T@pi)@upper)
                if not math.isfinite(bound) or bound>objective+1e-7:raise ValueError('complete LP crossed bound')
                out.update(status='lp_optimal',objective=objective,lower_bound=max(0,bound),
                           lp_certificate_slack=objective-bound)
            return out
        integer=np.array(M.integer,dtype=bool)
        if np.max(np.abs(x[integer]-np.rint(x[integer])),initial=0)>1e-6:
            raise ValueError('fractional pattern MIP incumbent')
        stage='restore';began=time.perf_counter()
        selected=[r for r,idx in zip(routes,model.rids) if x[idx]>.5]
        z=np.rint(x[model.zids]);qcap=np.zeros((inst.groups,inst.H))
        for a,value in zip(model.arcs,z):
            if a[1]==a[3]:qcap[a[0],a[2]]+=inst.stations*value
        net=assay_network(inst,selected,np.ones(len(selected)),qcap,deadline)
        if net['status']!=0 or net['error_bound']>1e-7:
            out['status']='restoration_unresolved';return out
        rounded=[]
        for edge,value in net['assignments']:
            if abs(value-round(value))>1e-7:raise ValueError('pattern assay matching not integral')
            if round(value):rounded.append((edge,round(value)))
        net={**net,'assignments':rounded}
        improved=sum(a[5]*v for a,v in zip(model.arcs,z))+sum(r.cost(inst) for r in selected)+net['cost']
        if improved>objective+1e-6 or (status==0 and abs(improved-objective)>1e-6):
            raise ValueError('pattern matching objective inconsistent')
        op=dict(status=0,objective=improved,z_values=z,arcs=model.arcs,route_values=np.ones(len(selected)))
        checked=integer_schedule(inst,selected,op,net,tol=1e-6)
        completion=time.perf_counter()-start
        if out['lower_bound']>improved+1e-6:raise ValueError('pattern bounds crossed after recovery')
        out.update(upper_bound=improved,objective=improved,schedule_check=checked,
                   candidate_completion_seconds=completion,candidate_within_budget=completion<=seconds,
                   _incumbent=dict(routes=selected,operation=op,network=net))
        if status==0:
            if math.ceil(out['lower_bound']-1e-7)!=round(improved):raise ValueError('pattern false optimal label')
            out['status']='optimal'
        return out
    except TimeoutError:
        out['status']=stage+'_time_limit'
        return out
    finally:
        elapsed=time.perf_counter()-began
        if stage=='enumeration':out['enumeration_seconds']=elapsed
        elif stage=='build':out['build_seconds']=elapsed
        elif stage=='solve' and not out['solve_seconds']:out['solve_seconds']=elapsed
        elif stage=='restore':out['restore_seconds']=elapsed
        out['seconds']=time.perf_counter()-start
        out['late_validation']=out['seconds']>seconds
        if out['status'] in ('optimal','lp_optimal','infeasible','lp_infeasible') and out['seconds']>seconds:
            out['status']='validation_time_limit'
