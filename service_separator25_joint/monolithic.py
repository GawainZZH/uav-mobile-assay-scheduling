"""Full same-domain MILP with the shared optimality-preserving separator rule."""
from contextlib import contextmanager
from review88_physical import milp as base
from .interface import certify,pack_limits
from .cases import split


@contextmanager
def strengthening(inst,boundaries,metadata):
    saved=base.build
    from .model import build as original
    def build(data,policy='J',deadline=None):
        m=original(data,policy,deadline)
        cert=certify(data,boundaries)
        cuts=[c['epoch'] for c in cert['cuts'] if c['valid']]
        removed=0
        for (g,i),col in m.zids.items():
            p,t,q,u=m.arcs[i]
            if p!=q and any(t<c<u for c in cuts):
                m.M.upper[col]=0;removed+=1
        effective=(0,*cuts,data.H);rows=0
        for stage,(local,ids) in enumerate(split(data,effective)):
            lo,hi=pack_limits(local);selected=set(ids)
            for g in range(data.groups):
                key=('stage_pack_max',g,stage);m.M.row(key,hi);rows+=1
                for i,r in enumerate(m.routes):
                    if r.job in selected:m.M.columns[m.rids[g,i]][key]=1
        m.A=m.M.sparse().tocsc()
        metadata.update(certification=cert,zeroed_cross_boundary_arcs=removed,stage_pack_rows=rows,
                        domains='same integer optimum; parked-cut dominance, not all-calendar identity')
        return m
    base.build=build
    try:yield
    finally:base.build=saved


def solve(inst,boundaries,seconds=60,policy='J',checkpoint=None):
    if policy not in ('J','QC'):raise ValueError('SCIP adapter supports J/QC; no implicit TP face')
    from .scip_backend import Proxy
    from fleet_sample_assay_certified.master import native_run
    import numpy as np
    import time,math
    began=time.perf_counter();deadline=began+seconds;metadata={}
    r=dict(method='PQ+',solver='SCIP',policy=policy,status='unresolved',optimal=False,
           objective=None,lower_bound=None,fleet=None,stages=[],native_calls=0)
    with strengthening(inst,boundaries,metadata):
        try:
            m=base.build(inst,policy,deadline-1);c=np.array(m.M.cost,dtype=float)
            rhs=np.array(m.M.rhs,dtype=float);eq=np.array(m.M.eq,dtype=bool)
            proxy=Proxy(m.A,c,np.array(m.M.upper),np.where(eq,rhs,-np.inf),rhs,np.array(m.M.integer))
            stage=dict(stage='joint',matrix=base.dump_matrix(m,c,None));r['stages'].append(stage)
            if checkpoint:checkpoint('prepared',r)
            status,sol,info,msg=native_run(proxy,deadline-1);r['native_calls']=1
            bound=float(info.mip_dual_bound);stage.update(status=status,message=msg,nodes=int(info.mip_node_count),
                lower_bound=bound if math.isfinite(bound) and abs(bound)<1e18 else None)
            if sol.value_valid:
                from .model import complete
                native_vector=np.array(sol.col_value)
                x=complete(m,native_vector,(0,)*inst.groups,(0,)*inst.groups);activity=m.A@x
                error=max(np.max(abs(activity[eq]-rhs[eq]),initial=0),np.max(activity[~eq]-rhs[~eq],initial=0),
                    np.max(-x,initial=0),np.max(x-np.array(m.M.upper),initial=0),np.max(abs(x-np.rint(x)),initial=0))
                if error>1e-6:raise ValueError(('SCIP matrix residual',error))
                fleet,physical=base.restore(m,x);value=float(physical['objective'])
                if abs(c@x-value)>1e-6:raise ValueError('SCIP physical objective')
                closed=status==0 and abs(bound-value)<=1e-6
                stage.update(matrix=base.dump_matrix(m,c,x),objective=value,fleet=fleet,optimal=closed,
                    native_vector=native_vector.tolist())
                r.update(status='optimal' if closed else 'feasible',objective=value,lower_bound=bound,
                    fleet=fleet,optimal=closed,proof_scope='native_numerical_integer_optimum' if closed else None)
            elif status==2:r.update(status='infeasible')
        except TimeoutError as e:r['reason']=str(e)
        except Exception as e:r.update(status='correctness_error',reason=repr(e))
    r['seconds']=time.perf_counter()-began
    if r['seconds']>seconds:r.update(status='workflow_time_limit',optimal=False,proof_scope=None)
    r.update(method='PQ+',shared_strengthening=metadata)
    return r
