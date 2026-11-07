"""Finite logical decomposition prototype; master solves use branch-and-bound.

This is a repeated-master implementation, not a single-tree branch-and-check
handler. Only verified complete fleets are global incumbents. Unresolved local
queries never justify infeasibility or acceptance of a master-only incumbent.
"""
from fractions import Fraction as F
import math,time
import numpy as np

from service_separator25_joint.solve import context,concatenate
from .contracts import histogram
from .common import compact,project_actions
from .matrix import encode,decode
from .native import solve_mip


def record_entries(d,objects):
    records=[]
    for key,k in objects.items():
        floor=min(d.floors[i]for i,t in enumerate(d.templates)if t==key)
        for ctx,r in sorted(k.solved.items()):
            pattern=histogram(zip(*ctx))
            if r['status']=='infeasible':records.append(dict(template=key,pattern=pattern,infeasible=True))
            else:
                bound=F(str(k.lower(*ctx)))
                if r.get('objective')is not None and bound>F(str(r['objective'])):
                    raise ValueError('local bound exceeds physical context incumbent')
                if bound>floor:records.append(dict(template=key,pattern=pattern,bound=bound))
    return records


def edges_for(d,actions):
    path=[]
    for k,stage in enumerate(actions):
        start=tuple(a[0]for a in stage);end=tuple(a[1]for a in stage);used=tuple(a[2]for a in stage)
        ctx,perm=context(start,end,used,True)
        path.append(dict(stage=k,kernel=d.templates[k],context=ctx,perm=perm,reorder=tuple(range(d.groups))))
    return tuple(path)


def solve_prepared(d,objects,source,anchor,deadline,*,master_solver=solve_mip):
    began=time.perf_counter();actions,parts=project_actions(d,source)
    best=sum(F(p['objective'])for p in parts);fleet=source;lower=F(anchor)
    if lower>best:raise ValueError('global bound exceeds validated physical incumbent')
    records=record_entries(d,objects);masters=[];queries=[];deliveries=[];reason=None
    def refresh(path):
        nonlocal best,fleet
        if all(objects[e['kernel']].solved.get(e['context'],{}).get('objective')is not None for e in path):
            restored=concatenate(d.inst,d.stages,objects,path)
            value=F(restored['check']['objective'])
            if value<best:
                best=value;fleet=restored;deliveries.append(dict(objective=best,fleet=fleet))
            if lower>best:raise ValueError('physical completion contradicts master lower bound')
    def result():
        return dict(status='optimal'if lower==best else'feasible',optimal=lower==best,
            objective=best,lower_bound=lower,fleet=fleet,cuts=records,
            master_calls=len(masters),local_calls=len(queries),master_events=masters,
            local_events=queries,deliveries=deliveries,reason=reason,seconds=time.perf_counter()-began,
            single_tree=False,global_upper_basis='physically validated complete fleet')
    # The same validated-price closure applies to the monolithic comparator.
    if lower==best:return result()
    for iteration in range(128):
        if time.perf_counter()>=deadline:reason='workflow_allowance';break
        try:master=compact(d,best,records=records,anchor=lower,deadline=deadline)
        except TimeoutError:reason='master_construction_allowance';break
        # The current physical incumbent is a valid master start under all valid cuts.
        source_actions,source_parts=project_actions(d,fleet)
        start=encode(master,source_actions,tuple(p['objective']for p in source_parts))
        try:out=master_solver(master,min(deadline,time.perf_counter()+10),initial=start)
        except TimeoutError:reason='master_allowance';break
        event={k:v for k,v in out.items()if k!='vector'}
        event.update(shape=list(master.A.shape),nonzeros=master.A.nnz,
            matrix=_dump(master,out.get('vector')))
        masters.append(event)
        if out['status']=='infeasible':raise ValueError('master infeasibility contradicts physical incumbent')
        if out.get('bound')is not None:lower=max(lower,F(math.ceil(out['bound']-1e-6)))
        if lower>best:raise ValueError('master bound exceeds physical incumbent')
        if lower==best:break
        if out.get('vector')is None:reason='master_unresolved';break
        decoded=decode(master,out['vector']);path=edges_for(d,decoded['actions'])
        refresh(path)
        if lower==best:break
        chosen=None
        for e,theta in zip(path,decoded['costs']):
            k=objects[e['kernel']];r=k.solved.get(e['context'],{})
            if not r.get('optimal',False)and (r.get('objective')is None or theta+1e-6<r['objective']):
                chosen=e;break
        if chosen is None:reason='master_bound_open';break
        k=objects[chosen['kernel']];ctx=chosen['context']
        before_lower=F(str(k.lower(*ctx)));before_upper=k.solved.get(ctx,{}).get('objective')
        q=time.perf_counter()
        try:r=k.solve(*ctx,min(deadline,q+5))
        except TimeoutError:
            records=record_entries(d,objects);refresh(path);reason='local_allowance';break
        queries.append(dict(
            kernel=chosen['kernel'],context=ctx,seconds=time.perf_counter()-q,result=r))
        records=record_entries(d,objects)
        refresh(path)
        if lower==best:break
        # A proof, stronger lower bound or cheaper feasible schedule must advance
        # this candidate. An exhausted local quantum is never a no-good cut.
        after_lower=math.inf if r['status']=='infeasible'else float(k.lower(*ctx))
        after_upper=r.get('objective')
        progress=(r['status']=='infeasible'or r.get('optimal',False)or after_lower>float(before_lower)+1e-6 or
                  after_upper is not None and (before_upper is None or after_upper<before_upper))
        if not progress:reason='local_unresolved';break
    else:reason='finite_prototype_iteration_cap'
    return result()


def _dump(data,vector):
    from review88_physical.milp import dump_matrix
    return dump_matrix(data,np.asarray(data.M.cost,dtype=float),vector)
