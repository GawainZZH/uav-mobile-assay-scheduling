from dataclasses import replace
from fractions import Fraction
from itertools import product
import hashlib
import json
import math
import time
from review88_physical.evidence import plain
from review88_physical.physical import validate_fleet
from .cases import split
from .interface import certify,canonical,pack_limits
from .kernel import Kernel


def signature(inst):
    d=plain(inst);d.pop('name')
    for j in d['jobs']:j.pop('name')
    return hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest()


def context(start,end,consumed,pooling):
    if pooling:
        order=tuple(sorted(range(len(start)),key=lambda i:(start[i],end[i],consumed[i])))
    else:order=tuple(range(len(start)))
    return (tuple(start[i] for i in order),tuple(end[i] for i in order),tuple(consumed[i] for i in order)),order


def master(inst,stages,objects,*,pooling=True,upper=False,deadline=math.inf):
    initial=tuple((0,inst.batteries) for _ in range(inst.groups));labels={initial:(0,())}
    total_states=1;total_edges=0
    for stage,(local,ids,key) in enumerate(stages):
        if time.perf_counter()>=deadline:raise TimeoutError('interface graph budget')
        new={};k=objects[key];lo,hi=pack_limits(local);n=len(ids)
        for state,(value,path) in labels.items():
            if time.perf_counter()>=deadline:raise TimeoutError('interface expansion budget')
            for consumed in product(range(lo,hi+1),repeat=inst.groups):
                if sum(consumed)!=n or any(c>b for c,(p,b) in zip(consumed,state)):continue
                ends=product(range(len(inst.pads)),repeat=inst.groups) if stage+1<len(stages) else [(0,)*inst.groups]
                for end in ends:
                    start=tuple(p for p,b in state);ctx,perm=context(start,end,consumed,pooling)
                    result=k.solved.get(ctx)
                    if upper:cost=math.inf if result is None or result.get('objective') is None else result['objective']
                    else:cost=k.lower(*ctx)
                    if not math.isfinite(cost):continue
                    nxt=tuple((p,b-c) for p,c,(_,b) in zip(end,consumed,state))
                    if pooling:nxt,reorder=canonical(nxt)
                    else:reorder=tuple(range(inst.groups))
                    total_edges+=1
                    edge=dict(stage=stage,kernel=key,context=ctx,perm=perm,reorder=reorder,
                              state=state,next_state=nxt,consumed=consumed)
                    candidate=(value+cost,path+(edge,))
                    if nxt not in new or candidate[0]<new[nxt][0]:new[nxt]=candidate
        labels=new;total_states+=len(labels)
        if not labels:return math.inf,(),dict(states=total_states,edges=total_edges)
    value,path=min(labels.values(),key=lambda r:r[0])
    return value,path,dict(states=total_states,edges=total_edges)


def concatenate(inst,stages,objects,path):
    plans=[dict(routes=[],trail=[],assays=[]) for _ in range(inst.groups)]
    physical_order=list(range(inst.groups));offset=0
    for (local,ids,key),edge in zip(stages,path):
        result=objects[key].solved[edge['context']];fleet=result['fleet']
        for native,canonical_group in enumerate(edge['perm']):
            target=plans[physical_order[canonical_group]];plan=fleet['group_plans'][native]
            for r in plan['routes']:
                target['routes'].append(replace(r,job=ids[r.job],launch=r.launch+offset,
                    recovery=r.recovery+offset,sample=r.sample+offset))
            target['trail'].extend((p,t+offset,q,u+offset) for p,t,q,u in plan['trail'])
            target['assays'].extend((ids[j],s+offset,e+offset,k+offset) for j,s,e,k in plan['assays'])
        physical_order=[physical_order[i] for i in edge['reorder']];offset+=local.H
    check=validate_fleet(inst,plans)
    if not check['valid']:raise ValueError(('concatenated fleet',check['errors']))
    return dict(group_plans=plans,check=check)


def solve(inst,boundaries,seconds=60,*,dual_transfer=True,pooling=True,policy='J',checkpoint=None):
    began=time.perf_counter();deadline=began+seconds;objects={};events=[]
    report=dict(method='SK' if dual_transfer else 'SK-no-dual',status='unresolved',objective=None,
                lower_bound=0,optimal=False,fleet=None,proof_scope=None)
    cert=certify(inst,boundaries);report['separators']=cert
    # Merge uncertified cuts; this can reduce to one full-horizon exact kernel.
    boundaries=(0,)+tuple(c['epoch'] for c in cert['cuts'] if c['valid'])+(inst.H,)
    report['effective_boundaries']=boundaries
    stages=[];lb=0;ub=math.inf;upper_path=();iterations=0
    try:
        for local,ids in split(inst,boundaries):
            key=signature(local)
            if key not in objects:objects[key]=Kernel(local,deadline-1,policy,dual_transfer)
            stages.append((local,ids,key))
        while time.perf_counter()<deadline-1:
            iterations+=1
            lb,path,size=master(inst,stages,objects,pooling=pooling,deadline=deadline-1)
            report.update(lower_bound=None if not math.isfinite(lb) else lb,interface=size,iterations=iterations)
            candidate,cp,_=master(inst,stages,objects,pooling=pooling,upper=True,deadline=deadline-1)
            if candidate<ub:ub=candidate;upper_path=cp
            if math.isfinite(ub) and lb>=ub:break
            if not path:
                if math.isinf(lb):
                    report.update(status='infeasible',infeasible_proved=True,
                        proof_scope='complete_interface_exclusion_by_native_numerical_kernels')
                break
            unknown=next((e for e in path if e['context'] not in objects[e['kernel']].solved or
                not objects[e['kernel']].solved[e['context']]['optimal']),None)
            if unknown is None:ub=lb;upper_path=path;break
            k=objects[unknown['kernel']];r=k.solve(*unknown['context'],deadline-1)
            events.append(dict(kernel=unknown['kernel'],**r))
            if checkpoint:checkpoint('kernel',dict(event=plain(events[-1]),master=report))
            if r['status']=='unresolved':break
        if math.isfinite(ub):
            fleet=concatenate(inst,stages,objects,upper_path)
            if fleet['check']['objective']!=ub:raise ValueError('composed objective')
            closed=lb>=ub
            report.update(status='optimal' if closed else 'feasible',objective=ub,fleet=fleet,
                          optimal=closed,proof_scope='exact_interface_and_native_numerical_kernels' if closed else None)
    except TimeoutError as e:report['reason']=str(e)
    except Exception as e:report.update(status='correctness_error',reason=repr(e))
    report.update(seconds=time.perf_counter()-began,native_calls=sum(k.native_calls for k in objects.values()),
        kernel_build_seconds=sum(k.build_seconds for k in objects.values()),kernel_count=len(objects),
        evaluated_contexts=sum(len(k.solved) for k in objects.values()),
        dual_price_count=sum(len(k.prices) for k in objects.values()),events=events,
        kernel_shapes={key:list(k.m.A.shape) for key,k in objects.items()},
        interface_path=upper_path if math.isfinite(ub) else ())
    if report['seconds']>seconds:report.update(optimal=False,proof_scope=None,status='workflow_time_limit')
    return report
