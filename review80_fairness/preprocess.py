"""Common necessary vehicle reachability before either model construction."""
from collections import Counter
from contextlib import contextmanager
import hashlib
import json
import time


def shortest_road(inst):
    n=len(inst.pads)
    d=[[0 if p==q else inst.road(p,q) for q in range(n)] for p in range(n)]
    for k in range(n):
        for p in range(n):
            for q in range(n):d[p][q]=min(d[p][q],d[p][k]+d[k][q])
    return d


def rejection(inst,r,dist):
    if dist[r.p][r.q]>r.recovery-r.launch:return 'launch_return'
    if dist[0][r.p]>r.launch:return 'depot_launch'
    if r.recovery+dist[r.q][0]>inst.H:return 'return_depot'
    return None


def filter_routes(inst,routes,deadline=None):
    began=time.perf_counter();dist=shortest_road(inst);kept=[];excluded=Counter();contexts=Counter()
    before_cross=after_cross=0
    for i,r in enumerate(routes):
        if i%128==0 and deadline is not None and time.perf_counter()>=deadline:
            raise TimeoutError('common exact vehicle reachability budget')
        before_cross+=r.p!=r.q
        reason=rejection(inst,r,dist)
        if reason:
            excluded[reason]+=1;contexts[(r.p,r.launch,r.q,r.recovery,reason)]+=1
        else:
            kept.append(r);after_cross+=r.p!=r.q
    return tuple(kept),dict(before=len(routes),after=len(kept),
        before_cross_pad=before_cross,after_cross_pad=after_cross,excluded=dict(excluded),
        rejected_contexts=[dict(p=p,launch=l,q=q,recovery=e,reason=c,sorties=n)
                           for (p,l,q,e,c),n in sorted(contexts.items())],
        shortest_road=dist,groups_in_call=inst.groups,filter_seconds=time.perf_counter()-began)


@contextmanager
def common_enumeration(events):
    from fleet_sample_assay_queue_branch import pattern
    from fleet_sample_assay_group import pricing
    originals=[]
    for module in (pattern,pricing):
        old=module.enumerate_routes;originals.append((module,old))
        def wrapped(inst,deadline,old=old):
            began=time.perf_counter();routes=old(inst,deadline)
            kept,record=filter_routes(inst,routes,deadline)
            record.update(total_enumeration_and_filter_seconds=time.perf_counter()-began)
            events.append(record)
            return kept
        module.enumerate_routes=wrapped
    try:yield
    finally:
        for module,old in originals:module.enumerate_routes=old
