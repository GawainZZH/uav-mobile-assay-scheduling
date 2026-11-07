"""One physical G1 matrix; original task labels stay fixed across cover masks."""
from dataclasses import replace
from fractions import Fraction as F
from math import lcm
from types import SimpleNamespace
import numpy as np

from review88_physical.types import Route
from service_separator25_joint.model import build,complete
from service_separator25_joint.kernel import validate_boundary
from revision37.carrier_normalization import normalize
from singleton_context_oracle import road_data


def checked_subset(inst,subset):
    subset=tuple(subset)
    if subset!=tuple(sorted(set(subset)))or any(type(j)is not int or not 0<=j<len(inst.jobs)for j in subset):
        raise ValueError('canonical original task subset required')
    if len(subset)>inst.batteries:raise ValueError('subset exceeds battery inventory')
    return subset


def subset_instance(inst,subset):
    subset=checked_subset(inst,subset);P=len(inst.pads);nodes=tuple(range(P))+tuple(P+j for j in subset)
    ft=None if inst.flight_times is None else tuple(tuple(inst.flight_times[a][b]for b in nodes)for a in nodes)
    return replace(inst,groups=1,jobs=tuple(inst.jobs[j]for j in subset),flight_times=ft),subset


def build_template(inst,deadline,*,crew_capacity=None):
    capacity=inst.drones if crew_capacity is None else crew_capacity
    if type(capacity)is not int or capacity<1:raise ValueError('positive original crew capacity required')
    m=build(replace(inst,groups=1),'J',deadline);m.M.eq[m.M.index['pack',0]]=True
    if capacity<inst.drones:
        for s in range(inst.H):m.M.row(('ground_crew',0,s),capacity)
        for i,r in enumerate(m.routes):
            for s in range(r.recovery,r.recovery+inst.ground_duration):m.M.columns[m.rids[0,i]][('ground_crew',0,s)]=1
        m.A=m.M.sparse().tocsc()
    cost=tuple(F(c)for c in m.M.cost);scale=lcm(*(v.denominator for v in cost))
    e=[]
    for j,job in enumerate(inst.jobs):
        values=[inst.flight_fixed+inst.air_cost*(r.recovery-r.launch)+job.weight*(r.recovery+inst.handoff_delay+inst.assay_duration-r.sample)
                for r in m.routes if r.job==j]
        e.append(min(values)if values else None)
    return SimpleNamespace(inst=inst,m=m,A=m.A,crew_capacity=capacity,eq=np.asarray(m.M.eq,dtype=bool),
        base=np.asarray(m.M.rhs,dtype=float),upper=np.asarray(m.M.upper,dtype=float),
        integer=np.asarray(m.M.integer,dtype=bool),cost_scale=scale,
        scaled_cost=np.asarray([int(c*scale)for c in cost],dtype=float),elementary_tasks=tuple(e),roads=road_data(inst))


def context_rhs(t,subset,start,end):
    subset=checked_subset(t.inst,subset)
    if any(type(p)is not int or not 0<=p<len(t.inst.pads)for p in (start,end)):raise ValueError('invalid pad')
    rhs=t.base.copy();selected=set(subset);M=t.m.M
    for j in range(len(t.inst.jobs)):rhs[M.index['cover',j]]=int(j in selected)
    rhs[M.index['pack',0]]=len(subset)
    for p in range(len(t.inst.pads)):
        rhs[M.index['flow',0,p,0]]=int(p==start)
        rhs[M.index['flow',0,p,t.inst.H]]=-int(p==end)
    return rhs


def elementary(t,subset,start,end):
    subset=checked_subset(t.inst,subset)
    if any(t.elementary_tasks[j]is None for j in subset):return None
    return sum((F(t.elementary_tasks[j])for j in subset),F(0))+t.inst.road_cost*t.roads['distance'][start][end]


def validate_subset(inst,plan,subset,start,end,*,crew_capacity=None):
    local,ids=subset_instance(inst,subset);mapping={j:i for i,j in enumerate(ids)}
    routes=[Route(**r)if isinstance(r,dict)else r for r in plan['routes']]
    if sorted(r.job for r in routes)!=list(ids):raise ValueError('physical subset cover differs')
    translated=dict(trail=plan['trail'],routes=[replace(r,job=mapping[r.job])for r in routes],
        assays=[(mapping[j],s,e,k)for j,s,e,k in plan['assays']])
    checked=validate_boundary(local,[translated],(start,),(end,))
    capacity=inst.drones if crew_capacity is None else crew_capacity
    if type(capacity)is not int or capacity<1:raise ValueError('positive original crew capacity required')
    if any(sum(r.recovery<=s<r.recovery+inst.ground_duration for r in routes)>capacity for s in range(inst.H)):
        checked.update(valid=False,errors=list(checked['errors'])+['original ground crew capacity exceeded'])
    return checked


def embed_plan(t,plan,subset,start,end):
    checked=validate_subset(t.inst,plan,subset,start,end,crew_capacity=t.crew_capacity)
    if not checked['valid']:raise ValueError(('invalid subset source',checked))
    m=t.m;x=np.zeros(len(t.upper))
    def key(r):return (r.p,r.launch,r.q,r.recovery,r.job,r.sample)
    rindex={key(r):i for i,r in enumerate(m.routes)};aindex={a:i for i,a in enumerate(m.arcs)}
    for raw in plan['routes']:
        r=Route(**raw)if isinstance(raw,dict)else raw
        x[m.rids[0,rindex[key(r)]]]=1
    for a in plan['trail']:x[m.zids[0,aindex[tuple(a)]]]=1
    for j,s,e,k in plan['assays']:x[m.assignments[0,j,s,e,k]]=1
    residual(t,x,context_rhs(t,subset,start,end),completed=True)
    if abs(np.asarray(m.M.cost)@x-float(checked['objective']))>1e-6:raise ValueError('source objective mismatch')
    return x,checked


def residual(t,vector,rhs,*,completed=False,relax=False):
    x=np.asarray(vector,dtype=float);a=t.A@x
    if x.shape!=t.upper.shape or not np.isfinite(x).all():raise ValueError('invalid subset vector')
    integer=np.zeros(len(x),dtype=bool)if relax else np.ones(len(x),dtype=bool)if completed else t.integer
    error=max(np.max(abs(a[t.eq]-rhs[t.eq]),initial=0),np.max(a[~t.eq]-rhs[~t.eq],initial=0),
        np.max(-x,initial=0),np.max(x-t.upper,initial=0),np.max(abs(x[integer]-np.rint(x[integer])),initial=0))
    if error>1e-6:raise ValueError(('actual subset RHS residual',float(error)))
    return float(error)


def restore_subset(t,raw,subset,start,end):
    rhs=context_rhs(t,subset,start,end);m=t.m;saved=list(m.M.rhs)
    try:
        m.M.rhs[:]=rhs.tolist();x,note=normalize(m,raw,rhs)
        x=complete(m,x,(start,),(end,))
    finally:m.M.rhs[:]=saved
    residual(t,x,rhs,completed=True)
    selected=[r for i,r in enumerate(m.routes)if x[m.rids[0,i]]>.5]
    arcs=[a for i,a in enumerate(m.arcs)if x[m.zids[0,i]]>.5]
    trail=[];node=(start,0)
    while node!=(end,t.inst.H):
        nxt=[a for a in arcs if a[:2]==node]
        if len(nxt)!=1:raise ValueError('subset carrier path')
        a=nxt[0];trail.append(a);node=a[2:]
        if len(trail)>len(arcs):raise ValueError('subset path cycle')
    if len(trail)!=len(arcs):raise ValueError('unused subset carrier arcs')
    plan=dict(routes=selected,trail=trail,assays=[(j,s,e,k)for(g,j,s,e,k),col in m.assignments.items()if x[col]>.5])
    checked=validate_subset(t.inst,plan,subset,start,end,crew_capacity=t.crew_capacity)
    if not checked['valid']or abs(float(checked['objective'])-np.asarray(m.M.cost)@x)>1e-6:
        raise ValueError(('subset physical objective mismatch',checked))
    return x,plan,checked,note
