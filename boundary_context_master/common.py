"""Common strengthened physical model, boundary coordinates and source lifts."""
from fractions import Fraction as F
from types import SimpleNamespace
import math
import time
import numpy as np

from review88_physical import milp as base
from service_separator25_joint.cases import split
from service_separator25_joint.interface import pack_limits
from service_separator25_joint.solve import signature
from referee_boundary21_v2.interface import certify
from revision22.strengthening import boundary_rule
from service_separator25_joint.monolithic import strengthening
from primal_transfer.incumbents import project_fleet
from .matrix import build,put_row,encode,decode,check_vector


def describe(inst,cuts,*,objects=None,deadline=math.inf):
    cuts=tuple(cuts);cert=certify(inst,cuts)
    if (0,)+tuple(c['epoch']for c in cert['cuts']if c['valid'])+(inst.H,)!=cuts:
        raise ValueError('all supplied boundaries must qualify')
    stages=[];floors=[];elementary=[];cache={}
    for local,ids in split(inst,cuts):
        key=signature(local);stages.append((local,ids,key))
        if key not in cache:
            if objects is not None:
                k=objects[key];cache[key]=(F(k.envelope['bound']),F(k.elementary))
            else:
                routes=base.independent_routes(local,deadline)
                value=sum(min(local.flight_fixed+local.air_cost*(r.recovery-r.launch)+
                    local.jobs[j].weight*(r.recovery+local.handoff_delay+local.assay_duration-r.sample)
                    for r in routes if r.job==j)for j in range(len(local.jobs)))
                cache[key]=(F(value),F(value))
        f,e=cache[key];floors.append(f);elementary.append(e)
    return SimpleNamespace(inst=inst,cuts=cuts,stages=tuple(stages),
        groups=inst.groups,pads=len(inst.pads),tasks=tuple(len(ids)for local,ids,key in stages),
        limits=tuple(pack_limits(local)for local,ids,key in stages),
        templates=tuple(key for local,ids,key in stages),floors=tuple(floors),elementary=tuple(elementary),
        road=tuple(tuple(inst.road_cost*inst.road(p,q)for q in range(len(inst.pads)))for p in range(len(inst.pads))),
        batteries=(inst.batteries,)*inst.groups)


def compact(d,upper,*,records=(),anchor=0,deadline=math.inf):
    return build(d.groups,d.pads,d.tasks,d.limits,d.batteries,d.floors,d.elementary,
                 d.road,d.templates,upper,records=records,anchor=anchor,deadline=deadline)


def stage_costs(d,m):
    result=[]
    for (local,ids,key),a,b in zip(d.stages,d.cuts,d.cuts[1:]):
        selected=set(ids);cost={}
        for (g,i),col in m.rids.items():
            if m.routes[i].job in selected and m.M.cost[col]:cost[col]=m.M.cost[col]
        for (g,j,s,e,k),col in m.assignments.items():
            if j in selected and m.M.cost[col]:cost[col]=m.M.cost[col]
        for (g,i),col in m.zids.items():
            p,t,q,u=m.arcs[i]
            if a<=t<u<=b and m.M.cost[col]:cost[col]=m.M.cost[col]
        result.append(cost)
    return tuple(result)


def build_base(d,deadline=math.inf):
    metadata={}
    with boundary_rule(True),strengthening(d.inst,d.cuts,metadata):
        m=base.build(d.inst,'J',deadline)
    if len(d.inst.jobs)==d.groups*d.inst.batteries:
        for g in range(d.groups):m.M.eq[m.M.index['pack',g]]=True
    for k,cost in enumerate(stage_costs(d,m)):
        put_row(m.M,('shared_stage_cost_floor',k),{j:-v for j,v in cost.items()},-d.floors[k])
    m.A=m.M.sparse().tocsc();m.common_strengthening=metadata
    return m


def augment(d,m,upper,*,records=(),anchor=0,deadline=math.inf):
    old_rows,old_columns=m.A.shape;costs=stage_costs(d,m)
    a=build(d.groups,d.pads,d.tasks,d.limits,d.batteries,d.floors,d.elementary,d.road,d.templates,
        upper,records=records,anchor=anchor,matrix=m.M,theta_objective=False,deadline=deadline)
    for k,ts in enumerate(a.types):
        for g in range(d.groups):
            for p in range(d.pads):
                for kind,c,side in [('start',d.cuts[k],0),('end',d.cuts[k+1],1)]:
                    entries={a.x[k,g,t]:1 for t in ts if t[side]==p}
                    for (gg,i),col in m.zids.items():
                        pp,t,q,u=m.arcs[i]
                        if gg==g and (pp==p and t==c if c<d.inst.H else q==p and u==c):
                            entries[col]=entries.get(col,0)-1
                    put_row(m.M,('cb','physical_pad',k,g,p,kind),entries,0,True)
            entries={a.x[k,g,t]:t[2]for t in ts}
            ids=set(d.stages[k][1])
            for (gg,i),col in m.rids.items():
                if gg==g and m.routes[i].job in ids:entries[col]=-1
            put_row(m.M,('cb','physical_pack',k,g),entries,0,True)
        entries={j:-v for j,v in costs[k].items()};entries[a.theta[k]]=1
        put_row(m.M,('cb','physical_cost',k),entries,0,True)
    a.A=m.M.sparse().tocsr();m.A=a.A.tocsc();a.physical=m;a.stage_costs=costs
    a.original_rows=old_rows;a.original_columns=old_columns
    return a


def project_actions(d,fleet):
    fragments=project_fleet(d.inst,d.cuts,fleet);actions=[]
    for (local,ids,key),a,b in zip(d.stages,d.cuts,d.cuts[1:]):
        stage=[]
        for plan in fleet['group_plans']:
            nodes={t:p for p,t,q,u in plan['trail']}
            nodes.update({u:q for p,t,q,u in plan['trail']})
            used=sum((r['job']if isinstance(r,dict)else r.job)in ids for r in plan['routes'])
            stage.append((nodes[a],nodes[b],used))
        actions.append(tuple(stage))
    return tuple(actions),fragments


def embed_source(d,a,fleet):
    actions,fragments=project_actions(d,fleet)
    vector=encode(a,actions,tuple(f['objective']for f in fragments));m=a.physical
    def route_key(r):
        rr=r if isinstance(r,dict)else vars(r)
        return tuple(rr[k]for k in ('p','launch','q','recovery','job','sample'))
    routes={route_key(r):i for i,r in enumerate(m.routes)};arcs={tuple(r):i for i,r in enumerate(m.arcs)}
    for g,plan in enumerate(fleet['group_plans']):
        for r in plan['routes']:vector[m.rids[g,routes[route_key(r)]]]=1
        for arc in plan['trail']:vector[m.zids[g,arcs[tuple(arc)]]]=1
        for j,s,e,k in plan['assays']:vector[m.assignments[g,j,s,e,k]]=1
    check_vector(a,vector)
    if abs(float(np.asarray(m.M.cost)@vector)-sum(float(f['objective'])for f in fragments))>1e-6:
        raise ValueError('source physical/stage cost mismatch')
    return vector,fragments,actions


def complete_and_restore(d,a,raw):
    from service_separator25_joint.model import complete
    from revision37.carrier_normalization import normalize
    decoded=decode(a,raw);m=a.physical;saved=list(m.M.upper)
    try:
        for k,c in enumerate(d.cuts[1:-1]):
            for (g,i),col in m.zids.items():
                p,t,q,u=m.arcs[i]
                if t==c and p!=decoded['actions'][k][g][1]:m.M.upper[col]=0
        vector,_=normalize(m,raw)
        vector=complete(m,vector,(0,)*d.groups,(0,)*d.groups)
        for k,cost in enumerate(a.stage_costs):
            vector[a.theta[k]]=sum(float(v)*vector[j]for j,v in cost.items())
    finally:
        m.M.upper[:]=saved
    check_vector(a,vector);fleet,physical=base.restore(m,vector)
    if abs(float(np.asarray(m.M.cost)@vector)-float(physical['objective']))>1e-6:
        raise ValueError('completed physical objective differs')
    return vector,fleet
