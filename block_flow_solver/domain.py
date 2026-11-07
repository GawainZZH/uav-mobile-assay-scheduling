"""Finite task/endpoint families and consumable-stock arcs, without event LP."""
from fractions import Fraction as F
from itertools import combinations,product
from types import SimpleNamespace as NS


def finite(*,groups,pads,budget,tasks,families):
    if any(type(x)is not int or x<1 for x in (groups,pads))or type(budget)is not int or budget<0:
        raise ValueError('positive group/pad counts and nonnegative integer stock required')
    tasks=tuple(tuple(s)for s in tasks);K=len(tasks);flat=sum(tasks,())
    if not K or len(set(flat))!=len(flat)or any(type(j)is not int for j in flat):raise ValueError('task partition required')
    fs=[];arcs=[];seen=set()
    for raw in families:
        f=dict(raw);k=f['stage'];S=tuple(f['subset']);p=f['start'];q=f['end'];key=(k,S,p,q)
        if key in seen or type(k)is not int or not 0<=k<K or S!=tuple(sorted(set(S)))or not set(S)<=set(tasks[k]):
            raise ValueError('unique canonical stage family required')
        if any(type(x)is not int or not 0<=x<pads for x in (p,q))or k==0 and p!=0 or k==K-1 and q!=0:
            raise ValueError('family endpoint domain differs')
        if len(S)>budget:raise ValueError('family exceeds original stock')
        seen.add(key);f.update(subset=S,lower=F(str(f['lower'])),upper=F(str(f['upper'])))
        if not 0<=f['lower']<=f['upper']:raise ValueError('finite nonnegative cost interval required')
        for name,default in [('feasible',False),('exact',False),('plan',None),('infeasible',False)]:f.setdefault(name,default)
        i=len(fs);fs.append(f)
        for b in ((budget,)if k==0 else range(len(S),budget+1)):
            arcs.append(dict(family=i,tail=(k,p,b),head=(k+1,q,b-len(S)),subset=S))
    return NS(groups=groups,pads=pads,budget=budget,tasks=tasks,families=fs,arcs=arcs,stages=None,inst=None,cuts=None)


def from_instance(inst,cuts):
    from referee_boundary21_v2.interface import certify
    from service_separator25_joint.cases import split
    from service_separator25_joint.interface import pack_limits
    from singleton_context_oracle import road_data,_trail
    from group_subset_oracle.model import validate_subset
    cuts=tuple(cuts)
    if not certify(inst,cuts)['valid']:raise ValueError('all boundaries must be qualified')
    stages=split(inst,cuts);families=[];P=len(inst.pads)
    for k,(local,ids)in enumerate(stages):
        lo,hi=pack_limits(local);roads=road_data(local)
        for size in range(lo,min(hi,inst.batteries,len(ids))+1):
            for selected in combinations(range(len(ids)),size):
                S=tuple(ids[j]for j in selected)
                cap=local.road_cost*local.H+size*(local.flight_fixed+local.air_cost*local.endurance)+local.H*sum((local.jobs[j].weight for j in selected),F(0))
                for p,q in product((0,)if k==0 else range(P),(0,)if k==len(stages)-1 else range(P)):
                    f=dict(stage=k,subset=S,local_subset=selected,start=p,end=q,lower=0,upper=cap)
                    if not S:
                        if roads['distance'][p][q]>local.H:f.update(upper=0,infeasible=True)
                        else:
                            plan=dict(routes=[],assays=[],trail=_trail(local,roads,p,q))
                            v=validate_subset(local,plan,(),p,q,crew_capacity=inst.drones)
                            if not v['valid']:raise ValueError('shortest empty plan invalid')
                            f.update(lower=v['objective'],upper=v['objective'],plan=plan,feasible=True,exact=True)
                    families.append(f)
    d=finite(groups=inst.groups,pads=P,budget=inst.batteries,tasks=tuple(ids for _,ids in stages),families=families)
    d.stages=stages;d.inst=inst;d.cuts=cuts
    return d
