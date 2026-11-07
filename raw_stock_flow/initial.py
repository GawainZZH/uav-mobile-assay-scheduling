"""Physics-only template initialization and an optional constructive incumbent."""
import time
from math import comb
from dataclasses import replace
from fractions import Fraction as F
from collections import Counter
from types import SimpleNamespace as NS
from review88_physical.evidence import instance
from review88_physical.types import Route
from review88_physical.physical import validate_fleet
from service_separator25_joint.cases import split
from service_separator25_joint.solve import signature
from service_separator25_joint.interface import pack_limits
from referee_boundary21_v2.interface import certify
from block_flow_solver.domain import from_instance
from block_flow_solver.matrix import build,row
from group_subset_oracle.model import validate_subset
from singleton_context_oracle import road_data,_single
from .result import currency


class DomainLimit(ValueError):pass
class IndependentTaskInfeasible(ValueError):pass


def clock(deadline):
    if time.perf_counter()>=deadline:raise TimeoutError('raw initialization deadline')


def family_lower(elementary,singletons,selected,road):
    E=sum((F(elementary[j])for j in selected),F(0))
    extra=max([F(road)]+[F(singletons[j])-F(elementary[j])for j in selected])
    return E+extra


def source_assignment(tasks,*,groups,budget):
    if any(len(S)>groups for S in tasks)or sum(map(len,tasks))>groups*budget:return None
    loads=[0]*groups;assignment=[[]for _ in range(groups)]
    for S in tasks:
        used=set()
        for j in S:
            choices=[g for g in range(groups)if g not in used and loads[g]<budget]
            if not choices:return None
            g=min(choices,key=lambda g:(loads[g],g));used.add(g);loads[g]+=1;assignment[g].append(j)
    return tuple(tuple(S)for S in assignment)


def initialize(raw,deadline,*,max_families=20000,max_arcs=200000):
    began=time.perf_counter()
    if set(raw)!={'instance','cuts'}:raise ValueError('raw instance/cuts only; no old answers, models, bounds or keys')
    original=instance(raw['instance']);inst,D,quantum=currency(original);cuts=tuple(raw['cuts'])
    if not certify(inst,cuts)['valid']:raise DomainLimit('supplied cuts do not qualify; merge before decomposition')
    stages=split(inst,cuts);P=len(inst.pads);K=len(stages);nf=na=0
    for k,(local,ids)in enumerate(stages):
        lo,hi=pack_limits(local);ends=(1 if k==0 else P)*(1 if k==K-1 else P)
        for size in range(lo,min(hi,inst.batteries,len(ids))+1):
            n=comb(len(ids),size)*ends;nf+=n;na+=n*(1 if k==0 else inst.batteries-size+1)
    if nf>max_families or na>max_arcs:raise DomainLimit(('stock family/arc limit',nf,na))
    clock(deadline);d=from_instance(inst,cuts);templates={};records=[];floors=[]
    from review88_physical.milp import independent_routes
    for k,(local,ids)in enumerate(d.stages):
        clock(deadline);key=signature(replace(local,groups=1))
        if key not in templates:
            routes=independent_routes(local,deadline);byjob=[[]for _ in ids]
            for r in routes:byjob[r.job].append(r)
            elementary=[]
            for j in range(len(ids)):
                if not byjob[j]:raise IndependentTaskInfeasible(('no physical sortie for required job',ids[j]))
                elementary.append(min(local.flight_fixed+local.air_cost*(r.recovery-r.launch)+
                    local.jobs[j].weight*(r.recovery+local.handoff_delay+local.assay_duration-r.sample)for r in byjob[j]))
            roads=road_data(local);single={}
            for j in range(len(ids)):
                for a in range(P):
                    for b in range(P):
                        clock(deadline);r=_single(local,j,a,b,byjob[j],roads,deadline)
                        if r['status']=='optimal':
                            checked=validate_subset(local,r['plan'],(j,),a,b,crew_capacity=inst.drones)
                            if not checked['valid']or checked['objective']!=r['objective']:raise ValueError('analytic singleton physical mismatch')
                        single[j,a,b]=r
            templates[key]=dict(elementary=tuple(elementary),single=single,roads=roads,route_count=len(routes))
        t=templates[key];floors.append(sum(t['elementary'],F(0)))
        for i,f in enumerate(d.families):
            if f['stage']!=k or not f['subset']:continue
            S=tuple(f['local_subset']);a=f['start'];b=f['end'];values=[t['single'][j,a,b]for j in S]
            if any(r['status']=='infeasible'for r in values):
                f.update(lower=F(0),upper=F(0),infeasible=True,exact=False,feasible=False,plan=None,proof_scope='removal_to_infeasible_singleton')
            elif len(S)==1:
                r=values[0];f.update(lower=r['objective'],upper=r['objective'],plan=r['plan'],feasible=True,exact=True,
                    proof_scope='analytic_singleton_complete_enumeration')
            else:
                singleton={j:t['single'][j,a,b]['objective']for j in S}
                f['lower']=family_lower(t['elementary'],singleton,S,local.road_cost*t['roads']['distance'][a][b])
                f['proof_scope']='task_intrinsic_costs_plus_one_shared_road_envelope'
            records.append(dict(family=i,stage=k,subset=f['subset'],start=a,end=b,lower=f['lower'],upper=f['upper'],infeasible=f['infeasible'],proof_scope=f['proof_scope']))
    assignment=source_assignment(d.tasks,groups=d.groups,budget=d.budget);source=None;paths=None
    index={(f['stage'],f['subset'],f['start'],f['end']):i for i,f in enumerate(d.families)}
    if assignment is not None:
        paths=[]
        for jobs in assignment:
            path=[]
            for k,S in enumerate(d.tasks):
                selected=tuple(j for j in S if j in jobs);i=index.get((k,selected,0,0))
                if i is None or not d.families[i]['exact']:path=None;break
                path.append(i)
            if path is None:paths=None;break
            paths.append(path)
        if paths is not None:
            plans=[]
            for path in paths:
                plan=dict(routes=[],trail=[],assays=[])
                for k,i in enumerate(path):
                    local,ids=d.stages[k];offset=cuts[k];p=d.families[i]['plan']
                    for rawroute in p['routes']:
                        r=Route(**rawroute)if isinstance(rawroute,dict)else rawroute
                        plan['routes'].append(replace(r,job=ids[r.job],launch=r.launch+offset,recovery=r.recovery+offset,sample=r.sample+offset))
                    plan['trail'].extend((a,t+offset,b,u+offset)for a,t,b,u in p['trail'])
                    plan['assays'].extend((ids[j],s+offset,e+offset,a+offset)for j,s,e,a in p['assays'])
                plans.append(plan)
            checked=validate_fleet(inst,plans)
            if not checked['valid']:raise ValueError(('constructed raw source invalid',checked['errors']))
            source=dict(group_plans=plans,check=checked)
    return NS(domain=d,original=original,currency_scale=D,quantum=quantum,source=source,paths=paths,
        floors=tuple(floors),floor=sum(floors,F(0)),records=records,templates=templates,historical_keys=0,
        family_count=nf,arc_count=na,seconds=time.perf_counter()-began)


def master(initial):
    d=initial.domain;upper=None if initial.source is None else initial.source['check']['objective']
    m=build(d,cutoff=upper)
    for k,floor in enumerate(initial.floors):m.rows.append(row({m.theta[k]:-1},-floor,False,'analytic_stage_floor'))
    if initial.source is None:return m,None
    x=[F(0)]*len(m.cost);counts=Counter();costs=Counter();stagecost=Counter()
    arcs={(tuple(a['tail']),a['family']):i for i,a in enumerate(d.arcs)}
    for path in initial.paths:
        left=d.budget;pad=0
        for k,i in enumerate(path):
            b=m.bands[i];ai=arcs[(k,pad,left),i];x[m.flow_columns[ai]]+=1
            counts[k,(b['start'],b['end'],len(b['subset']))]+=1;costs[i]+=b['upper'];stagecost[k]+=b['upper']
            pad=b['end'];left-=len(b['subset'])
    for i,j in m.z.items():x[j]=costs[i]
    for k,j in m.theta.items():x[j]=stagecost[k]
    for(k,t,r),j in m.u.items():x[j]=int(counts[k,t]>=r)
    return m,x
