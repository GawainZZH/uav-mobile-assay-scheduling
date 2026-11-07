"""Shared physical candidate and polynomial baseline envelope; no optimizer."""
import time
from dataclasses import replace
from fractions import Fraction as F
from raw_stock_flow.initial import initialize,DomainLimit,family_lower,clock
from raw_stock_flow.result import currency
from review88_physical.evidence import instance,plain
from review88_physical.types import Route
from review88_physical.physical import validate_fleet
from service_separator25_joint.cases import split
from service_separator25_joint.solve import signature
from group_subset_oracle.model import validate_subset


def prepare(raw,witness,method,deadline,*,max_families=20000,max_arcs=200000):
    began=time.perf_counter()
    if method not in('F','S'):raise ValueError('fixed F/S methods')
    try:
        a=initialize(raw,deadline,max_families=max_families,max_arcs=max_arcs);lite=False
    except DomainLimit as exc:
        if method!='S' or not isinstance(exc.args[0],tuple) or exc.args[0][0]!='stock family/arc limit':raise
        original=instance(raw['instance']);work,D,quantum=currency(original)
        # split restores each local pack capacity to its task count. Only
        # global subset materialization is limited by this preprocessing view.
        small=dict(instance=plain(replace(original,batteries=1)),cuts=raw['cuts'])
        a=initialize(small,deadline);lite=True
        a.original=original;a.currency_scale=D;a.quantum=quantum
        a.domain.inst=work;a.domain.budget=work.batteries;a.domain.stages=split(work,a.domain.cuts)
        a.domain.arcs=[];a.arc_count=0;a.source=None;a.paths=None
        if any(signature(replace(local,groups=1))not in a.templates for local,ids in a.domain.stages):
            raise ValueError('singleton view changed local physics')
    # Both methods normalize the same supplied candidate. An incidental default
    # assignment from the full enumerator must not give F a different incumbent.
    a.source=None;a.paths=None
    a.compact_fallback=lite;d=a.domain
    checked=validate_fleet(a.original,witness['group_plans'])
    if not checked['valid']:raise ValueError('shared candidate is not physical')
    index={(f['stage'],tuple(f['subset']),f['start'],f['end']):i for i,f in enumerate(d.families)}
    paths=[]
    for plan in witness['group_plans']:
        path=[]
        for k,(local,ids)in enumerate(d.stages):
            clock(deadline);lo,hi=d.cuts[k:k+2];mapping={j:i for i,j in enumerate(ids)}
            trail=[tuple(v)for v in plan['trail']if lo<=v[1]<v[3]<=hi]
            if not trail or trail[0][1]!=lo or trail[-1][3]!=hi:raise ValueError('candidate must expose every qualified boundary')
            p,q=trail[0][0],trail[-1][2]
            routes=[]
            for rawroute in plan['routes']:
                r=Route(**rawroute)if isinstance(rawroute,dict)else rawroute
                if r.job in mapping:routes.append(replace(r,job=mapping[r.job],launch=r.launch-lo,recovery=r.recovery-lo,sample=r.sample-lo))
            selected=tuple(sorted(r.job for r in routes));globalS=tuple(ids[j]for j in selected)
            sub=dict(routes=routes,trail=[(p,t-lo,q,u-lo)for p,t,q,u in trail],
                assays=[(mapping[j],s-lo,e-lo,t-lo)for j,s,e,t in plan['assays']if j in mapping])
            check=validate_subset(local,sub,selected,p,q,crew_capacity=d.inst.drones)
            if not check['valid']:raise ValueError(('candidate block invalid',check))
            key=(k,globalS,p,q)
            if key not in index:
                if not lite:raise ValueError('physical candidate missing from full family domain')
                t=a.templates[signature(replace(local,groups=1))]
                singles={j:t['single'][j,p,q]['objective']for j in selected}
                lower=family_lower(t['elementary'],singles,selected,local.road_cost*t['roads']['distance'][p][q])
                index[key]=len(d.families);d.families.append(dict(stage=k,subset=globalS,local_subset=selected,start=p,end=q,
                    lower=lower,upper=check['objective'],feasible=False,exact=False,infeasible=False,plan=None))
            i=index[key];f=d.families[i]
            if f['infeasible']or check['objective']<f['lower']:raise ValueError('shared candidate contradicts valid lower bound')
            if not f['feasible']or check['objective']<f['upper']:
                f.update(upper=check['objective'],plan=sub,feasible=True,exact=f['lower']==check['objective'],
                    proof_scope='analytic_lower_and_verified_shared_candidate')
            path.append(i)
        paths.append(path)
    plans=[]
    for path in paths:
        plan=dict(routes=[],assays=[],trail=[])
        for k,i in enumerate(path):
            local,ids=d.stages[k];offset=d.cuts[k];p=d.families[i]['plan']
            for rr in p['routes']:
                r=Route(**rr)if isinstance(rr,dict)else rr
                plan['routes'].append(replace(r,job=ids[r.job],launch=r.launch+offset,recovery=r.recovery+offset,sample=r.sample+offset))
            plan['assays'].extend((ids[j],s+offset,e+offset,t+offset)for j,s,e,t in p['assays'])
            plan['trail'].extend((p,t+offset,q,u+offset)for p,t,q,u in p['trail'])
        plans.append(plan)
    check=validate_fleet(d.inst,plans)
    if not check['valid']:raise ValueError(('normalized candidate invalid',check))
    if a.source is None or check['objective']<a.source['check']['objective']:
        a.source=dict(group_plans=plans,check=check);a.paths=paths
    a.family_count=len(d.families);a.seconds=time.perf_counter()-began
    return a


def compact_rows(q,a):
    """E-C <= 0 plus endpoint/task activated singleton/road increments."""
    rows=[];K=len(a.domain.stages);P=len(q.inst.pads)
    for k,(local,ids)in enumerate(a.domain.stages):
        t=a.templates[signature(replace(local,groups=1))]
        for g in range(q.inst.groups):
            base={q.s[k,g,j]:F(v)for j,v in enumerate(t['elementary'])}
            for col,v in q.costs[k,g].items():base[col]=base.get(col,F(0))-v
            rows.append(dict(entries=base,rhs=0,kind='compact_intrinsic_floor'))
            for p in ((0,)if k==0 else range(P)):
                for z in ((0,)if k==K-1 else range(P)):
                    ends=(()if k==0 else(q.y[k,g,p],))+(()if k==K-1 else(q.y[k+1,g,z],))
                    def add(active,delta,kind):
                        if delta<0:raise ValueError('negative singleton increment')
                        if not delta:return
                        e=dict(base)
                        for col in active:e[col]=e.get(col,F(0))+delta
                        rows.append(dict(entries=e,rhs=delta*(len(active)-1),kind=kind))
                    add(ends,local.road_cost*t['roads']['distance'][p][z],'compact_shared_road')
                    for j in range(len(ids)):
                        r=t['single'][j,p,z];active=ends+(q.s[k,g,j],)
                        if r['status']=='infeasible':
                            rows.append(dict(entries={c:1 for c in active},rhs=len(active)-1,kind='compact_singleton_exclusion'))
                        else:add(active,F(r['objective'])-t['elementary'][j],'compact_singleton_increment')
    return rows
