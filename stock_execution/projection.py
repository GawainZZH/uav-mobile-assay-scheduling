"""Exact checks supporting the compact-bound path-lift argument."""
from fractions import Fraction as F
from dataclasses import replace
from collections import defaultdict
from raw_stock_flow.initial import family_lower
from service_separator25_joint.solve import signature


def certify(a):
    checked=excluded=0
    for f in a.domain.families:
        local,ids=a.domain.stages[f['stage']];t=a.templates[signature(replace(local,groups=1))]
        S=tuple(f['local_subset']);p=f['start'];q=f['end']
        impossible=any(t['single'][j,p,q]['status']=='infeasible'for j in S)
        if impossible:
            if not f['infeasible']:raise ValueError('missing singleton exclusion in flow family')
            excluded+=1;continue
        if f['infeasible']:excluded+=1;continue
        lower=family_lower(t['elementary'],{j:t['single'][j,p,q]['objective']for j in S},S,
            local.road_cost*t['roads']['distance'][p][q])
        if F(f['lower'])<lower:raise ValueError('flow lower fails compact column inequality')
        checked+=1
    return dict(valid=True,families_checked=checked,excluded_families=excluded,
        scope='compact analytic lower and singleton exclusions only',new_optimizer_calls=0)


def lift(m,vector):
    """A common per-group fractional path mixture, not a physical calendar.

    Exact matrix feasibility includes all retained inequalities. Distribute
    each nonempty-family cost in proportion to its arc flow; empty cost is
    fixed. Decompose the G-unit stock flow, then use weight/G for every group.
    """
    x=[F(str(v))for v in vector];d=m.domain;G=d.groups;K=len(d.tasks)
    if len(x)!=len(m.cost)or any(v<0 or v>u for v,u in zip(x,m.upper)):raise ValueError('flow vector domain')
    for r in m.rows:
        v=sum((c*x[j]for j,c in r['entries'].items()),F(0))-r['rhs']
        if v!=0 if r['eq']else v>0:raise ValueError('flow row violation')
    weight={i:sum((x[j]for j in m.family_columns[i]),F(0))for i in range(len(m.bands))}
    unit={i:(x[m.z[i]]/weight[i]if weight[i]else F(0))if i in m.z else F(f['lower'])
        for i,f in enumerate(m.bands)}
    left=[x[j]for j in m.flow_columns];adj=defaultdict(list)
    for ai,arc in enumerate(d.arcs):adj[arc['tail']].append(ai)
    source=(0,0,d.budget);paths=[]
    while any(left[i]>0 for i in adj[source]):
        node=source;path=[]
        for k in range(K):
            choices=[i for i in adj[node]if left[i]>0]
            if not choices:raise ValueError('nonconserving positive stock flow')
            i=choices[0];path.append(i);node=d.arcs[i]['head']
        w=min(left[i]for i in path)
        for i in path:left[i]-=w
        paths.append(dict(weight=w,per_group_weight=w/G,arcs=path,
            unit_stage_costs=[unit[d.arcs[i]['family']]for i in path]))
    if any(left)or sum((p['weight']for p in paths),F(0))!=G:raise ValueError('stock path decomposition incomplete')
    stages=[]
    for k,ids in enumerate(d.tasks):
        fs=[i for i,f in enumerate(m.bands)if f['stage']==k]
        own={j:sum((weight[i]for i in fs if j in m.bands[i]['subset']),F(0))/G for j in ids}
        starts={p:sum((weight[i]for i in fs if m.bands[i]['start']==p),F(0))/G for p in range(d.pads)}
        ends={p:sum((weight[i]for i in fs if m.bands[i]['end']==p),F(0))/G for p in range(d.pads)}
        cost=sum((weight[i]*unit[i]for i in fs),F(0))/G
        if cost!=x[m.theta[k]]/G:raise ValueError('cost disaggregation differs')
        stages.append(dict(ownership=own,start=starts,end=ends,cost=cost))
    return dict(group_count=G,paths=paths,stages=stages,scope='fractional analytic feature lift; not original physical feasibility')
