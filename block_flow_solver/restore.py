"""Stock paths to a physically validated original fleet."""
from dataclasses import replace
from fractions import Fraction as F
from block_group_flow.contracts import decompose
from review88_physical.types import Route
from review88_physical.physical import validate_fleet
from .bands import exact


def reconstruct(m,x):
    d=m.domain
    if d.inst is None:raise ValueError('physical reconstruction requires an original instance')
    flows=[];arcs=[]
    for a,j in zip(d.arcs,m.flow_columns):
        v=round(float(x[j]))
        if abs(v-float(x[j]))>1e-6:raise ValueError('fractional stock path')
        i=a['family'];b=m.bands[i]
        if v and not exact(b):raise ValueError('selected family lacks exact feasible local plan')
        flows.append(v);arcs.append(dict(a,cost=b['upper']))
    proof=decompose(arcs,flows,groups=d.groups,pads=d.pads,budget=d.budget,tasks=d.tasks)
    plans=[]
    for path in proof['paths']:
        out=dict(routes=[],trail=[],assays=[])
        for k,ai in enumerate(path):
            local,ids=d.stages[k];offset=d.cuts[k];plan=m.bands[d.arcs[ai]['family']]['plan']
            for raw in plan['routes']:
                r=Route(**raw)if isinstance(raw,dict)else raw
                out['routes'].append(replace(r,job=ids[r.job],launch=r.launch+offset,recovery=r.recovery+offset,sample=r.sample+offset))
            out['trail'].extend((p,t+offset,q,u+offset)for p,t,q,u in plan['trail'])
            out['assays'].extend((ids[j],s+offset,e+offset,a+offset)for j,s,e,a in plan['assays'])
        plans.append(out)
    physical=validate_fleet(d.inst,plans)
    if not physical['valid']or F(str(physical['objective']))!=proof['objective']:
        raise ValueError('physical fleet invalid or costs differ from stock paths')
    return dict(objective=proof['objective'],physical=physical,proof=proof,fleet=dict(group_plans=plans,check=physical))
