"""Complete old strong rows, fewer scheduling integers and explicit ownership."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
import numpy as np
from vehicle_root_ablation.model import append_rows,adapter,assert_inventory


def extend(base,physical_columns,ownership):
    if not 0<physical_columns<=base['shape'][1]:raise ValueError('physical column inventory')
    seen=set()
    for cols in ownership.values():
        if any(not 0<=j<physical_columns or j in seen for j in cols):raise ValueError('ownership supports must be disjoint physical events')
        seen.update(cols)
    out={k:list(v)if isinstance(v,list)else v for k,v in base.items()};out['shape']=list(base['shape']);out['x']=None
    out['integer'][:physical_columns]=[0]*physical_columns;s={}
    for key,cols in sorted(ownership.items()):
        j=out['shape'][1];s[key]=j;out['shape'][1]+=1
        out['cost'].append(0);out['upper'].append(1);out['integer'].append(1)
        i=out['shape'][0];out['shape'][0]+=1;out['rhs'].append(0);out['equalities'].append(True)
        for col,value in [(j,1)]+[(c,-1)for c in cols]:
            out['row'].append(i);out['col'].append(col);out['data'].append(value)
    return out,s


def lift(vector,positions,ownership):
    out=list(vector)+[0]*len(positions)
    for key,col in positions.items():out[col]=sum(vector[j]for j in ownership[key])
    return out


def check_costs(base,costs,theta,cost_rows):
    total=[F(0)]*base['shape'][1];expected={k:{theta[k]:F(1)}for k in theta}
    for (k,g),row in costs.items():
        for j,v in row.items():
            v=F(str(v))
            if v<0:raise ValueError('negative costs not covered by conditional-row proof')
            total[j]+=v;expected[k][j]=expected[k].get(j,F(0))-v
    actual={k:{}for k in theta};reverse={row:k for k,row in cost_rows.items()}
    for i,j,v in zip(base['row'],base['col'],base['data']):
        if i in reverse:
            row=actual[reverse[i]];row[j]=row.get(j,F(0))+F(str(v))
    for k,i in cost_rows.items():
        if not base['equalities'][i]or base['rhs'][i]!=0 or actual[k]!=expected[k]:raise ValueError('group/stage theta identity differs')
    for j,c in enumerate(base['cost']):
        if base['upper'][j]>0 and total[j]!=F(str(c)):raise ValueError('group/block global cost identity differs')


def prepare_problem(payload,record,deadline):
    from review88_physical.evidence import instance
    from review88_physical.milp import dump_matrix
    from boundary_context_master.common import describe,build_base,augment,stage_costs
    inst=instance(payload['instance']);cuts=tuple(payload['cuts']);common=record['common']
    d=describe(inst,cuts,deadline=deadline)
    d.floors=tuple(F(str(record['envelopes'][key]['bound']))for key in d.templates)
    m=build_base(d,deadline)
    a=augment(d,m,common['source_objective'],records=common['shared_cost_rows'],deadline=deadline)
    base=common['augmented_strong_LP']['matrix'];fresh=dump_matrix(a,np.asarray(a.M.cost),None)
    assert_inventory(base,fresh)
    ownership={};costs={};physical_group={}
    for mapping in (m.zids,m.rids,m.assignments):
        for key,col in mapping.items():physical_group[col]=key[0]
    if set(physical_group)!=set(range(common['original_shape'][1])):raise ValueError('unmapped physical columns')
    for k,(local,ids,template)in enumerate(d.stages):
        for g in range(inst.groups):
            for j,original_j in enumerate(ids):
                ownership[k,g,j]=tuple(col for (gg,i),col in m.rids.items()if gg==g and m.routes[i].job==original_j)
            costs[k,g]={col:F(str(v))for col,v in a.stage_costs[k].items()if physical_group[col]==g}
    check_costs(base,costs,a.theta,{k:a.M.index['cb','physical_cost',k]for k in a.theta})
    anchored=append_rows(base,[dict(name='common_anchor',entries={j:-1 for j in a.theta.values()},rhs=-common['common_integer_anchor'])])
    matrix,s=extend(anchored,common['original_shape'][1],ownership)
    source=lift(common['source_vector'],s,ownership)
    from vehicle_root_ablation.model import exact_vector
    if exact_vector(matrix,source,integer_flags=True)!=F(str(common['source_objective'])):raise ValueError('extended source mismatch')
    return NS(inst=inst,d=d,a=a,m=m,base=base,matrix=matrix,s=s,ownership=ownership,costs=costs,
              source_vector=source,source=payload['fleet'],anchor=common['common_integer_anchor'],
              adapter=adapter(matrix),physical_columns=common['original_shape'][1])


def decode(problem,vector):
    features=[]
    for k,(local,ids,key)in enumerate(problem.d.stages):
        for g in range(problem.inst.groups):
            selected=[t for t in problem.a.types[k]if vector[problem.a.x[k,g,t]]>.5]
            if len(selected)!=1:raise ValueError('unique integer endpoint/consumption choice required')
            p,q,v=selected[0]
            values=[vector[problem.s[k,g,j]]for j in range(len(ids))]
            if any(abs(x-round(x))>1e-6 for x in values):raise ValueError('fractional ownership at integer interface')
            S=tuple(j for j,x in enumerate(values)if x>.5)
            if len(S)!=v:raise ValueError('pack/ownership mismatch')
            cost=sum(float(c)*vector[j]for j,c in problem.costs[k,g].items())
            features.append(dict(stage=k,group=g,template=key,subset=S,start=p,end=q,used=v,cost=cost,
                                 choice=problem.a.x[k,g,(p,q,v)]))
    return features
