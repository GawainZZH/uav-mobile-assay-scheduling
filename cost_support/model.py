"""Exact continuous extension, retaining every old strengthened row."""
from copy import copy,deepcopy
from fractions import Fraction as F
from math import isfinite
import time
from vehicle_root_ablation.model import adapter


def numeric_adapter(matrix):
    # Keep exact coefficients in the saved representation; SciPy requires a
    # numeric execution view instead of an object array of Fractions.
    numeric=dict(matrix);numeric['data']=[float(F(str(c)))for c in matrix['data']]
    return adapter(numeric)


def _validated_costs(p):
    M=p.matrix;K=len(p.d.stages);G=p.inst.groups
    if set(p.costs)!={(k,g)for k in range(K)for g in range(G)}:raise ValueError('complete group-block cost inventory required')
    costs={};total={};expected={k:{p.a.theta[k]:F(1)}for k in range(K)}
    for key,raw in p.costs.items():
        k,g=key;row={}
        for j,value in raw.items():
            c=F(str(value))
            if type(j)is not int or not 0<=j<p.physical_columns or c<0:raise ValueError('nonnegative original physical group costs required')
            if c:
                row[j]=c;total[j]=total.get(j,F(0))+c;expected[k][j]=expected[k].get(j,F(0))-c
        costs[key]=row
    for j,c in enumerate(M['cost']):
        if M['upper'][j]>0 and total.get(j,F(0))!=F(str(c)):raise ValueError('physical group costs do not partition current objective')
    theta={col:k for k,col in p.a.theta.items()};candidates=set()
    for i,j,value in zip(M['row'],M['col'],M['data']):
        if j in theta and F(str(value))==1 and M['equalities'][i]and F(str(M['rhs'][i]))==0:candidates.add(i)
    actual={i:{}for i in candidates}
    for i,j,value in zip(M['row'],M['col'],M['data']):
        if i in actual:actual[i][j]=actual[i].get(j,F(0))+F(str(value))
    actual={i:{j:c for j,c in row.items()if c}for i,row in actual.items()}
    matched={}
    for k,row in expected.items():
        clean={j:c for j,c in row.items()if c};hits=[i for i,r in actual.items()if r==clean]
        if len(hits)!=1:raise ValueError('unique original stage-cost identity required')
        U=M['upper'][p.a.theta[k]]
        if not isfinite(float(U))or F(str(U))<0:raise ValueError('finite nonnegative stage-cost upper bound required')
        matched[k]=hits[0]
    return costs,matched


def compress(parent,deadline=float('inf')):
    if hasattr(parent,'kappa'):raise ValueError('cost extension cannot be applied twice')
    started=time.perf_counter();costs,identities=_validated_costs(parent)
    q=copy(parent);q.cost_parent=parent;q.parent_columns=len(parent.matrix['cost']);q.physical_costs=costs
    q.stage_cost_identity_rows=identities;q.kappa={};q.costs={}
    q.matrix={key:list(value)if isinstance(value,list)else value for key,value in parent.matrix.items()}
    q.matrix['shape']=list(parent.matrix['shape'])
    for key,row in sorted(costs.items()):
        k,g=key;col=q.matrix['shape'][1];q.kappa[key]=col;q.matrix['shape'][1]+=1
        q.matrix['cost'].append(0);q.matrix['upper'].append(parent.matrix['upper'][parent.a.theta[k]]);q.matrix['integer'].append(0)
        i=q.matrix['shape'][0];q.matrix['shape'][0]+=1;q.matrix['rhs'].append(0);q.matrix['equalities'].append(True)
        for j,c in [(col,F(1))]+[(j,-c)for j,c in row.items()]:
            q.matrix['row'].append(i);q.matrix['col'].append(j);q.matrix['data'].append(c)
        q.costs[key]={col:F(1)}
    q.adapter=numeric_adapter(q.matrix);q.source_vector=lift(q,parent.source_vector)
    from boundary_context_master.matrix import check_vector
    check_vector(q.adapter,q.source_vector)
    if time.perf_counter()>deadline:raise TimeoutError('cost-support extension allowance')
    q.cost_extension_seconds=time.perf_counter()-started
    return q


def lift(q,vector):
    n=len(vector)
    if n==q.base['shape'][1]:
        from position_ownership.model import lift as position_lift
        old=position_lift(q.cost_parent,vector)
    elif n==q.parent_columns:old=list(vector)
    else:raise ValueError('original or parent vector dimensions required')
    out=old+[0]*len(q.kappa)
    for key,col in q.kappa.items():out[col]=sum((c*F(str(old[j]))for j,c in q.physical_costs[key].items()),F(0))
    return out


def project(q,vector):
    if len(vector)!=len(q.matrix['cost']):raise ValueError('cost-extended vector dimensions required')
    return list(vector[:q.parent_columns])


def substitute_cost_row(q,row,key):
    if row.get('kind')!='cost'or key not in q.physical_costs:raise ValueError('original conditional cost row required')
    entries={int(j):F(str(c))for j,c in row['entries'].items()if F(str(c))}
    expected=q.physical_costs[key]
    physical={j:-c for j,c in entries.items()if j<q.physical_columns}
    if physical!=expected:raise ValueError('old row physical cost coefficients differ')
    out=deepcopy(row);out['entries']={j:c for j,c in entries.items()if j>=q.physical_columns}
    col=q.kappa[key]
    if col in out['entries']:raise ValueError('old row unexpectedly contains a cost-total coordinate')
    out['entries'][col]=F(-1)
    return out
