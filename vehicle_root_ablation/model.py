"""Saved strong matrix plus append-only rows, with complete source premises."""
from fractions import Fraction as F
from itertools import product
from types import SimpleNamespace as NS
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix


def integer(v):
    q=F(str(v))
    if q.denominator!=1:raise ValueError('integer root matrix required')
    return q.numerator


def append_rows(base,rows):
    out={key:list(value)if isinstance(value,list)else value for key,value in base.items()}
    out['shape']=list(base['shape']);out['x']=None
    for row in rows:
        i=out['shape'][0]
        for j,value in sorted(row['entries'].items()):
            if not 0<=j<out['shape'][1]:raise ValueError('new column forbidden')
            a=integer(value)
            if a:out['row'].append(i);out['col'].append(j);out['data'].append(a)
        out['rhs'].append(integer(row['rhs']));out['equalities'].append(False);out['shape'][0]+=1
    return out


def vehicle_row(theta,boundary,anchor,conditional,groups,name):
    gap=F(conditional)-F(anchor)
    if gap<0:raise ValueError('conditional bound below anchor')
    entries={j:F(-1)for j in theta}
    for j in boundary:entries[j]=entries.get(j,F(0))+gap
    return dict(name=name,entries=entries,rhs=-(F(conditional)-groups*gap))


def point_row(theta,support,groups,bound,floor,name):
    gap=F(bound)-F(floor)
    if gap<0:raise ValueError('context bound below floor')
    entries={j:gap for j in support};entries[theta]=F(-1)
    return dict(name=name,entries=entries,rhs=-(F(bound)-groups*gap))


def exact_vector(matrix,vector,*,integer_flags=False):
    x=[v if isinstance(v,F)else F(str(v))for v in vector]
    if len(x)!=matrix['shape'][1]:raise ValueError('vector dimensions')
    activity=[F(0)]*matrix['shape'][0]
    for i,j,a in zip(matrix['row'],matrix['col'],matrix['data']):
        if x[j]:activity[i]+=integer(a)*x[j]
    for i,(a,b,eq)in enumerate(zip(activity,matrix['rhs'],matrix['equalities'])):
        if (a!=F(str(b))if eq else a>F(str(b))):raise ValueError(('exact root residual',i,str(a-F(str(b)))))
    for j,(v,u,flag)in enumerate(zip(x,matrix['upper'],matrix['integer'])):
        if not 0<=v<=F(str(u))or integer_flags and flag and v.denominator!=1:raise ValueError(('root domain',j))
    return sum((F(str(c))*v for c,v in zip(matrix['cost'],x)),F(0))


def assert_inventory(expected,actual):
    """Exact algebraic equality, independent of COO storage order or splitting."""
    for key in ('shape','rhs','equalities','cost','upper','integer'):
        if expected[key]!=actual[key]:raise ValueError('root inventory differs: '+key)
    def coefficients(m):
        if not len(m['row'])==len(m['col'])==len(m['data']):raise ValueError('malformed COO inventory')
        out={}
        for i,j,a in zip(m['row'],m['col'],m['data']):
            if not (0<=i<m['shape'][0]and 0<=j<m['shape'][1]):raise ValueError('COO coordinate outside matrix')
            out[i,j]=out.get((i,j),0)+integer(a)
        return {key:a for key,a in out.items()if a}
    if coefficients(expected)!=coefficients(actual):raise ValueError('solved matrix coefficients differ')


def adapter(matrix):
    A=coo_matrix((matrix['data'],(matrix['row'],matrix['col'])),shape=matrix['shape']).tocsc()
    return NS(A=A,M=NS(cost=matrix['cost'],upper=matrix['upper'],rhs=matrix['rhs'],
        eq=matrix['equalities'],integer=matrix['integer']))


def build(bundle,method):
    from review88_physical.evidence import instance
    from service_separator25_joint.cases import split
    from service_separator25_joint.solve import signature
    from service_separator25_joint.interface import pack_limits
    from boundary_context_master.contracts import histogram,order_bits
    from group_subset_oracle.trial import checked_price_matrix
    from boundary_context_trial import audit_helper
    if method not in ('BASE','LP','INT'):raise ValueError('fixed method only')
    old=bundle['old_raw'];base=old['common']['augmented_strong_LP']['matrix']
    inst=instance(bundle['global_input']['instance']);cuts=tuple(bundle['global_input']['cuts']);stages=list(split(inst,cuts))
    G=inst.groups;P=len(inst.pads);K=len(stages);column=old['common']['original_shape'][1]
    if (G,P,K)!=(5,3,4):raise ValueError('fixed G5 instance differs')
    xs={};us={};theta={};domains=[];floor=[];audit=audit_helper()
    oldLP=old['common']['augmented_strong_LP']
    certified_anchor=audit.check_price(base,dict(price=oldLP['price'],scale=1,
        exact_dual=oldLP['exact_bound'],bound=oldLP['integer_bound']))
    for k,(local,ids)in enumerate(stages):
        lo,hi=pack_limits(local)
        if (lo,hi)!=(0,4):raise ValueError('type bounds differ')
        ts=tuple(product((0,)if k==0 else range(P),(0,)if k==K-1 else range(P),range(lo,hi+1)))
        domains.append(ts)
        for g,t in product(range(G),ts):xs[k,g,t]=column;column+=1
        for t,r in product(ts,range(1,G+1)):us[k,t,r]=column;column+=1
        theta[k]=column;column+=1
        e=old['envelopes'][signature(local)];floor.append(audit.check_price(e['matrix'],e))
    if column!=base['shape'][1]:raise ValueError('compact coordinate inventory differs')
    origin=old['common']['original_shape'][1];stride=origin//G
    if origin%G:raise ValueError('group stride')
    arcs=tuple((p,t,q,t+d)for p in range(P)for t in range(inst.H)for q in range(P)
        if t+(d:=1 if p==q else inst.road(p,q))<=inst.H)
    pos=[i for i,a in enumerate(arcs)if a[:2]==(0,cuts[1])]
    boundary=tuple(g*stride+i for g in range(G)for i in pos)
    A=coo_matrix((base['data'],(base['row'],base['col'])),shape=base['shape']).tocsr()
    start=old['common']['original_shape'][0]+684;mapping=[]
    for g in range(G):
        row=start+g*(2*P+1)+1;expected={xs[0,g,(0,0,v)]:F(1)for v in range(5)}
        expected.update({g*stride+i:F(-1)for i in pos});a,b=A.indptr[row:row+2]
        actual={int(A.indices[z]):F(str(A.data[z]))for z in range(a,b)}
        if actual!=expected or not base['equalities'][row]or base['rhs'][row]!=0:raise ValueError('physical B mapping differs')
        mapping.append(row)
    cost_rows=[start+k*(G*(2*P+1)+1)+G*(2*P+1)for k in range(K)]
    aggregate=np.asarray(A[cost_rows].sum(axis=0)).ravel()
    expected=-np.asarray(base['cost'],dtype=float)
    for j in theta.values():expected[j]+=1
    active=np.asarray(base['upper'])>0
    if np.max(abs(aggregate[active]-expected[active]),initial=0)>0 or any(not base['equalities'][i]or base['rhs'][i]!=0 for i in cost_rows):
        raise ValueError('stage/global cost identity differs')
    cap=old['common']['original_shape'][0]+682
    a,b=A.indptr[cap:cap+2]
    if {int(A.indices[z]):integer(A.data[z])for z in range(a,b)}!={j:1 for j in theta.values()}or base['rhs'][cap]!=1108:
        raise ValueError('same objective cutoff missing')
    L=F(old['common']['common_integer_anchor'])
    if L!=964:raise ValueError('common anchor differs')
    if certified_anchor!=L:raise ValueError('common anchor price does not certify the stated bound')
    rows=[dict(name='explicit_common_anchor',entries={j:-1 for j in theta.values()},rhs=-L)]
    records=list(old['cuts'])
    single=sum((F(str(bundle[f'price_{1<<j}']['objective']))for j in range(4)),F(0))
    first_key=old['local_events'][8]['kernel'];single_pattern=histogram(zip(*old['local_events'][8]['context']))
    records.append(dict(template=first_key,pattern=single_pattern,bound=single))
    mixed=next(r for r in bundle['group_result']['contexts']if r['source_query']==1)
    if not mixed['optimal']:raise ValueError('new317 context not certified')
    records.append(dict(template=old['local_events'][1]['kernel'],pattern=histogram(zip(*mixed['context'])),bound=F(str(mixed['lower_bound']))))
    for ri,record in enumerate(records):
        pattern=tuple((tuple(t),n)for t,n in record['pattern']);bits=order_bits(pattern)
        for k,(local,ids)in enumerate(stages):
            if signature(local)!=record['template']or any(t not in domains[k]for t,n in pattern):continue
            support=tuple(us[k,t,level]for t,level in bits)
            if record.get('infeasible'):rows.append(dict(name=f'legacy_forbid_{ri}_{k}',entries={j:1 for j in support},rhs=G-1))
            else:rows.append(point_row(theta[k],support,G,max(floor[k],F(str(record['bound']))),floor[k],f'context_{ri}_{k}'))
    Acols=tuple(xs[0,g,(0,0,v)]for g in range(G)for v in (0,1));T=single+sum(floor[1:],F(0));gap=T-L
    entries={j:F(-1)for j in theta.values()};entries.update({j:gap for j in Acols})
    rows.append(dict(name='existing_singleton_OR',entries=entries,rhs=-(T-G*gap)))
    model=bundle['group_model'];alphaR=(52,62,82,59);alphaI=(52,62,96,56);proofs=[]
    from dataclasses import replace
    if signature(replace(stages[0][0],groups=1))!=signature(instance(model['instance'])):
        raise ValueError('priced task data differ from original first block')
    if model['crew_capacity']!=bundle['group_input']['crew_capacity']:raise ValueError('original crew capacity differs')
    if model['cost_scale']!=1:raise ValueError('integer objective grid differs')
    for mask in range(16):
        r=bundle[f'price_{mask}'];S=tuple(j for j in range(4)if mask>>j&1)
        if tuple(r['subset'])!=S or not r.get('optimal'):raise ValueError('subset price not complete')
        if r.get('analytic'):rbound=F(str(r['lower_bound']))
        else:rbound=audit.check_price(checked_price_matrix(model,r,S),r['LP'])
        ibound=F(str(r['lower_bound']))
        if sum(alphaR[j]for j in S)>rbound or sum(alphaI[j]for j in S)>ibound:raise ValueError('coverage cost premise fails')
        proofs.append(dict(mask=mask,rational_lower=rbound,integer_lower=ibound))
    cR=F(sum(alphaR))+sum(floor[1:],F(0));cI=F(sum(alphaI))+sum(floor[1:],F(0))
    baseline_rows=len(rows)
    if method in ('LP','INT'):rows.append(vehicle_row(tuple(theta.values()),boundary,L,cR,G,'new_rational_vehicle_cost'))
    if method=='INT':rows.append(vehicle_row(tuple(theta.values()),boundary,L,cI,G,'new_native_vehicle_cost'))
    out=append_rows(base,rows)
    metadata=dict(method=method,original_shape=base['shape'],shape=out['shape'],baseline_added_rows=baseline_rows,
        appended_rows=rows,theta=theta,boundary_columns=boundary,mapping_rows=mapping,cost_mapping_rows=cost_rows,
        cutoff=1108,anchor=L,remaining_floors=floor[1:],rational_conditional=cR,native_conditional=cI,
        subset_price_checks=proofs,new_variables=0)
    return out,metadata
