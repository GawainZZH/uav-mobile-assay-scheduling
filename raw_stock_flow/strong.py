"""Fresh complete strengthened physical MILP; no saved numerical models."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
from itertools import product
import numpy as np
from service_separator25_joint.solve import signature
from boundary_context_master.common import describe,build_base,augment,project_actions
from boundary_context_master.matrix import encode,check_vector
from review88_physical.milp import dump_matrix,restore
from review88_physical.types import Route
from block_price_propagation.full_sync import full_rows
from block_price_tree.functional import affine_rows
from cost_support.model import numeric_adapter


def road_lower(inst):
    from singleton_context_oracle import road_data
    return tuple(tuple(inst.road_cost*v for v in row)for row in road_data(inst)['distance'])


def append(matrix,rows):
    out={k:list(v)if isinstance(v,list)else v for k,v in matrix.items()};out['shape']=list(matrix['shape']);out['x']=None
    for r in rows:
        i=out['shape'][0];out['shape'][0]+=1;out['rhs'].append(F(str(r['rhs'])));out['equalities'].append(bool(r.get('eq',False)))
        for j,c in r['entries'].items():
            j=int(j);c=F(str(c))
            if not 0<=j<out['shape'][1]:raise ValueError('shared physical coordinate outside model')
            if c:out['row'].append(i);out['col'].append(j);out['data'].append(c)
    return out


def build(initial,deadline):
    inst=initial.domain.inst;K=len(initial.domain.tasks);G=inst.groups;P=len(inst.pads)
    objects={signature(local):NS(envelope={'bound':floor},elementary=floor)
        for(local,ids),floor in zip(initial.domain.stages,initial.floors)}
    d=describe(inst,initial.domain.cuts,objects=objects,deadline=deadline)
    # A physical reservation table need not be metric; only the lower-bound
    # table uses shortest paths. Actual physical arcs retain their input costs.
    d.road=road_lower(inst)
    cap=G*inst.H*inst.road_cost+len(inst.jobs)*(inst.flight_fixed+inst.air_cost*inst.endurance)+inst.H*sum((j.weight for j in inst.jobs),F(0))
    upper=initial.source['check']['objective']if initial.source else cap
    m=build_base(d,deadline);physical_columns=len(m.M.cost);a=augment(d,m,max(F(1),upper),anchor=initial.floor,deadline=deadline)
    matrix=dump_matrix(m,np.asarray(m.M.cost,dtype=float),None);old_columns=len(matrix['cost']);matrix['integer'][physical_columns:]=[0]*(old_columns-physical_columns)
    s={};y={};costs={};physical_costs={};ownership={};positions={};rows=[]
    def variable(bound,integer=False):
        j=matrix['shape'][1];matrix['shape'][1]+=1;matrix['cost'].append(0);matrix['upper'].append(F(bound));matrix['integer'].append(int(integer));return j
    for k,(local,ids,key)in enumerate(d.stages):
        for g,j in product(range(G),range(len(ids))):
            cols=tuple(col for(gg,i),col in m.rids.items()if gg==g and m.routes[i].job==ids[j])
            col=variable(1);s[k,g,j]=col;ownership[k,g,j]=cols
            rows.append(dict(entries={col:1,**{c:-1 for c in cols}},rhs=0,eq=True,kind='raw_ownership_definition'))
    for c,g,p in product(range(1,K),range(G),range(P)):
        cols=tuple(col for(gg,i),col in m.zids.items()if gg==g and m.arcs[i][0]==p and m.arcs[i][1]==d.cuts[c])
        col=variable(1,True);y[c,g,p]=col;positions[c,g,p]=cols
        rows.append(dict(entries={col:1,**{c:-1 for c in cols}},rhs=0,eq=True,kind='raw_position_definition'))
    for k,(local,ids,key)in enumerate(d.stages):
        lo,hi=d.cuts[k:k+2];selected=set(ids)
        for g in range(G):
            c={}
            for(gg,i),col in m.zids.items():
                p,t,b,u=m.arcs[i]
                if gg==g and lo<=t<u<=hi and m.M.cost[col]:c[col]=F(m.M.cost[col])
            for(gg,i),col in m.rids.items():
                if gg==g and m.routes[i].job in selected and m.M.cost[col]:c[col]=F(m.M.cost[col])
            for(gg,j,sample,e,start),col in m.assignments.items():
                if gg==g and j in selected and m.M.cost[col]:c[col]=F(m.M.cost[col])
            C=local.H*local.road_cost+len(ids)*(local.flight_fixed+local.air_cost*local.endurance)+local.H*sum((j.weight for j in local.jobs),F(0))
            col=variable(min(C,upper));costs[k,g]={col:F(1)};physical_costs[k,g]=c
            rows.append(dict(entries={col:1,**{j:-v for j,v in c.items()}},rhs=0,eq=True,kind='raw_cost_definition'))
    matrix=append(matrix,rows)
    q=NS(inst=inst,d=d,m=m,a=a,matrix=matrix,s=s,y=y,costs=costs,physical_costs=physical_costs,
        ownership=ownership,positions=positions,old_columns=old_columns,physical_columns=physical_columns,initial=initial)
    learned=full_rows(q,initial.domain.families);q.matrix=append(q.matrix,learned)
    if initial.source and upper==0:q.matrix=append(q.matrix,[dict(entries={j:c for j,c in enumerate(q.matrix['cost'])if c},rhs=0,kind='zero_cost_source_cutoff')])
    q.adapter=numeric_adapter(q.matrix);q.analytic_rows=learned;q.source_vector=embed(q,initial.source)if initial.source else None
    return q


def embed(q,fleet):
    actions,fragments=project_actions(q.d,fleet);old=encode(q.a,actions,[f['objective']for f in fragments]);m=q.m
    def routekey(r):
        r=Route(**r)if isinstance(r,dict)else r
        return(r.p,r.launch,r.q,r.recovery,r.job,r.sample)
    routes={routekey(r):i for i,r in enumerate(m.routes)};arcs={tuple(a):i for i,a in enumerate(m.arcs)}
    for g,plan in enumerate(fleet['group_plans']):
        for r in plan['routes']:old[m.rids[g,routes[routekey(r)]]]=1
        for a in plan['trail']:old[m.zids[g,arcs[tuple(a)]]]=1
        for j,s,e,k in plan['assays']:old[m.assignments[g,j,s,e,k]]=1
    x=list(old)+[F(0)]*(q.matrix['shape'][1]-len(old))
    for key,col in q.s.items():x[col]=sum(F(str(old[j]))for j in q.ownership[key])
    for key,col in q.y.items():x[col]=sum(F(str(old[j]))for j in q.positions[key])
    for key,mapping in q.costs.items():
        col=next(iter(mapping));x[col]=sum((v*F(str(old[j]))for j,v in q.physical_costs[key].items()),F(0))
    check_vector(q.adapter,x);return x


def complete(q,vector):
    from service_separator25_joint.model import complete as carrier
    check_vector(q.adapter,vector);m=q.m;saved=list(m.M.upper)
    try:
        for c,g,p in q.y:
            if vector[q.y[c,g,p]]>.5:
                for(gg,i),col in m.zids.items():
                    a,t,b,u=m.arcs[i]
                    if gg==g and t==q.d.cuts[c]and a!=p:m.M.upper[col]=0
        raw=np.asarray(vector[:q.old_columns],dtype=float).copy();services=list(m.rids.values())+list(m.assignments.values())
        raw[services]=np.rint(raw[services]);done=carrier(m,raw,(0,)*q.inst.groups,(0,)*q.inst.groups)
        fleet,physical=restore(m,done)
    finally:m.M.upper[:]=saved
    actual=embed(q,fleet)
    if float(physical['objective'])>sum(float(c)*float(v)for c,v in zip(q.matrix['cost'],vector))+1e-5:
        raise ValueError('physical completion increased objective')
    return fleet,actual


def synchronize(q,bands,prices):
    rows=full_rows(q,bands)+affine_rows(q,prices)
    return dict(matrix=append(q.matrix,rows),new_rows=rows)


def coordinates(q):return dict(s=list(q.s.items()),y=list(q.y.items()),costs=list(q.costs.items()))


def transfer(saved,cuts,learned):
    from review88_physical.evidence import instance
    from service_separator25_joint.cases import split
    from block_price_propagation.full_sync import export
    inst=instance(saved['working_instance']);features=saved['feature_coordinates']
    target=NS(inst=inst,d=NS(stages=split(inst,cuts)),
        s={tuple(k):j for k,j in features['s']},y={tuple(k):j for k,j in features['y']},
        costs={tuple(k):{int(j):F(str(v))for j,v in e.items()}for k,e in features['costs']})
    rows=full_rows(target,learned['matrix']['families'])+affine_rows(target,learned['shared_prices'])
    return export(saved['matrix'],rows,saved['source_vector'])
