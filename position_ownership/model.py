"""Saved strong matrix with continuous old auxiliaries and sufficient y/s bits."""
from fractions import Fraction as F
from itertools import product
from types import SimpleNamespace as NS
import time
from service_separator25_joint.cases import split
from service_separator25_joint.solve import signature
from service_separator25_joint.interface import pack_limits
from referee_boundary21_v2.interface import certify
from responsibility_branch.model import check_costs
from vehicle_root_ablation.model import adapter
from .support import recover_indices


def _row(A,matrix,row,expected):
    a,b=A.indptr[row:row+2];actual={int(j):float(v)for j,v in zip(A.indices[a:b],A.data[a:b])if v}
    if actual!={j:float(v)for j,v in expected.items()if v}or not matrix['equalities'][row]or matrix['rhs'][row]!=0:
        raise ValueError(('saved auxiliary mapping differs',row))


def build(inst,cuts,base,physical_shape,original_augmented_rows,source,deadline=float('inf'),checkpoint=None):
    started=time.perf_counter();cuts=tuple(cuts);cert=certify(inst,cuts)
    if (0,)+tuple(c['epoch']for c in cert['cuts']if c['valid'])+(inst.H,)!=cuts:raise ValueError('all boundaries must qualify')
    support=recover_indices(inst,cuts,base,*physical_shape,deadline=deadline,checkpoint=checkpoint)
    A=support.A;G=inst.groups;P=len(inst.pads);stages=tuple((local,ids,signature(local))for local,ids in split(inst,cuts));K=len(stages)
    limits=tuple(pack_limits(local)for local,ids,key in stages);tasks=tuple(len(ids)for local,ids,key in stages)
    x={};u={};theta={};types=[];column=physical_shape[1]
    for k,(lo,hi)in enumerate(limits):
        ts=tuple(product((0,)if k==0 else range(P),(0,)if k==K-1 else range(P),range(lo,hi+1)));types.append(ts)
        for g,t in product(range(G),ts):x[k,g,t]=column;column+=1
        for t,r in product(ts,range(1,G+1)):u[k,t,r]=column;column+=1
        theta[k]=column;column+=1
    if column!=base['shape'][1]:raise ValueError('saved auxiliary column inventory differs')
    if any(base['integer'][j]!=1 for j in (*x.values(),*u.values()))or any(base['integer'][j]for j in theta.values()):
        raise ValueError('unexpected saved auxiliary domains')
    block=G*(2*P+1)+1;mapping_start=original_augmented_rows-K*block
    if mapping_start<physical_shape[0]:raise ValueError('physical/auxiliary row inventory differs')
    ownership={};positions={};costs={};cost_rows={};node_out={};node_in={}
    for i,(p,t,q,v)in enumerate(support.arcs):node_out.setdefault((p,t),[]).append(i);node_in.setdefault((q,v),[]).append(i)
    for k,(local,ids,key)in enumerate(stages):
        a,b=cuts[k:k+2];selected=set(ids)
        for g in range(G):
            for p in range(P):
                for kind,c,side in ((0,a,0),(1,b,1)):
                    indices=node_out.get((p,c),())if c<inst.H else node_in.get((p,c),())
                    e={x[k,g,t]:1 for t in types[k]if t[side]==p};e.update({support.zids[g,i]:-1 for i in indices})
                    _row(A,base,mapping_start+k*block+g*(2*P+1)+2*p+kind,e)
                if k<K-1:positions[k+1,g,p]=tuple(support.zids[g,i]for i in node_out[p,b])
            e={x[k,g,t]:t[2]for t in types[k]}
            for j,original_j in enumerate(ids):
                cols=tuple(col for(gg,i),col in support.rids.items()if gg==g and support.routes[i].job==original_j)
                ownership[k,g,j]=cols;e.update({col:-1 for col in cols})
            _row(A,base,mapping_start+k*block+g*(2*P+1)+2*P,e)
            c={}
            for(gg,i),col in support.zids.items():
                pp,t,q,v=support.arcs[i]
                if gg==g and a<=t<v<=b and base['cost'][col]:c[col]=F(str(base['cost'][col]))
            for(gg,i),col in support.rids.items():
                if gg==g and support.routes[i].job in selected and base['cost'][col]:c[col]=F(str(base['cost'][col]))
            for(gg,j,s,e,t),col in support.assignments.items():
                if gg==g and j in selected and base['cost'][col]:c[col]=F(str(base['cost'][col]))
            costs[k,g]=c
        cost_rows[k]=mapping_start+k*block+G*(2*P+1)
    check_costs(base,costs,theta,cost_rows)
    M=NS(cost=list(base['cost']),upper=list(base['upper']),rhs=list(base['rhs']),eq=list(base['equalities']),integer=list(base['integer']))
    support.M=M;support.A=A.tocsc()
    d=NS(inst=inst,cuts=cuts,stages=stages,templates=tuple(key for local,ids,key in stages))
    old=NS(M=M,A=A,x=x,u=u,theta=theta,types=tuple(types),groups=G,pads=P,tasks=tasks,limits=limits,
           batteries=(inst.batteries,)*G,physical=support)
    matrix={key:list(v)if isinstance(v,list)else v for key,v in base.items()};matrix['shape']=list(base['shape']);matrix['x']=None
    matrix['integer']=[0]*column;scols={};ycols={}
    for mapping,target in ((ownership,scols),(positions,ycols)):
        for key,indices in sorted(mapping.items()):
            col=matrix['shape'][1];target[key]=col;matrix['shape'][1]+=1
            matrix['cost'].append(0);matrix['upper'].append(1);matrix['integer'].append(1)
            row=matrix['shape'][0];matrix['shape'][0]+=1;matrix['rhs'].append(0);matrix['equalities'].append(True)
            for j,c in [(col,1)]+[(j,-1)for j in indices]:matrix['row'].append(row);matrix['col'].append(j);matrix['data'].append(c)
    p=NS(inst=inst,d=d,a=old,m=support,base=base,matrix=matrix,s=scols,y=ycols,ownership=ownership,
         position_support=positions,costs=costs,adapter=adapter(matrix),source=source,physical_columns=physical_shape[1],anchor=F(0))
    from boundary_context_master.common import embed_source
    from boundary_context_master.matrix import check_vector
    vector,_,_=embed_source(d,old,source);p.source_vector=lift(p,vector);check_vector(p.adapter,p.source_vector)
    p.build_seconds=time.perf_counter()-started
    if checkpoint:checkpoint(dict(component='complete_interface',seconds=p.build_seconds,integers=sum(matrix['integer'])))
    if time.perf_counter()>deadline:raise TimeoutError('position interface validation deadline')
    return p


def lift(p,vector):
    if len(vector)!=p.base['shape'][1]:raise ValueError('original vector dimensions')
    out=list(vector)+[0]*(len(p.s)+len(p.y))
    for mapping,columns in ((p.ownership,p.s),(p.position_support,p.y)):
        for key,col in columns.items():out[col]=sum(vector[j]for j in mapping[key])
    return out


def decode(p,vector):
    G=p.inst.groups;K=len(p.d.stages);P=len(p.inst.pads)
    positions={}
    for c in range(1,K):
        for g in range(G):
            values=[vector[p.y[c,g,a]]for a in range(P)]
            if any(abs(v-round(v))>1e-6 for v in values):raise ValueError('fractional boundary position')
            selected=[a for a,v in enumerate(values)if v>.5]
            if len(selected)!=1:raise ValueError('one boundary position required')
            positions[c,g]=selected[0]
    features=[]
    for k,(local,ids,key)in enumerate(p.d.stages):
        for g in range(G):
            values=[vector[p.s[k,g,j]]for j in range(len(ids))]
            if any(abs(v-round(v))>1e-6 for v in values):raise ValueError('fractional task ownership')
            S=tuple(j for j,v in enumerate(values)if v>.5);a=0 if k==0 else positions[k,g];b=0 if k==K-1 else positions[k+1,g]
            lo,hi=p.a.limits[k]
            if not lo<=len(S)<=hi:raise ValueError('ownership count outside original resource bounds')
            endpoints=(()if k==0 else(p.y[k,g,a],))+(()if k==K-1 else(p.y[k+1,g,b],))
            features.append(dict(stage=k,group=g,template=key,subset=S,start=a,end=b,used=len(S),choice=endpoints,
                cost=sum(float(c)*vector[j]for j,c in p.costs[k,g].items())))
    return features
