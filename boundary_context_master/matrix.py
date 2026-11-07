"""Sparse boundary-choice formulation and explicit integer-vector decoding."""
from fractions import Fraction as F
from itertools import product
from types import SimpleNamespace
import math
import time
import numpy as np

from fleet_sample_assay.master import Matrix
from .contracts import histogram,order_bits,cost_row,resource_path,dimensions


def put_row(M,key,entries,rhs=0,eq=False):
    M.row(key,float(rhs),eq)
    for col,value in entries.items():
        if value:M.columns[col][key]=float(value)


def build(groups,pads,tasks,limits,batteries,floors,elementary,road,templates,upper,
          *,anchor=0,records=(),matrix=None,theta_objective=True,deadline=math.inf):
    tasks=tuple(tasks);limits=tuple(limits);batteries=tuple(batteries)
    floors=tuple(F(v)for v in floors);elementary=tuple(F(v)for v in elementary)
    templates=tuple(templates);K=len(tasks);expected=dimensions(groups,pads,limits)
    if not K==len(limits)==len(floors)==len(elementary)==len(templates) or len(batteries)!=groups:
        raise ValueError('inconsistent complete master domains')
    if any(type(n)is not int or n<0 for n in tasks)or any(b<0 for b in batteries)or min(floors)<0:
        raise ValueError('invalid task/stock/floor domain')
    if len(road)!=pads or any(len(r)!=pads or min(r)<0 for r in road)or upper<=0:
        raise ValueError('nonnegative complete road matrix and finite upper bound required')
    if not math.isfinite(float(upper)):raise ValueError('finite incumbent upper bound required')
    M=Matrix()if matrix is None else matrix;start_columns=len(M.cost)
    x={};u={};theta={};types=[]
    for k,(lo,hi)in enumerate(limits):
        if time.perf_counter()>=deadline:raise TimeoutError('compact matrix construction deadline')
        ts=tuple(product((0,)if k==0 else range(pads),
                         (0,)if k==K-1 else range(pads),range(lo,hi+1)))
        types.append(ts)
        for g,t in product(range(groups),ts):x[k,g,t]=M.column({},upper=1,integer=True)
        for t,r in product(ts,range(1,groups+1)):u[k,t,r]=M.column({},upper=1,integer=True)
        theta[k]=M.column({},cost=1 if theta_objective else 0,upper=float(upper))
    for k,ts in enumerate(types):
        for g in range(groups):
            put_row(M,('cb','choice',k,g),{x[k,g,t]:1 for t in ts},1,True)
        put_row(M,('cb','task',k),{x[k,g,t]:t[2]for g,t in product(range(groups),ts)},tasks[k],True)
        for t in ts:
            entries={x[k,g,t]:1 for g in range(groups)}
            entries.update({u[k,t,r]:-1 for r in range(1,groups+1)})
            put_row(M,('cb','count',k,t),entries,0,True)
            for r in range(1,groups):
                put_row(M,('cb','order',k,t,r),{u[k,t,r+1]:1,u[k,t,r]:-1})
        put_row(M,('cb','floor',k),{theta[k]:-1},-floors[k])
        entries={x[k,g,t]:F(road[t[0]][t[1]])for g,t in product(range(groups),ts)}
        entries[theta[k]]=-1
        put_row(M,('cb','road_floor',k),entries,-elementary[k])
    for k in range(K-1):
        for g,p in product(range(groups),range(pads)):
            entries={x[k,g,t]:1 for t in types[k]if t[1]==p}
            entries.update({x[k+1,g,t]:-1 for t in types[k+1]if t[0]==p})
            put_row(M,('cb','flow',k,g,p),entries,0,True)
    for g in range(groups):
        put_row(M,('cb','stock',g),{x[k,g,t]:t[2]for k in range(K)for t in types[k]},batteries[g])
    put_row(M,('cb','objective_cut'),{theta[k]:1 for k in range(K)},upper)
    put_row(M,('cb','anchor'),{theta[k]:-1 for k in range(K)},-F(anchor))
    for ri,rec in enumerate(records):
        if 'template'not in rec or 'pattern'not in rec or not rec.get('infeasible',False)and 'bound'not in rec:
            raise ValueError('certified context lower bound or infeasibility required')
        pattern=tuple((tuple(t),int(c))for t,c in rec['pattern']);bits=order_bits(pattern)
        if len(bits)!=groups:raise ValueError('context group count differs')
        for k,key in enumerate(templates):
            if key!=rec['template'] or any(t not in types[k]for t,c in pattern):continue
            if rec.get('infeasible',False):
                put_row(M,('cb','forbid',ri,k),{u[k,t,r]:1 for t,r in bits},groups-1)
            else:
                row=cost_row(pattern,max(floors[k],F(rec['bound'])),floors[k])
                entries={u[k,t,r]:row['coefficient']for t,r in row['support']}
                entries[theta[k]]=-1
                put_row(M,('cb','cost',ri,k),entries,-row['constant'])
    if len(M.cost)-start_columns!=expected['total']:raise ValueError('coordinate inventory differs')
    return SimpleNamespace(M=M,A=M.sparse().tocsr(),x=x,u=u,theta=theta,types=tuple(types),
        groups=groups,pads=pads,tasks=tasks,limits=limits,batteries=batteries,
        templates=templates,floors=floors,elementary=elementary,road=road,upper=upper,
        anchor=F(anchor),start_columns=start_columns)


def encode(m,actions,costs):
    r=resource_path(actions,m.tasks,m.limits,m.batteries,m.pads)
    if len(costs)!=len(m.tasks):raise ValueError('complete stage costs required')
    vector=np.zeros(len(m.M.cost))
    for k,(stage,h,cost)in enumerate(zip(actions,r['contexts'],costs)):
        for g,t in enumerate(stage):vector[m.x[k,g,tuple(t)]]=1
        for t,level in order_bits(h):vector[m.u[k,t,level]]=1
        vector[m.theta[k]]=float(cost)
    return vector


def check_vector(m,vector,*,integer=True,tolerance=1e-6):
    vector=np.asarray(vector,dtype=float)
    if vector.shape!=(len(m.M.cost),)or not np.isfinite(vector).all():
        raise ValueError('invalid finite matrix vector')
    eq=np.asarray(m.M.eq,dtype=bool);rhs=np.asarray(m.M.rhs,dtype=float);a=m.A@vector
    flags=np.asarray(m.M.integer,dtype=bool)
    error=float(max(np.max(abs(a[eq]-rhs[eq]),initial=0),
        np.max(a[~eq]-rhs[~eq],initial=0),np.max(-vector,initial=0),
        np.max(vector-np.asarray(m.M.upper),initial=0),
        np.max(abs(vector[flags]-np.rint(vector[flags])),initial=0)if integer else 0))
    if error>tolerance:raise ValueError(('matrix residual',error))
    return error


def decode(m,vector):
    check_vector(m,vector);actions=[]
    for k,ts in enumerate(m.types):
        stage=[]
        for g in range(m.groups):
            chosen=[t for t in ts if vector[m.x[k,g,t]]>.5]
            if len(chosen)!=1:raise ValueError('unique integer group action required')
            stage.append(chosen[0])
        actions.append(tuple(stage))
    resource_path(actions,m.tasks,m.limits,m.batteries,m.pads)
    return dict(actions=tuple(actions),costs=tuple(float(vector[m.theta[k]])for k in range(len(actions))))
