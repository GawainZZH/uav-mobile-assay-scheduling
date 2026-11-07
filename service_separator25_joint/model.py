"""Carrier-event activation and constructive elimination of vehicle binaries."""
import math
import numpy as np
from review88_physical.milp import build as original_build


def build(inst,policy='J',deadline=None):
    m=original_build(inst,policy,deadline);M=m.M
    outgoing={};waiting={}
    for i,(p,t,b,u) in enumerate(m.arcs):
        outgoing.setdefault((p,t),[]).append(i)
        if p==b:waiting.setdefault(t,[]).append(i)
    for g in range(inst.groups):
        for i,r in enumerate(m.routes):
            for kind,node in [('launch',(r.p,r.launch)),('return',(r.q,r.recovery))]:
                key=('event_presence',kind,g,i);M.row(key)
                M.columns[m.rids[g,i]][key]=1
                for ai in outgoing[node]:M.columns[m.zids[g,ai]][key]=-1
        active={}
        for (gg,j,s,e,k),col in m.assignments.items():
            if gg!=g:continue
            for t in range(k,k+inst.assay_duration):active.setdefault((j,t),[]).append(col)
        for (j,t),cols in active.items():
            key=('sample_parked',g,j,t);M.row(key)
            for col in cols:M.columns[col][key]=1
            for ai in waiting[t]:M.columns[m.zids[g,ai]][key]=-1
    for col in m.zids.values():M.integer[col]=0
    m.A=M.sparse().tocsc();return m


def complete(m,vector,start,end):
    """Fix service events; min-cost unit carrier paths preserve every event."""
    x=np.asarray(vector,dtype=float).copy();old=np.array(m.M.cost)@x
    services=[*m.rids.values(),*m.assignments.values()]
    if np.max(abs(x[services]-np.rint(x[services])),initial=0)>1e-6:raise ValueError('fractional service event')
    for g in range(m.inst.groups):
        required={0:start[g],m.inst.H:end[g]};active=np.zeros(m.inst.H,dtype=int)
        pending=np.zeros(m.inst.H,dtype=int)
        def require(t,p):
            if t in required and required[t]!=p:raise ValueError('incompatible required carrier nodes')
            required[t]=p
        for i,r in enumerate(m.routes):
            if x[m.rids[g,i]]>.5:
                require(r.launch,r.p);require(r.recovery,r.q)
                for t in range(r.recovery,r.recovery+m.inst.ground_duration+1):require(t,r.q)
        for (gg,j,s,e,k),col in m.assignments.items():
            if gg==g and x[col]>.5:
                active[k:k+m.inst.assay_duration]=1
                if m.policy=='QC':pending[e:k+m.inst.assay_duration]=1
        required_prefix=np.zeros(m.inst.H+2,dtype=int)
        for t in required:required_prefix[t+1]=1
        required_prefix=np.cumsum(required_prefix);active_prefix=np.r_[0,np.cumsum(active)]
        outgoing={}
        for i,(p,t,b,u) in enumerate(m.arcs):
            if m.M.upper[m.zids[g,i]]==0:continue
            if t in required and required[t]!=p:continue
            if u in required and required[u]!=b:continue
            if required_prefix[u]-required_prefix[t+1]>0:continue
            if p!=b and active_prefix[u]-active_prefix[t]>0:continue
            if p!=b and pending[t]:continue
            outgoing.setdefault((p,t),[]).append(i)
        distance={(start[g],0):0};previous={}
        for t in range(m.inst.H):
            for p in range(len(m.inst.pads)):
                node=(p,t)
                if node not in distance:continue
                for i in outgoing.get(node,()):
                    a=m.arcs[i];nxt=(a[2],a[3]);cost=distance[node]+m.M.cost[m.zids[g,i]]
                    if cost<distance.get(nxt,math.inf):distance[nxt]=cost;previous[nxt]=(node,i)
        node=(end[g],m.inst.H)
        if node not in distance:raise ValueError('no integral carrier completion')
        for (gg,i),col in m.zids.items():
            if gg==g:x[col]=0
        while node!=(start[g],0):
            node,i=previous[node];x[m.zids[g,i]]=1
    if np.array(m.M.cost)@x>old+1e-6:raise ValueError('carrier completion increased cost')
    return x
