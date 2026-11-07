"""Recover original physical indices; never instantiate a physical MILP."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
import time
from scipy.sparse import coo_matrix
from review88_physical.milp import independent_routes


def recover_indices(inst,cuts,matrix,physical_rows,physical_columns,deadline=float('inf'),checkpoint=None):
    started=time.perf_counter();G=inst.groups;P=len(inst.pads);H=inst.H;N=len(inst.jobs)
    if physical_columns%G:raise ValueError('physical group stride differs')
    stride=physical_columns//G;basic=P*(3*H+1)+2*H+1
    A=coo_matrix((matrix['data'],(matrix['row'],matrix['col'])),shape=matrix['shape']).tocsr()
    def clock(label):
        if checkpoint:checkpoint(dict(component=label,seconds=time.perf_counter()-started))
        if time.perf_counter()>deadline:raise TimeoutError('saved support validation deadline: '+label)
    def check(row,expected,rhs=0,eq=False):
        lo,hi=A.indptr[row:row+2];actual={int(j):float(v)for j,v in zip(A.indices[lo:hi],A.data[lo:hi])if v}
        if actual!={j:float(v)for j,v in expected.items()if v}or matrix['rhs'][row]!=rhs or bool(matrix['equalities'][row])!=eq:
            raise ValueError(('saved physical role differs',row))
    # Read the actual saved DAG via its two flow incidences, then verify order.
    road_count=0
    while road_count<stride and not matrix['integer'][road_count]:road_count+=1
    C=A.tocsc();arcs=[]
    for col in range(road_count):
        lo,hi=C.indptr[col:col+2];nodes={}
        for row,value in zip(C.indices[lo:hi],C.data[lo:hi]):
            if row>=basic:continue
            p,t=divmod(int(row),3*H+1)
            if p<P and t<=H:nodes[int(value)]=(p,t)
        if set(nodes)!={-1,1}:raise ValueError('saved road flow incidence missing')
        p,t=nodes[1];q,u=nodes[-1];arcs.append((p,t,q,u))
    expected=tuple((p,t,q,t+d)for p in range(P)for t in range(H)for q in range(P)
        if t+(d:=1 if p==q else inst.road(p,q))<=H)
    if tuple(arcs)!=expected:raise ValueError('saved road DAG ordering/domain differs')
    routes=independent_routes(inst,deadline);variants=sorted({(r.job,r.sample,r.recovery)for r in routes})
    assays=[]
    for j,s,e in variants:
        latest=min(H-inst.assay_duration,s+inst.jobs[j].shelf,inst.jobs[j].deadline-inst.assay_duration)
        assays.extend((j,s,e,k)for k in range(e+inst.handoff_delay,latest+1))
    if road_count+len(routes)+len(assays)!=stride:raise ValueError('saved service column inventory differs')
    zids={(g,i):g*stride+i for g in range(G)for i in range(road_count)}
    rids={(g,i):g*stride+road_count+i for g in range(G)for i in range(len(routes))}
    assignments={(g,*a):g*stride+road_count+len(routes)+i for g in range(G)for i,a in enumerate(assays)}
    arc_index={a:i for i,a in enumerate(arcs)}
    for g in range(G):
        for i,(p,t,q,u)in enumerate(arcs):
            j=zids[g,i];cost=0 if p==q else inst.road_cost*(u-t)
            upper=0 if p!=q and any(t<c<u for c in cuts[1:-1])else 1
            if F(str(matrix['cost'][j]))!=cost or matrix['upper'][j]!=upper or matrix['integer'][j]:raise ValueError('saved road column domain/cost differs')
        for i,r in enumerate(routes):
            j=rids[g,i];cost=inst.flight_fixed+inst.air_cost*(r.recovery-r.launch)
            if F(str(matrix['cost'][j]))!=cost or matrix['upper'][j]!=1 or not matrix['integer'][j]:raise ValueError('saved flight column domain/cost differs')
        for a in assays:
            j,s,e,k=a;col=assignments[g,*a];cost=inst.jobs[j].weight*(k+inst.assay_duration-s)
            if F(str(matrix['cost'][col]))!=cost or matrix['upper'][col]!=1 or not matrix['integer'][col]:raise ValueError('saved assay column domain/cost differs')
    clock('DAG_and_service_indices')
    outgoing={};incoming={}
    for i,(p,t,q,u)in enumerate(arcs):outgoing.setdefault((p,t),[]).append(i);incoming.setdefault((q,u),[]).append(i)
    for g in range(G):
        for p in range(P):
            for t in range(H+1):
                row=g*basic+p*(3*H+1)+t;entries={zids[g,i]:1 for i in outgoing.get((p,t),())}
                for i in incoming.get((p,t),()):entries[zids[g,i]]=entries.get(zids[g,i],0)-1
                check(row,entries,int(p==0 and t==0)-int(p==0 and t==H),True)
        check((g+1)*basic-1,{rids[g,i]:1 for i in range(len(routes))},inst.batteries,N==G*inst.batteries)
    for j in range(N):check(G*basic+j,{rids[g,i]:1 for g in range(G)for i,r in enumerate(routes)if r.job==j},1,True)
    supply_start=G*basic+N;block=len(variants)+len(routes)*inst.ground_duration
    for g in range(G):
        for vi,v in enumerate(variants):
            entries={rids[g,i]:-1 for i,r in enumerate(routes)if(r.job,r.sample,r.recovery)==v}
            entries.update({assignments[g,*a]:1 for a in assays if a[:3]==v})
            check(supply_start+g*block+vi,entries,0,True)
        row=supply_start+g*block+len(variants)
        for i,r in enumerate(routes):
            for t in range(r.recovery,r.recovery+inst.ground_duration):
                check(row,{rids[g,i]:1,zids[g,arc_index[r.q,t,r.q,t+1]]:-1});row+=1
        for t in range(H):
            check(g*basic+P*(3*H+1)+2*t,{rids[g,i]:1 for i,r in enumerate(routes)if r.launch<=t<r.recovery+inst.ground_duration},inst.drones)
            entries={assignments[g,*a]:1 for a in assays if a[3]<=t<a[3]+inst.assay_duration}
            entries.update({zids[g,arc_index[p,t,p,t+1]]:-inst.stations for p in range(P)})
            check(g*basic+P*(3*H+1)+2*t+1,entries)
    clock('flow_cover_supply_ground_occupancy')
    active={}
    for a in assays:
        j,s,e,k=a
        for t in range(k,k+inst.assay_duration):active.setdefault((j,t),[]).append(a)
    event_start=supply_start+G*block;event_stride=2*len(routes)+len(active);event_rows=[]
    for g in range(G):
        for i,r in enumerate(routes):
            for kind,node in enumerate(((r.p,r.launch),(r.q,r.recovery))):
                row=event_start+g*event_stride+2*i+kind;entries={rids[g,i]:1}
                entries.update({zids[g,ai]:-1 for ai in outgoing[node]});check(row,entries);event_rows.append(row)
        for i,((j,t),items)in enumerate(active.items()):
            entries={assignments[g,*a]:1 for a in items}
            entries.update({zids[g,arc_index[p,t,p,t+1]]:-1 for p in range(P)})
            check(event_start+g*event_stride+2*len(routes)+i,entries)
    K=len(cuts)-1
    if event_start+G*event_stride+K*G+K!=physical_rows:raise ValueError('physical row schema differs')
    clock('event_presence_and_sample_parking')
    return NS(inst=inst,arcs=tuple(arcs),routes=routes,zids=zids,rids=rids,assignments=assignments,
              physical_columns=physical_columns,physical_rows=physical_rows,event_rows=event_rows,
              cover_rows=tuple(range(G*basic,G*basic+N)),A=A,basic_rows=basic,
              seconds=time.perf_counter()-started,policy='J')
