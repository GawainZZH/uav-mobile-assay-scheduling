"""Operational gradients generated without any numerical optimizer."""
from math import ceil
from random import Random
import time
from review88_physical.types import Instance,Job,Route
from review88_physical.evidence import plain
from review88_physical.physical import validate_fleet
from referee_boundary21_v2.interface import certify


def configurations():
    base=dict(G=4,P=3,K=4,n=4,B=6,extended=False)
    return[dict(base,name='A')]+[dict(base,**change)for change in(
        dict(name='G3',G=3),dict(name='G6',G=6),dict(name='P2',P=2),dict(name='P4',P=4),
        dict(name='K2',K=2),dict(name='K6',K=6),dict(name='N2',n=2),dict(name='N6',n=6),
        dict(name='B4',B=4),dict(name='M16',extended=True),dict(name='M24',K=6,extended=True))]


def candidate(inst,roles,rounds,period=64):
    began=time.perf_counter();G=inst.groups;loads=[0]*G;regions=[[None]*G for _ in range(rounds)]
    used=[[set(),set()]for _ in range(rounds)];assignment=[None]*len(roles);nodes=0
    def visit(i):
        nonlocal nodes
        nodes+=1
        if nodes>50000:raise ValueError('predeclared constructive assignment cap')
        if i==len(roles):return True
        r,phase,region=roles[i];seen=set()
        for g in sorted(range(G),key=lambda g:(loads[g],g)):
            state=(loads[g],regions[r][g],g in used[r][phase])
            if state in seen:continue
            seen.add(state)
            if loads[g]>=inst.batteries or g in used[r][phase]or regions[r][g]not in(None,region):continue
            old=regions[r][g];regions[r][g]=region;used[r][phase].add(g);loads[g]+=1;assignment[i]=g
            if visit(i+1):return True
            assignment[i]=None;loads[g]-=1;used[r][phase].remove(g);regions[r][g]=old
        return False
    if not visit(0):raise ValueError('fixed operating configuration has no constructive shared candidate')
    plans=[dict(routes=[],assays=[],trail=[])for _ in range(G)];P=len(inst.pads)
    for g,plan in enumerate(plans):
        for r in range(rounds):
            ids=[j for j,(rr,phase,region)in enumerate(roles)if rr==r and assignment[j]==g]
            start=r*period;end=(r+1)*period
            if not ids:
                plan['trail'].extend((0,t,0,t+1)for t in range(start,end));continue
            region=roles[ids[0]][2];pad=0 if region==0 else P-1
            arrive=start+inst.road(0,pad)
            if pad:plan['trail'].append((0,start,pad,arrive))
            local=[]
            for j in ids:
                sample=inst.jobs[j].lo;launch=sample-inst.fly(pad,P+j);recover=sample+inst.sampling_duration+inst.fly(P+j,pad)
                if launch<arrive:raise ValueError('candidate cannot reach its first launch')
                route=Route(pad,launch,pad,recover,j,sample);plan['routes'].append(route);local.append(route)
            available=start
            for route in sorted(local,key=lambda v:(v.recovery,v.job)):
                assay=max(route.recovery+inst.handoff_delay,available);available=assay+inst.assay_duration
                plan['assays'].append((route.job,route.sample,route.recovery,assay))
            leave=max(available,max(v.recovery+inst.ground_duration for v in local))
            plan['trail'].extend((pad,t,pad,t+1)for t in range(arrive,leave))
            back=leave+inst.road(pad,0)
            if back>end:raise ValueError('shared candidate does not finish its nominal service wave')
            if pad:plan['trail'].append((pad,leave,0,back))
            plan['trail'].extend((0,t,0,t+1)for t in range(back,end))
    check=validate_fleet(inst,plans)
    if not check['valid']:raise ValueError(('shared candidate invalid',inst.name,check['errors']))
    return dict(group_plans=plans,check=check),dict(seconds=time.perf_counter()-began,assignment_nodes=nodes,
        native_optimizer_calls=0,optimization_objective_used=False)


def make(config,seed):
    began=time.perf_counter();G,P,K,n,B=(config[k]for k in('G','P','K','n','B'));rng=Random(seed)
    locations={2:(0,3),3:(0,1,3),4:(0,1,2,3)}[P];pads=tuple((300*p,0)for p in locations)
    waters=tuple((300*p,100)for p in range(4));jobs=[];roles=[];period=64;H=period*K
    for r in range(K):
        parity=(r+seed)%2
        for phase in range(2):
            for j in range(n//2):
                region=(j+parity)%2;site=(0 if region==0 else 3)if phase==0 else(1 if region==0 else 2)
                sample=r*period+16+14*phase+j+rng.randint(-2,2)
                deadline=min(H,r*period+(68 if config['extended']else 50))
                jobs.append(Job(waters[site],sample,sample+2,40 if config['extended']else 20,deadline,
                    4 if rng.randrange(3)==0 else 1,f'r{r}-p{phase}-{j}'));roles.append((r,phase,region))
    coordinates=pads+tuple(j.xy for j in jobs)
    ft=tuple(tuple(0 if a==b else 1+ceil((abs(x[0]-y[0])+abs(x[1]-y[1]))/300)
        for b,y in enumerate(coordinates))for a,x in enumerate(coordinates))
    roads=tuple(tuple(4*abs(a-b)for b in locations)for a in locations)
    inst=Instance(f'RSC_{config["name"]}_s{seed}',H,pads,tuple(jobs),groups=G,drones=2,batteries=B,
        stations=1,assay_duration=6,ground_duration=4,handoff_delay=1,endurance=7,
        flight_fixed=2,air_cost=8,road_cost=1,flight_times=ft,road_times=roads)
    nominal=tuple(range(0,H+1,period));proposed=certify(inst,nominal);cuts=nominal
    while True:
        proof=certify(inst,cuts)
        if proof['valid']:break
        valid={c['epoch']for c in proof['cuts']if c['valid']};new=(0,)+tuple(c for c in cuts[1:-1]if c in valid)+(H,)
        if new==cuts:raise ValueError('boundary normalization did not make progress')
        cuts=new
    witness,construction=candidate(inst,roles,K)
    sizes=[sum(a<=j.lo<b for j in jobs)for a,b in zip(cuts,cuts[1:])]
    metadata=dict(configuration=config['name'],seed=seed,parameters=config,proposed_cuts=len(nominal)-2,
        qualified_cuts=len(cuts)-2,largest_block=max(sizes),blocks=len(sizes),block_sizes=sizes,
        total_tasks=len(jobs),total_stock=G*B,road_diameter=proposed['maximum_shortest_road'],
        qualification=proposed,final_qualification=proof,candidate_construction=construction,
        generator_seconds=time.perf_counter()-began)
    return dict(raw=plain(dict(instance=inst,cuts=cuts)),candidate=plain(witness),metadata=plain(metadata))
