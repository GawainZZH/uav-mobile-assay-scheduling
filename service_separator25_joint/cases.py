from dataclasses import replace
from review88_physical.types import Instance, Job


def campaign(rounds, *, period=60, assay=6, road=4, drones=2, stations=1,
             width=0, phase_tasks=3, transfer=1, stock_extra=0, seed=0, alternating=False, phase_gap=15,phase_offset=12):
    pads=((0,0),(300,0));waters=((0,100),(300,100))
    jobs=[]
    for r in range(rounds):
        for phase in range(2):
            cohort=phase_tasks[phase] if isinstance(phase_tasks,(tuple,list)) else phase_tasks
            location=1-phase if alternating and r%2 else phase
            for k in range(cohort):
                s=r*period+phase_offset+phase_gap*phase+k
                w=4 if (k+phase+seed)%3==0 else 1
                jobs.append(Job(waters[location],s,s+width,20,
                    r*period+min(50,period-1),w,f'r{r}-p{phase}-{k}'))
    nodes=pads+tuple(j.xy for j in jobs)
    def flight(a,b):
        if a==b:return 0
        # Same engineering reservation matrix as the released G1 family.
        if a<2 and b>=2:return 2 if nodes[b][0]==pads[a][0] else 3
        if b<2 and a>=2:return 2 if nodes[a][0]==pads[b][0] else 3
        return 2 if nodes[a][0]==nodes[b][0] else 3
    ft=tuple(tuple(flight(a,b) for b in range(len(nodes))) for a in range(len(nodes)))
    return Instance(f'S25_r{rounds}_n{len(jobs)}_q{assay}_d{drones}_m{stations}_w{width}_R{road}_s{seed}',
        rounds*period,pads,tuple(jobs),groups=2,drones=drones,
        batteries=(len(jobs)+1)//2+stock_extra,stations=stations,
        assay_duration=assay,ground_duration=4,handoff_delay=transfer,
        endurance=7,flight_fixed=2,air_cost=8,road_cost=1,
        flight_times=ft,road_times=((0,road),(road,0)))


def split(inst,boundaries):
    out=[]
    for a,b in zip(boundaries,boundaries[1:]):
        ids=[j for j,job in enumerate(inst.jobs) if a<=job.lo<b]
        nodes=list(range(len(inst.pads)))+[len(inst.pads)+j for j in ids]
        ft=tuple(tuple(inst.fly(x,y) for y in nodes) for x in nodes)
        jobs=tuple(replace(inst.jobs[j],lo=inst.jobs[j].lo-a,hi=inst.jobs[j].hi-a,
                           deadline=inst.jobs[j].deadline-a) for j in ids)
        out.append((replace(inst,name=inst.name+f'-{a}-{b}',H=b-a,jobs=jobs,
                            batteries=len(ids),flight_times=ft),tuple(ids)))
    return tuple(out)
