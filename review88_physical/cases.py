"""Predeclared release grid, eight assignment bridges and four correctness cases."""
from dataclasses import replace
from math import ceil,hypot
from .types import Instance,Job,Route


def _instance(name,gap,offset,road,sites,epochs,q,weight,groups=1,drones=1,batteries=6):
    pads=((0,0),(gap,0));waters=((0,offset),(gap,offset))
    jobs=tuple(Job(waters[p],s,s,20,90,weight if j==0 else 1,str(j))
               for j,(p,s) in enumerate(zip(sites,epochs)))
    nodes=pads+tuple(j.xy for j in jobs)
    ft=tuple(tuple(0 if a==b else 1+ceil(hypot(x[0]-y[0],x[1]-y[1])/300)
                   for b,y in enumerate(nodes)) for a,x in enumerate(nodes))
    return Instance(name,90,pads,jobs,groups=groups,drones=drones,batteries=batteries,
        assay_duration=q,ground_duration=4,handoff_delay=1,endurance=7,flight_fixed=2,
        air_cost=8,road_cost=1,flight_times=ft,road_times=((0,road),(road,0)))


def operations_cases():
    return tuple(_instance(f'A_G{g}_q{q}_w{w}_ell1',gap,offset,road,
        (0,1,1,0,0,1),(2,16,27,44,55,68),q,w)
        for g,(gap,offset,road) in enumerate(((300,100,4),(500,150,5)),1)
        for q in (6,8,10) for w in (1,4))


def bridge_cases():
    out=[]
    for n in (6,8):
        sites=(0,1,1,0,0,1) if n==6 else (0,1,1,1,0,0,1,1)
        epochs=(2,16,27,44,55,68) if n==6 else (2,16,27,37,44,55,61,68)
        for q in (6,10):
            for d in (1,2):out.append(_instance(f'B_n{n}_q{q}_D{d}',500,150,5,
                sites,epochs,q,4,groups=2,drones=d,batteries=n//2))
    return tuple(out)


def cases():return operations_cases()+bridge_cases()


def witness(inst):
    groups=[dict(routes=[],trail=[(0,t,0,t+1) for t in range(inst.H)],assays=[])
            for _ in range(inst.groups)]
    if inst.name.startswith('T'):
        for j,job in enumerate(inst.jobs):
            g=j%inst.groups;s=job.lo;l=s-inst.fly(0,len(inst.pads)+j)
            e=s+inst.sampling_duration+inst.fly(len(inst.pads)+j,0)
            groups[g]['routes'].append(Route(0,l,0,e,j,s))
            groups[g]['assays'].append((j,s,e,e+inst.handoff_delay))
    else:
        for j,job in enumerate(inst.jobs):
            g=j%inst.groups;s=job.lo;l=s-inst.fly(0,len(inst.pads)+j)
            e=s+inst.sampling_duration+inst.fly(len(inst.pads)+j,0)
            groups[g]['routes'].append(Route(0,l,0,e,j,s))
            groups[g]['assays'].append((j,s,e,e+inst.handoff_delay))
    return dict(group_plans=groups)


def tiny_cases():
    a=Instance('T1_G1_ell0',12,((0,0),),(Job((1,0),1,1,8,12,2),),
        batteries=1,assay_duration=2,ground_duration=4,endurance=3)
    b=replace(a,name='T2_G1_ell1',handoff_delay=1)
    c=Instance('T3_G2_ell1',15,((0,0),),(Job((1,0),1,1,9,15,2),Job((1,0),6,6,8,15)),
        groups=2,batteries=1,assay_duration=3,ground_duration=4,endurance=3,handoff_delay=1)
    return a,b,c,replace(c,name='T4_G2_D2_ell1',drones=2)
