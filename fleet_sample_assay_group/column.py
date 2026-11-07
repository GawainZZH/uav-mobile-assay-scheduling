"""Decode anonymous pending-sample actions and independently replay a group."""
from dataclasses import replace
from fractions import Fraction as F
from fleet_sample_assay.physical import vehicle_arcs
from fleet_sample_assay.check import integer_schedule
from .exact import rational


def decode(inst,routes,actions,required,forbidden):
    selected=[routes[i] for a in actions for i in (a[3] if a[0]=='wait' else a[5])]
    covered=tuple(sorted(j for r in selected for j,s in r.visits))
    if len(covered)!=len(set(covered)) or not set(required)<=set(covered) or set(forbidden)&set(covered):
        raise ValueError('decoded coverage/branch requirements')
    recover={};pending=[];assign=[];trail=[];node=(0,0)
    for r in selected:recover.setdefault(r.recovery,[]).append(r)
    for a in actions:
        kind,t,p=a[:3]
        if (p,t)!=node:raise ValueError('decoded vehicle path discontinuity')
        for r in recover.get(t,[]):
            if r.q!=p:raise ValueError('decoded recovery location')
            for j,s in r.visits:
                e,d=r.window(inst,j,s);pending.append((inst.jobs[j].weight,d,j,s,e))
        if kind=='move':
            q,u=a[3:5];arc=(0,p,t,q,u,(u-t)*inst.road_cost)
        else:
            for weight,expiry in a[4]:
                candidates=[x for x in pending if x[:2]==(weight,expiry)]
                if not candidates:raise ValueError('anonymous assay cannot be decoded')
                chosen=min(candidates);pending.remove(chosen);w,d,j,s,e=chosen
                assign.append(((0,j,s,e,t),1))
            arc=(0,p,t,p,t+1,0)
        trail.append(arc);node=arc[3],arc[4]
    if node!=(0,inst.H) or pending:raise ValueError('decoded unfinished group schedule')
    nodes={(0,0)}|{(a[3],a[4]) for a in trail}
    if any((r.p,r.launch) not in nodes or (r.q,r.recovery) not in nodes for r in selected):
        raise ValueError('decoded endpoint must be a vehicle path node')
    road=sum((rational(a[5]) for a in trail),F(0));flight=sum((rational(r.cost(inst)) for r in selected),F(0))
    assay=sum((rational(inst.jobs[j].weight*(k+1)) for (g,j,s,e,k),v in assign),F(0))
    cost=road+flight+assay;all_arcs=vehicle_arcs(inst);chosen=set(trail)
    op=dict(status=0,objective=float(cost),z_values=[int(a in chosen) for a in all_arcs],arcs=all_arcs,route_values=[1]*len(selected))
    net=dict(status=0,cost=float(assay),assignments=assign)
    if covered:
        order={j:i for i,j in enumerate(covered)};sub=replace(inst,jobs=tuple(inst.jobs[j] for j in covered))
        remapped=[replace(r,visits=tuple((order[j],s) for j,s in r.visits)) for r in selected]
        remap_net={**net,'assignments':[((g,order[j],s,e,k),v) for (g,j,s,e,k),v in assign]}
        check=integer_schedule(sub,remapped,op,remap_net,tol=1e-6)
    else:
        if selected or assign:raise ValueError('empty column has service')
        if any(a not in all_arcs for a in trail):raise ValueError('empty path has nonphysical arc')
        check=dict(valid=True,cost=float(cost),flights=0,vehicle_paths=1)
    return dict(covered=covered,cost=float(cost),exact_cost=(str(cost.numerator),str(cost.denominator)),
        routes=tuple(selected),trail=tuple(trail),operation=op,network=net,check=check)
