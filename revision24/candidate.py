"""Exact physical candidate transport; no optimization or carrier repair."""
from copy import deepcopy
from fractions import Fraction
import numpy as np
from review88_physical import physical,milp


def snapshot_prefix(events,seconds):
    copied=deepcopy(events)
    for event in copied:
        event['lower_bound']=0
        if event['released']>seconds:raise ValueError('candidate outside prefix')
    return dict(status='complete',seconds=seconds,primal_events=copied,
        fleets=[event['fleet']for event in copied],future_information=False)


def validate_candidate(inst,fleet):
    check=physical.validate_fleet(inst,fleet['group_plans'])
    if not check['valid']:raise ValueError(('candidate physical',check['errors']))
    saved=fleet.get('check',{})
    for key in ['objective','road_cost','flight_cost','age_cost']:
        if key in saved and Fraction(str(saved[key]))!=check[key]:
            raise ValueError(('candidate cost mismatch',key,saved[key],check[key]))
    return deepcopy(dict(fleet,check=check)),check


def _route_key(route):
    row=route if isinstance(route,dict)else vars(route)
    return tuple(row[k]for k in ('p','launch','q','recovery','job','sample'))


def fleet_vector(m,fleet):
    fleet,check=validate_candidate(m.inst,fleet)
    routes={_route_key(route):i for i,route in enumerate(m.routes)}
    arcs={tuple(arc):i for i,arc in enumerate(m.arcs)}
    if len(routes)!=len(m.routes)or len(arcs)!=len(m.arcs):raise ValueError('ambiguous candidate support')
    x=np.zeros(len(m.M.cost))
    try:
        for g,plan in enumerate(fleet['group_plans']):
            for route in plan['routes']:x[m.rids[g,routes[_route_key(route)]]]=1
            for arc in plan['trail']:x[m.zids[g,arcs[tuple(arc)]]]=1
            for j,s,e,k in plan['assays']:x[m.assignments[g,j,s,e,k]]=1
    except KeyError as exc:raise ValueError(('candidate missing support',exc.args))from exc
    rhs=np.asarray(m.M.rhs);eq=np.asarray(m.M.eq,dtype=bool);activity=m.A@x
    error=max(np.max(abs(activity[eq]-rhs[eq]),initial=0),np.max(activity[~eq]-rhs[~eq],initial=0),
        np.max(-x,initial=0),np.max(x-np.asarray(m.M.upper),initial=0),np.max(abs(x-np.rint(x)),initial=0))
    if not np.all(np.isfinite(x))or error>1e-6:raise ValueError(('candidate matrix',error))
    if abs(float(np.asarray(m.M.cost)@x)-float(check['objective']))>1e-6:raise ValueError('candidate matrix cost')
    restored,recheck=milp.restore(m,x)
    for original,new in zip(fleet['group_plans'],restored['group_plans']):
        if (sorted(_route_key(r)for r in original['routes'])!=sorted(_route_key(r)for r in new['routes'])
            or list(map(tuple,original['trail']))!=list(map(tuple,new['trail']))
            or sorted(map(tuple,original['assays']))!=sorted(map(tuple,new['assays']))):
            raise ValueError('candidate changed in matrix roundtrip')
    if recheck['objective']!=check['objective']:raise ValueError('candidate restored cost')
    return x,restored,recheck


def admit(model,variables,m,fleet):
    raw,restored,check=fleet_vector(m,fleet)
    sol=model.createOrigSol()
    for var,value in zip(variables,raw):model.setSolVal(sol,var,float(value))
    feasible=bool(model.checkSol(sol,printreason=False,completely=True,original=True))
    if not feasible:
        model.freeSol(sol)
        raise ValueError('native rejected candidate feasibility')
    stored=bool(model.addSol(sol,free=True))
    return dict(feasible=feasible,stored=stored,objective=float(check['objective']),
        exact_objective=str(check['objective']),same_physical_schedule=True)
