"""Project public physical incumbents; never import source bounds or trees.

Physical projection does not by itself validate extra branch/shared rows. A
fragment must pass its current kernel's complete matrix before it is admitted.
Group permutation is for the homogeneous Instance domain, not heterogeneous
equipment, fixed named responsibilities, or differing verification calendars.
"""
from dataclasses import replace
from fractions import Fraction
import hashlib
import json

from review88_physical.types import Route
from review88_physical.physical import validate_fleet
from review88_physical.evidence import plain
from service_separator25_joint.cases import split


def _signature(inst):
    data=plain(inst);data.pop('name')
    for job in data['jobs']:job.pop('name')
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()


def project_fleet(inst,cuts,fleet):
    """Return complete joint-stage primal fragments and exact additive costs."""
    cuts=tuple(cuts)
    if len(cuts)<2 or any(type(c) is not int for c in cuts) or cuts[0]!=0 or cuts[-1]!=inst.H or any(
            a>=b for a,b in zip(cuts,cuts[1:])):
        raise ValueError('ordered integer cuts must span the entire horizon')
    plans=fleet['group_plans'];check=validate_fleet(inst,plans)
    if not check['valid']:raise ValueError(('full primal is not physically valid',check['errors']))
    for field in ('objective','road_cost','flight_cost','age_cost'):
        saved=fleet.get('check',{}).get(field)
        if saved is not None and Fraction(str(saved))!=check[field]:
            raise ValueError('saved primal cost disagrees with physical replay')
    routes=[];positions=[]
    for plan in plans:
        rr=tuple(Route(**r) if isinstance(r,dict) else r for r in plan['routes'])
        nodes={}
        for p,t,q,u in plan['trail']:nodes[t]=p;nodes[u]=q
        if any(c not in nodes for c in cuts):raise ValueError('primal must be parked at each cut')
        for c in cuts[1:-1]:
            if any(r.launch<c<r.recovery+inst.ground_duration for r in rr):
                raise ValueError('flight or ground obligation crosses a cut')
            if any(e<=c<k+inst.assay_duration for j,s,e,k in plan['assays']):
                raise ValueError('pending sample or assay crosses a cut')
        routes.append(rr);positions.append(nodes)
    if any(not any(a<=job.lo<b for job in inst.jobs) for a,b in zip(cuts,cuts[1:])):
        raise ValueError('empty task blocks are not kernel templates')

    # This independent boundary replayer adds/subtracts real external depot access.
    from service_separator25_joint.kernel import validate_boundary
    fragments=[]
    for stage,((a,b),(local,ids)) in enumerate(zip(zip(cuts,cuts[1:]),split(inst,cuts))):
        remap={j:i for i,j in enumerate(ids)};local_plans=[];triples=[]
        for g,plan in enumerate(plans):
            selected=[r for r in routes[g] if r.job in remap]
            if any(not a<=r.launch<r.recovery+inst.ground_duration<=b for r in selected):
                raise ValueError('selected task flight/ground is outside its block')
            assays=[(remap[j],s-a,e-a,k-a) for j,s,e,k in plan['assays'] if j in remap]
            if any(not 0<=k<k+inst.assay_duration<=b-a for j,s,e,k in assays):
                raise ValueError('selected assay is outside its block')
            local_plans.append(dict(
                routes=tuple(replace(r,job=remap[r.job],launch=r.launch-a,
                                     sample=r.sample-a,recovery=r.recovery-a) for r in selected),
                trail=tuple((p,t-a,q,u-a) for p,t,q,u in plan['trail'] if a<=t<u<=b),
                assays=tuple(assays)))
            triples.append((positions[g][a],positions[g][b],len(selected)))
        order=tuple(sorted(range(inst.groups),key=lambda g:triples[g]))
        context=tuple(tuple(triples[g][i] for g in order) for i in range(3))
        normalized=tuple(local_plans[g] for g in order)
        local_check=validate_boundary(local,normalized,context[0],context[1])
        if not local_check['valid']:raise ValueError(('fragment physical replay',local_check['errors']))
        fragments.append(dict(stage=stage,kernel=_signature(local),instance=local,context=context,
            objective=local_check['objective'],fleet=dict(group_plans=normalized,check=local_check),
            native_calls=0,source_lower_bound=None))
    if sum(p['objective'] for p in fragments)!=check['objective']:
        raise ValueError('primal fragment costs do not sum to the full physical cost')
    return tuple(fragments)


def _route_key(route):
    r=Route(**route) if isinstance(route,dict) else route
    return r.p,r.launch,r.q,r.recovery,r.job,r.sample


def embed_fragment(kernel,fragment):
    """Check a physical fragment against every row of this fixed kernel."""
    import numpy as np
    from service_separator25_joint.kernel import validate_boundary
    m=kernel.m;ctx=fragment['context']
    if fragment['kernel']!=_signature(kernel.inst):
        raise ValueError('fragment and kernel templates differ')
    if len(ctx)!=3 or any(len(v)!=m.inst.groups for v in ctx) or any(
            type(v) is not int for values in ctx for v in values):
        raise ValueError('invalid integer context')
    start,end,consumed=ctx
    if any(not 0<=p<len(m.inst.pads) for p in start+end) or any(c<0 for c in consumed):
        raise ValueError('context outside pad or pack domain')
    plans=fragment['fleet']['group_plans']
    checked=validate_boundary(kernel.inst,plans,start,end)
    if not checked['valid']:raise ValueError(('invalid fragment boundary replay',checked['errors']))
    if checked['objective']!=fragment['objective']:
        raise ValueError('fragment cost disagrees with boundary replay')
    routes={_route_key(r):i for i,r in enumerate(m.routes)}
    arcs={tuple(arc):i for i,arc in enumerate(m.arcs)}
    if len(routes)!=len(m.routes) or len(arcs)!=len(m.arcs):
        raise ValueError('ambiguous kernel support')
    x=np.zeros(len(kernel.cost),dtype=float)
    try:
        for g,plan in enumerate(plans):
            for r in plan['routes']:x[m.rids[g,routes[_route_key(r)]]]=1
            for arc in plan['trail']:x[m.zids[g,arcs[tuple(arc)]]]=1
            for j,s,e,k in plan['assays']:x[m.assignments[g,j,s,e,k]]=1
    except KeyError as exc:
        raise ValueError(('fragment has no current matrix column',exc.args)) from exc
    rhs=np.asarray(kernel.rhs(*ctx),dtype=float);eq=kernel.eq;activity=m.A@x
    error=max(np.max(abs(activity[eq]-rhs[eq]),initial=0),
              np.max(activity[~eq]-rhs[~eq],initial=0),
              np.max(-x,initial=0),np.max(x-kernel.upper,initial=0),
              np.max(abs(x-np.rint(x)),initial=0),
              abs(float(kernel.cost@x)-float(checked['objective'])))
    if not np.isfinite(x).all() or not np.isfinite(rhs).all() or error>1e-6:
        raise ValueError(('fragment matrix residual',error))
    return dict(vector=x,rhs=rhs,row_lower=np.where(eq,rhs,-np.inf),
                objective=checked['objective'],context=ctx,
                fleet=dict(group_plans=plans,check=checked),native_output=False)


def seed_kernel(kernel,fragment):
    """Retain a checked upper bound; source information never raises a lower bound."""
    import math
    embedded=embed_fragment(kernel,fragment);ctx=embedded['context']
    upper=float(embedded['objective']);lower=float(kernel.lower(*ctx));old=kernel.solved.get(ctx)
    if not math.isfinite(lower) or lower>upper+1e-6 or old and old['status']=='infeasible':
        raise ValueError('valid fragment contradicts the existing context bound')
    current=dict(status='unresolved',optimal=False,lower_bound=lower,objective=upper,
        context=ctx,fleet=embedded['fleet'],vector=embedded['vector'].tolist(),native_vector=None,
        seconds=0,nodes=0,incumbent_origin='shared_primal_fragment')
    if old and old.get('objective') is not None and old['objective']<=upper:
        current=old.copy();current['lower_bound']=max(current.get('lower_bound',0),lower)
    if current['lower_bound']>current['objective']+1e-6:
        raise ValueError('context bound exceeds the best retained primal')
    if current['lower_bound']>=current['objective']:
        current.update(status='optimal',optimal=True)
    kernel.mip.register_seed(embedded['row_lower'],embedded['rhs'],embedded['vector'],upper)
    kernel.solved[ctx]=current;kernel.attempts.setdefault(ctx,0)
    return dict(context=ctx,source_upper=upper,retained_upper=current['objective'],
                retained_lower=current['lower_bound'],closed=current['optimal'],native_calls=0)
