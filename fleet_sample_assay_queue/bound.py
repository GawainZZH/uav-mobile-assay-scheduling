"""Integer-valid batch cost row and exact finite-box certificates."""
import copy,math,time
from dataclasses import replace
from fractions import Fraction as F
import numpy as np
from scipy.sparse import csr_matrix,vstack

def queue_penalty(weights,stations):
    if type(stations) is not int or stations<1:raise ValueError('positive integer stations')
    ww=[F(w) for w in weights]
    if any(w<0 for w in ww):raise ValueError('negative weight')
    return sum((w*(i//stations) for i,w in enumerate(sorted(ww,reverse=True))),F(0))

def route_bound(inst,r):
    weights=[inst.jobs[j].weight for j,s in r.visits]
    return (r.recovery+1)*sum(weights)+queue_penalty(weights,inst.stations)

def strengthen_pattern(inst,model):
    """Copy metadata, retain all old coordinates/rows and the cost equality."""
    M=copy.deepcopy(model.M);M.row(('batch_queue',))
    ids=list(model.rids)+[model.theta]
    values=[float(route_bound(inst,r)) for r in model.routes]+[-1.]
    row=csr_matrix((values,([0]*len(ids),ids)),shape=(1,model.A.shape[1]))
    return replace(model,M=M,A=vstack((model.A,row),format='csc'))

def strengthen_compact(inst,B,trips,theta):
    """K<=2,m=1 only: pair adjacency is the exact within-trip queue penalty."""
    if inst.stations!=1 or inst.bottles>2:raise ValueError('compact queue map requires K<=2,m=1')
    terms=[(theta,-1)]
    for tr in trips.values():
        terms.extend((idx,inst.jobs[j].weight*(e+1)) for (j,s,e),idx in tr['U'].items())
        terms.extend((idx,min(inst.jobs[i].weight,inst.jobs[j].weight)) for (i,j),idx in tr['A'].items())
    B.row(terms)
    return B

def exact_certificate(A,c,upper,lower,rhs,x,pi,deadline=math.inf):
    """Any rational row prices give a valid lower bound on this finite box.

    Rows must be equality or <=, with integer matrix coefficients. Primal
    rational reconstruction is admitted only after exact feasibility checks.
    No LP solve, correction optimization, or tolerance admission occurs here.
    """
    A=A.tocsr();c=np.asarray(c);upper=np.asarray(upper);lower=np.asarray(lower);rhs=np.asarray(rhs)
    x=np.asarray(x);pi=np.asarray(pi)
    if A.shape!=(len(rhs),len(c)) or len(x)!=len(c) or len(pi)!=len(rhs) or len(upper)!=len(c) or len(lower)!=len(rhs):raise ValueError('certificate dimensions')
    if not all(np.all(np.isfinite(a)) for a in (A.data,c,upper,rhs,x,pi)) or np.min(upper,initial=0)<0:raise ValueError('certificate finite box')
    if np.any(A.data!=np.rint(A.data)):raise ValueError('exact certificate requires integer matrix')
    eq=lower==rhs
    if np.any(~eq & ~np.isneginf(lower)):raise ValueError('only equality and upper rows supported')
    q=lambda value:F(float(value)).limit_denominator(1000000)
    costs=[F(float(v)) for v in c];uu=[F(float(v)) for v in upper];rr=[F(float(v)) for v in rhs]
    xx=[q(v) for v in x];pp=[q(v) if equality else min(F(0),q(v)) for v,equality in zip(pi,eq)]
    reduced=costs.copy();row_cost=F(0);bad=[]
    for j,(v,u) in enumerate(zip(xx,uu)):
        if not 0<=v<=u:bad.append(('box',j))
    for i in range(A.shape[0]):
        if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('exact certificate total budget')
        activity=F(0);price=pp[i]
        for t in range(A.indptr[i],A.indptr[i+1]):
            j=int(A.indices[t]);a=int(A.data[t])
            if xx[j]:activity+=a*xx[j]
            if price:reduced[j]-=a*price
        row_cost+=rr[i]*price
        if (activity!=rr[i] if eq[i] else activity>rr[i]):bad.append(('row',i))
    bound=row_cost+sum((min(F(0),v)*u for v,u in zip(reduced,uu)),F(0))
    objective=sum((a*v for a,v in zip(costs,xx)),F(0))
    if not bad and bound>objective:raise ValueError('exact certificate crossed bounds')
    return dict(dual_lower=str(bound),primal_upper=str(objective) if not bad else None,
        primal_feasible=not bad,primal_violations=bad[:10],primal_violation_count=len(bad),
        reconstruction_denominator_limit=1000000,exact=True)
