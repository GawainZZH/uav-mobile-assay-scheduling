"""Finite LP prices and positive-support q-family certification."""
import time
from fractions import Fraction as F
import numpy as np
from .matrix import sparse,activity,band_rows
from .bands import exact,update


def check(m,x,*,integer=False,tolerance=1e-7):
    if len(x)!=len(m.cost)or not np.isfinite(np.asarray(x,dtype=float)).all():raise ValueError('finite master vector required')
    err=max([0.]+[-float(v)for v in x]+[float(v-u)for v,u in zip(x,m.upper)])
    for r in m.rows:
        a=activity(r,x);err=max(err,abs(a)if r['eq']else a)
    if integer:err=max(err,max((abs(float(x[j])-round(float(x[j])))for j,v in enumerate(m.integer)if v),default=0))
    if err>tolerance:raise ValueError(('stock-flow primal residual',err))
    return err


def rational_certificate(m,x,dual):
    y=[F(str(v)).limit_denominator(10**6)if r['eq']else min(F(0),F(str(v)).limit_denominator(10**6))for r,v in zip(m.rows,dual)]
    reduced=list(m.cost)
    for r,v in zip(m.rows,y):
        for j,c in r['entries'].items():reduced[j]-=v*c
    constant=sum((u*min(F(0),v)for u,v in zip(m.upper,reduced)),F(0))
    lower=constant+sum((v*r['rhs']for r,v in zip(m.rows,y)),F(0))
    primal=[F(0)if abs(float(v))<1e-9 else F(str(v)).limit_denominator(10**6)for v in x]
    valid=all(0<=v<=u for v,u in zip(primal,m.upper))
    for r in m.rows:
        a=sum((c*primal[j]for j,c in r['entries'].items()),F(0))-r['rhs']
        if a!=0 if r['eq']else a>0:valid=False
    value=sum((c*v for c,v in zip(m.cost,primal)),F(0))if valid else None
    if value is not None and lower>value:raise ValueError('rational LP weak-duality contradiction')
    return dict(exact_lower=lower,exact_primal=primal if valid else None,exact_primal_objective=value,
        exact_optimal=valid and value==lower,price=dict(weights=y,box_constant=constant))


def solve_lp(m,deadline):
    from fleet_sample_assay_certified.master import highs,pass_lp,native_run
    h=highs();started=time.perf_counter()
    try:
        rhs=np.asarray([float(r['rhs'])for r in m.rows]);eq=np.asarray([r['eq']for r in m.rows])
        pass_lp(h,sparse(m),np.asarray(m.cost,dtype=float),np.asarray(m.upper,dtype=float),np.where(eq,rhs,-np.inf),rhs)
        code,sol,info,message=native_run(h,deadline)
        result=dict(status=code,message=message,vector=None,raw_objective=None,row_count=len(m.rows),column_count=len(m.cost))
        if code==0 and sol.value_valid and sol.dual_valid:
            x=list(map(float,sol.col_value));residual=check(m,x)
            result.update(vector=x,raw_objective=sum(float(c)*v for c,v in zip(m.cost,x)),residual=residual,
                **rational_certificate(m,x,sol.row_dual))
        result['seconds']=time.perf_counter()-started
        return result
    finally:h.clear()


def support(m,x):
    return {i:sum(x[j]for j in js)for i,js in m.family_columns.items()if sum(x[j]for j in js)>0}


def close_root(m,query,validate,deadline):
    events=[];queried=set();iterations=[];started=time.perf_counter()
    while time.perf_counter()<deadline:
        lp=solve_lp(m,deadline);iterations.append(lp)
        if lp['status']!=0:
            return dict(status='infeasible'if lp['status']==2 else'unresolved',lp=lp,events=events,iterations=iterations,positive_support_exact=False)
        # Strict closure uses an exactly feasible primal, not a numerical support cutoff.
        x=lp['exact_primal']if lp['exact_primal']is not None else lp['vector']
        selected=support(m,x);unknown=[i for i in selected if not exact(m.bands[i])]
        if not unknown:
            return dict(status='closed'if lp['exact_optimal']else'numerical_support_closed',lp=lp,events=events,iterations=iterations,
                positive_support_exact=True,seconds=time.perf_counter()-started,
                all_support_proof_scopes={str(i):m.bands[i].get('proof_scope','exact_road')for i in selected})
        i=next((i for i in sorted(unknown,key=lambda i:(-float(selected[i]),i))if i not in queried),None)
        if i is None:return dict(status='unresolved',reason='positive family already queried without exact closure',lp=lp,events=events,iterations=iterations,positive_support_exact=False)
        queried.add(i);old=dict(m.bands[i]);r=query(i);update(m.bands[i],r,lambda p:validate(i,p))
        events.append(dict(family=i,record=r))
        if m.bands[i]['infeasible']or old['lower']!=m.bands[i]['lower']or old['upper']!=m.bands[i]['upper']:
            m.rows.extend(band_rows(m,i));continue
        if not exact(m.bands[i]):
            return dict(status='unresolved',reason='queried family yielded no interval progress',lp=lp,events=events,iterations=iterations,positive_support_exact=False)
    return dict(status='unresolved',reason='root deadline',events=events,iterations=iterations,positive_support_exact=False)
