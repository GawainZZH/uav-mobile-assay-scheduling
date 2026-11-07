"""Node master overlay with shared standard triple subset-row inequalities."""
import time
import warnings
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
from .master import Matrix,_linprog,vehicle_arcs,alpha,zero_cut,route_rows as base_route_rows
from .branch_pricing import actions
from itertools import combinations
from .guarded import _dump


def route_rows(inst,r,cuts,explicit=False,branches=(),use_sri=True):
    col=base_route_rows(inst,r,cuts,explicit)
    atoms=set(actions(r))
    for i,(key,value) in enumerate(branches):
        if key[0]=='action' and key[1] in atoms:col['branch',i]=1.0
    if use_sri:
        visited={j for j,s in r.visits}
        for S in combinations(range(len(inst.jobs)),3):
            coef=len(visited.intersection(S))//2
            if coef:col[('sri',)+S]=coef
    return col


def solve_master(inst,routes,*,cuts=(),explicit=False,phase=2,integer=False,deadline=None,branches=(),use_sri=True):
    if explicit and phase != 2:
        raise ValueError('explicit control only for phase two')
    if phase not in (1,2):raise ValueError('phase')
    began=time.perf_counter()
    if phase==1:cuts=()
    M=Matrix();arcs=vehicle_arcs(inst)
    for g in range(inst.groups):
        for p in range(len(inst.pads)):
            for t in range(inst.H+1):
                rhs=int(p==0 and t==0)-int(p==0 and t==inst.H)
                M.row(('flow',g,p,t),rhs,True)
        for p in range(len(inst.pads)):
            for t in range(inst.H):
                M.row(('launch',g,p,t));M.row(('recover',g,p,t))
        for k in range(inst.H):
            M.row(('drone',g,k),inst.drones);M.row(('battery',g,k),inst.batteries)
        if not explicit:
            for a in range(inst.H):
                for b in range(a,inst.H):M.row(('hall',g,a,b))
    for j in range(len(inst.jobs)):M.row(('cover',j),1,True)
    variants=sorted({(r.g,j,s,r.recovery) for r in routes for j,s in r.visits}) if explicit else []
    if explicit:
        for v in variants:M.row(('supply',)+v,0,True)
        for g in range(inst.groups):
            for k in range(inst.H):M.row(('slot',g,k))
        M.row(('assay_cost',),0,True)
    else:
        for n in range(len(cuts)):M.row(('cost',n))
    for i,(key,value) in enumerate(branches):M.row(('branch',i),value,True)
    if use_sri:
        for S in combinations(range(len(inst.jobs)),3):M.row(('sri',)+S,1)
    zids=[]
    for ai,(g,p,t,q,u,cost) in enumerate(arcs):
        entries={('flow',g,p,t):1,('flow',g,q,u):-1,
                 ('launch',g,p,t):-inst.drones,('recover',g,p,t):-inst.drones}
        if p==q:
            if explicit:entries[('slot',g,t)]=-inst.stations
            else:
                for a in range(t+1):
                    for b in range(t,inst.H):entries[('hall',g,a,b)]=-inst.stations
                for n,beta in enumerate(cuts):entries[('cost',n)]=inst.stations*beta[g][t]
        for i,(key,value) in enumerate(branches):
            if key==('vehicle',ai):entries[('branch',i)]=1.0
        zids.append(M.column(entries,cost if phase==2 else 0,upper=1,integer=integer))
    rids=[M.column(route_rows(inst,r,cuts,explicit,branches,use_sri),r.cost(inst) if phase==2 else 0,integer=integer) for r in routes]
    theta=M.column({('assay_cost',):1} if explicit else {('cost',n):-1 for n in range(len(cuts))},int(phase==2))
    aids=[M.column({('cover',j):1},int(phase==1),upper=np.inf if phase==1 else 0) for j in range(len(inst.jobs))]
    if explicit:
        for g,j,s,e in variants:
            d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
            for k in range(e,d+1):
                M.column({('supply',g,j,s,e):1,('slot',g,k):1,('assay_cost',):-inst.jobs[j].weight*(k+1)})
    if phase==1:
        M.upper[theta]=0
        for key,iseq in zip(M.keys,M.eq):
            if iseq:aids.append(M.column({key:1},1))
            aids.append(M.column({key:-1},1))
    A=M.sparse();rhs=np.array(M.rhs);eq=np.array(M.eq,dtype=bool)
    left=60 if deadline is None else deadline-time.perf_counter()
    if left<=0:raise TimeoutError('master build exhausted budget')
    solve_start=time.perf_counter()
    if integer:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',message='Unrecognized options detected.*threads')
            result=milp(np.array(M.cost),integrality=np.array(M.integer),bounds=Bounds(np.zeros(len(M.cost)),M.upper),
                constraints=LinearConstraint(A,np.where(eq,rhs,-np.inf),rhs),options={'time_limit':left,'threads':1})
    else:
        result=_linprog(M.cost,A[eq],rhs[eq],A[~eq],rhs[~eq],list(zip(np.zeros(len(M.cost)),M.upper)),left)
    rec=dict(status=int(result.status),message=str(result.message),objective=None,seconds=time.perf_counter()-began,
             solve_seconds=time.perf_counter()-solve_start,rows=A.shape[0],columns=A.shape[1],nonzeros=A.nnz)
    rec.update(_matrix=A,_rhs=rhs,_eq=eq,_cost=np.array(M.cost),_upper=np.array(M.upper),_keys=tuple(M.keys),_x=result.x,_native_objective=result.fun,_rids=rids,_zids=zids,_theta=theta,_aids=aids)
    if result.x is None or (not integer and result.status!=0):return rec
    rec.update(objective=float(result.fun),route_values=result.x[rids],z_values=result.x[zids],theta=float(result.x[theta]),
               artificial=result.x[aids],arcs=arcs)
    q=np.zeros((inst.groups,inst.H))
    for x,(g,p,t,qpad,u,cost) in zip(rec['z_values'],arcs):
        if p==qpad:q[g,t]+=inst.stations*x
    rec['q']=q
    if not integer:
        dual={}
        for i,val in zip(np.flatnonzero(eq),result.eqlin.marginals):dual[M.keys[i]]=float(val)
        for i,val in zip(np.flatnonzero(~eq),result.ineqlin.marginals):dual[M.keys[i]]=float(val)
        rec['dual']=dual
    return rec


def solve_node(inst,routes,*,cuts=(),branches=(),phase=2,integer=False,deadline=None,use_sri=True):
    start=time.perf_counter()
    r=solve_master(inst,routes,cuts=cuts,branches=branches,phase=phase,integer=integer,deadline=deadline,use_sri=use_sri)
    A=r['_matrix'];coo=A.tocoo()
    raw={k:v for k,v in r.items() if k!='_matrix'}
    raw['_matrix_triplets']=(coo.row,coo.col,coo.data,A.shape)
    _dump('node',dict(instance=inst,routes=routes,cuts=cuts,branches=branches,phase=phase,integer=integer,use_sri=use_sri,result=raw))
    if r['_x'] is None or (not integer and r['status']!=0):return r
    x=r['_x'];eq=r['_eq'];rhs=r['_rhs'];c=r['_cost'];ub=r['_upper']
    if not np.all(np.isfinite(x)) or not np.isfinite(r['objective']):raise ValueError('nonfinite node point')
    activity=A@x
    err=max(np.max(np.abs(activity[eq]-rhs[eq]),initial=0),np.max(np.maximum(0,activity[~eq]-rhs[~eq]),initial=0),np.max(np.maximum(0,-x),initial=0),np.max(np.maximum(0,x-ub),initial=0),abs(np.dot(c,x)-r['objective']))
    if err>1e-7:raise ValueError('node primal residual')
    if not integer and r['status']==0:
        pi=np.array([r['dual'][k] for k in r['_keys']])
        if not np.all(np.isfinite(pi)):raise ValueError('nonfinite node dual')
        pi[~eq]=np.minimum(0,pi[~eq]);r['dual']={k:float(v) for k,v in zip(r['_keys'],pi)}
        rc=c-A.T@pi
        box=ub.copy();C0=2*inst.groups+len(inst.jobs)+sum(value for key,value in branches)
        box[~np.isfinite(box)]=C0 if phase==1 else 0
        box[r['_rids']]=1+C0 if phase==1 else 1
        box[r['_theta']]=0 if phase==1 else inst.H*sum(j.weight for j in inst.jobs)
        if np.max(x-box,initial=0)>1e-7:raise ValueError('point violates derived finite box')
        lower=float(np.dot(rhs,pi)+np.dot(np.minimum(0,rc),box))
        if abs(r['objective']-lower)>1e-7:raise ValueError('node finite-box dual gap')
        r['verified_dual_bound']=lower;r['route_box_upper']=1+C0 if phase==1 else 1
    r['seconds_with_validation']=time.perf_counter()-start
    return r
