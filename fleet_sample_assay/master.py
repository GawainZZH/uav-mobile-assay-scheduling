"""Same-row explicit assay LP and exact projected route master for tiny controls."""
import time
import warnings
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import coo_matrix
from .physical import vehicle_arcs


def zero_cut(inst):
    return tuple(tuple(0.0 for k in range(inst.H)) for g in range(inst.groups))


def alpha(inst,beta,g,j,s,e):
    d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
    if e>d:
        raise ValueError('empty assay interval')
    return min(inst.jobs[j].weight*(k+1)-beta[g][k] for k in range(e,d+1))


def route_rows(inst,r,cuts,explicit=False):
    out={('cover',j):1.0 for j,s in r.visits}
    out[('launch',r.g,r.p,r.launch)]=1.0
    out[('recover',r.g,r.q,r.recovery)]=1.0
    for k in r.drone_slots():out[('drone',r.g,k)]=1.0
    for k in r.battery_slots(inst):out[('battery',r.g,k)]=1.0
    for j,s in r.visits:
        if explicit:
            out[('supply',r.g,j,s,r.recovery)]=-1.0
        else:
            lo,hi=r.window(inst,j,s)
            for a in range(lo+1):
                for b in range(hi,inst.H):
                    key=('hall',r.g,a,b)
                    out[key]=out.get(key,0)+1.0
            for n,beta in enumerate(cuts):
                key=('cost',n)
                out[key]=out.get(key,0)+alpha(inst,beta,r.g,j,s,r.recovery)
    return out


class Matrix:
    def __init__(self):
        self.keys=[];self.index={};self.rhs=[];self.eq=[]
        self.columns=[];self.cost=[];self.upper=[];self.integer=[]

    def row(self,key,rhs=0.0,eq=False):
        if key in self.index:raise ValueError('duplicate row')
        self.index[key]=len(self.keys);self.keys.append(key);self.rhs.append(rhs);self.eq.append(eq)

    def column(self,entries,cost=0.0,upper=np.inf,integer=False):
        self.columns.append(entries);self.cost.append(cost);self.upper.append(upper);self.integer.append(int(integer))
        return len(self.cost)-1

    def sparse(self):
        rows=[];cols=[];data=[]
        for i,col in enumerate(self.columns):
            for key,value in col.items():
                if value:
                    rows.append(self.index[key]);cols.append(i);data.append(value)
        return coo_matrix((data,(rows,cols)),shape=(len(self.keys),len(self.columns))).tocsr()


def _linprog(c,Aeq,beq,Aub,bub,bounds,seconds):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected.*threads')
        return linprog(c,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=bounds,
            method='highs',options={'time_limit':max(0.01,seconds),'threads':1})


def solve_master(inst,routes,*,cuts=(),explicit=False,phase=2,integer=False,deadline=None):
    if explicit and phase != 2:
        raise ValueError('explicit control only for phase two')
    if phase not in (1,2):raise ValueError('phase')
    began=time.perf_counter();M=Matrix();arcs=vehicle_arcs(inst)
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
    zids=[]
    for g,p,t,q,u,cost in arcs:
        entries={('flow',g,p,t):1,('flow',g,q,u):-1,
                 ('launch',g,p,t):-inst.drones,('recover',g,p,t):-inst.drones}
        if p==q:
            if explicit:entries[('slot',g,t)]=-inst.stations
            else:
                for a in range(t+1):
                    for b in range(t,inst.H):entries[('hall',g,a,b)]=-inst.stations
                for n,beta in enumerate(cuts):entries[('cost',n)]=inst.stations*beta[g][t]
        zids.append(M.column(entries,cost if phase==2 else 0,upper=1,integer=integer))
    rids=[M.column(route_rows(inst,r,cuts,explicit),r.cost(inst) if phase==2 else 0,integer=integer) for r in routes]
    theta=M.column({('assay_cost',):1} if explicit else {('cost',n):-1 for n in range(len(cuts))},int(phase==2))
    aids=[M.column({('cover',j):1},int(phase==1),upper=np.inf if phase==1 else 0) for j in range(len(inst.jobs))]
    if explicit:
        for g,j,s,e in variants:
            d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
            for k in range(e,d+1):
                M.column({('supply',g,j,s,e):1,('slot',g,k):1,('assay_cost',):-inst.jobs[j].weight*(k+1)})
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
    if result.status!=0:return rec
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


def assay_network(inst,routes,values,q,deadline=None):
    began=time.perf_counter();supply={}
    for r,x in zip(routes,values):
        if x < -1e-7:raise ValueError('negative route supply')
        if x<=1e-10:continue
        for j,s in r.visits:
            key=r.g,j,s,r.recovery
            supply[key]=supply.get(key,0)+max(0.0,float(x))
    if not supply:
        return dict(status=0,cost=0.0,beta=zero_cut(inst),seconds=time.perf_counter()-began)
    variants=sorted(supply);edges=[];c=[];row=[];col=[];caprow=[]
    for i,(g,j,s,e) in enumerate(variants):
        d=min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)
        for k in range(e,d+1):
            idx=len(edges);edges.append((g,j,s,e,k));c.append(inst.jobs[j].weight*(k+1))
            row.append(i);col.append(idx);caprow.append(g*inst.H+k)
    if not edges:return dict(status=2,cost=None,beta=None,seconds=time.perf_counter()-began)
    Aeq=coo_matrix((np.ones(len(edges)),(row,col)),shape=(len(variants),len(edges))).tocsr()
    Aub=coo_matrix((np.ones(len(edges)),(caprow,col)),shape=(inst.groups*inst.H,len(edges))).tocsr()
    left=60 if deadline is None else deadline-time.perf_counter()
    if left<=0:raise TimeoutError('network budget')
    result=_linprog(c,Aeq,[supply[v] for v in variants],Aub,np.maximum(0,q).ravel(),(0,None),left)
    rec=dict(status=int(result.status),cost=None,beta=None,seconds=time.perf_counter()-began)
    if result.status==0:
        beta=np.array(result.ineqlin.marginals).reshape((inst.groups,inst.H))
        if np.max(beta)>1e-7:raise AssertionError('capacity dual sign')
        rec.update(cost=float(result.fun),beta=tuple(tuple(min(0,float(b)) for b in g) for g in beta),
                   assignments=[(edge,float(x)) for edge,x in zip(edges,result.x) if x>1e-9])
    return rec
