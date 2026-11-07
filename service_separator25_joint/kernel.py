"""Branch-and-bound kernels sharing a matrix and exact transferable LP prices."""
from dataclasses import replace
from fractions import Fraction
import math
import time
import numpy as np
from review88_physical import milp as legacy
from review88_physical.physical import validate_fleet
from review88_physical.types import Route
from fleet_sample_assay_certified.master import highs,pass_lp,native_run,check
from .dual import price


def restore(m,x,start,end):
    x=np.asarray(x);plans=[]
    if np.max(abs(x-np.rint(x)),initial=0)>1e-6:raise ValueError('noninteger kernel vector')
    for g in range(m.inst.groups):
        routes=[r for i,r in enumerate(m.routes) if x[m.rids[g,i]]>.5]
        arcs=[a for i,a in enumerate(m.arcs) if x[m.zids[g,i]]>.5]
        trail=[];node=(start[g],0)
        while node!=(end[g],m.inst.H):
            selected=[a for a in arcs if a[:2]==node]
            if len(selected)!=1:raise ValueError('kernel path continuity')
            a=selected[0];trail.append(a);node=a[2:]
            if len(trail)>len(arcs):raise ValueError('kernel path cycle')
        if len(trail)!=len(arcs):raise ValueError('extra kernel path')
        assays=[(j,s,e,k) for (gg,j,s,e,k),col in m.assignments.items() if gg==g and x[col]>.5]
        if m.policy=='QC' and any(p!=b and any(e<=t<k+m.inst.assay_duration for j,s,e,k in assays) for p,t,b,u in trail):
            raise ValueError('QC kernel departure')
        plans.append(dict(routes=routes,trail=trail,assays=assays))
    physical=validate_boundary(m.inst,plans,start,end)
    if not physical['valid']:raise ValueError(physical['errors'])
    if abs(float(physical['objective'])-np.array(m.M.cost)@x)>1e-6:raise ValueError('kernel cost')
    return dict(group_plans=plans,check=physical)


def validate_boundary(inst,plans,start,end):
    """Add real depot access outside the kernel; reuse the independent validator."""
    before=max(inst.road(0,p) for p in start);after=max(inst.road(p,0) for p in end)
    jobs=tuple(replace(j,lo=j.lo+before,hi=j.hi+before,deadline=j.deadline+before) for j in inst.jobs)
    extended=replace(inst,H=before+inst.H+after,jobs=jobs)
    concrete=[];extra=Fraction(0)
    for g,plan in enumerate(plans):
        s,e=start[g],end[g];a=inst.road(0,s);b=inst.road(e,0)
        trail=[]
        if s!=0:trail.append((0,0,s,a))
        trail.extend((s,t,s,t+1) for t in range(a,before))
        trail.extend((p,t+before,q,u+before) for p,t,q,u in plan['trail'])
        finish=before+inst.H
        if e!=0:trail.append((e,finish,0,finish+b))
        trail.extend((0,t,0,t+1) for t in range(finish+b,extended.H))
        routes=[replace(r,launch=r.launch+before,recovery=r.recovery+before,sample=r.sample+before)
                if not isinstance(r,dict) else Route(**dict(r,launch=r['launch']+before,
                    recovery=r['recovery']+before,sample=r['sample']+before)) for r in plan['routes']]
        assays=[(j,s+before,e+before,k+before) for j,s,e,k in plan['assays']]
        concrete.append(dict(routes=routes,assays=assays,trail=trail));extra+=inst.road_cost*(a+b)
    result=validate_fleet(extended,concrete)
    result['objective']-=extra;result['road_cost']-=extra
    result['validation_access_cost']=extra
    return result


class Kernel:
    def __init__(self,inst,deadline,policy='J',dual_transfer=True):
        began=time.perf_counter();self.inst=inst;self.policy=policy
        from .model import build
        self.m=m=build(inst,policy,deadline)
        for g in range(inst.groups):m.M.eq[m.M.index['pack',g]]=True
        self.cost=np.array(m.M.cost,dtype=float);self.upper=np.array(m.M.upper,dtype=float)
        self.eq=np.array(m.M.eq,dtype=bool);self.base=np.array(m.M.rhs,dtype=float)
        self.prices=[];self.solved={};self.events=[];self.dual_transfer=dual_transfer
        self.projected=[];self.bound_cache={}
        self.lp=highs();pass_lp(self.lp,m.A,self.cost,self.upper,np.where(self.eq,self.base,-np.inf),self.base)
        from .scip_backend import Proxy
        self.mip=Proxy(m.A,self.cost,self.upper,np.where(self.eq,self.base,-np.inf),self.base,np.array(m.M.integer))
        self.build_seconds=time.perf_counter()-began;self.native_calls=0
        self.elementary=sum(min(inst.flight_fixed+inst.air_cost*(r.recovery-r.launch)+
            inst.jobs[j].weight*(r.recovery+inst.handoff_delay+inst.assay_duration-r.sample)
            for r in m.routes if r.job==j) for j in range(len(inst.jobs)))

    def rhs(self,start,end,consumed):
        rhs=self.base.copy()
        for g in range(self.inst.groups):
            for p in range(len(self.inst.pads)):
                rhs[self.m.M.index['flow',g,p,0]]=int(p==start[g])
                rhs[self.m.M.index['flow',g,p,self.inst.H]]=-int(p==end[g])
            rhs[self.m.M.index['pack',g]]=consumed[g]
        return rhs

    def lower(self,start,end,consumed):
        key=(tuple(start),tuple(end),tuple(consumed))
        if key in self.solved:
            r=self.solved[key]
            return math.inf if r['status']=='infeasible' else r.get('objective',0) if r['optimal'] else r['lower_bound']
        road=sum(self.inst.road_cost*self.inst.road(a,b) for a,b in zip(start,end))
        if key in self.bound_cache and self.bound_cache[key][0]==len(self.prices):return self.bound_cache[key][1]
        # Project once onto the only RHS coordinates changed by a boundary.
        while len(self.projected)<len(self.prices):
            p=self.prices[len(self.projected)];constant=p.constant
            variable=set()
            for g in range(self.inst.groups):
                variable.add(self.m.M.index['pack',g])
                variable.update(self.m.M.index['flow',g,pad,t] for pad in range(len(self.inst.pads)) for t in (0,self.inst.H))
            constant+=sum(int(y)*int(b) for i,(y,b) in enumerate(zip(p.weights,self.base)) if i not in variable)
            self.projected.append((p,constant))
        best=int(self.elementary+road)
        for p,constant in self.projected:
            numerator=constant
            for g in range(self.inst.groups):
                numerator+=int(p.weights[self.m.M.index['flow',g,start[g],0]])
                numerator-=int(p.weights[self.m.M.index['flow',g,end[g],self.inst.H]])
                numerator+=int(p.weights[self.m.M.index['pack',g]])*consumed[g]
            best=max(best,-(-numerator//p.scale))
        self.bound_cache[key]=(len(self.prices),best)
        return best

    def solve(self,start,end,consumed,deadline):
        key=(tuple(start),tuple(end),tuple(consumed))
        if key in self.solved and (self.solved[key]['optimal'] or self.solved[key]['status']=='infeasible'):
            return self.solved[key]
        began=time.perf_counter();rhs=self.rhs(*key);ids=np.arange(len(rhs),dtype=np.int32)
        lower=np.where(self.eq,rhs,-np.inf)
        if self.dual_transfer:
            check(self.lp.changeRowsBounds(len(rhs),ids,lower,rhs));status,sol,info,msg=native_run(self.lp,deadline);self.native_calls+=1
            if status==0 and sol.dual_valid:
                p=price(self.m.A,self.cost,self.upper,self.eq,sol.row_dual);self.prices.append(p)
            elif status==2:
                out=dict(status='infeasible',optimal=False,objective=None,lower_bound=None,fleet=None,context=key,seconds=time.perf_counter()-began)
                self.solved[key]=out;self.events.append(out);return out
        check(self.mip.changeRowsBounds(len(rhs),ids,lower,rhs));status,sol,info,msg=native_run(self.mip,deadline);self.native_calls+=1
        value=None;fleet=None
        if sol.value_valid:
            from .model import complete
            vector=complete(self.m,np.asarray(sol.col_value),start,end);a=self.m.A@vector
            err=max(np.max(abs(a[self.eq]-rhs[self.eq]),initial=0),np.max(a[~self.eq]-rhs[~self.eq],initial=0),
                np.max(-vector,initial=0),np.max(vector-self.upper,initial=0),np.max(abs(vector-np.rint(vector)),initial=0))
            if err>1e-6:raise ValueError(('kernel residual',err))
            fleet=restore(self.m,vector,start,end);value=int(fleet['check']['objective'])
        bound=float(info.mip_dual_bound)
        closed=status==0 and value is not None and abs(bound-value)<=1e-6
        out=dict(status='optimal' if closed else 'infeasible' if status==2 else 'unresolved',
            optimal=closed,objective=value,lower_bound=math.ceil(bound-1e-6) if math.isfinite(bound) else 0,
            fleet=fleet,context=key,seconds=time.perf_counter()-began,nodes=int(info.mip_node_count),
            vector=None if not sol.value_valid else vector.tolist(),
            native_vector=None if not sol.value_valid else list(sol.col_value))
        self.solved[key]=out;self.events.append(out);return out
