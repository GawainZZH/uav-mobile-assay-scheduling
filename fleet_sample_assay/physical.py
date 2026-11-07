"""Bounded-clock physical model and endpoint-conditioned exact label pricing."""
from dataclasses import dataclass
from math import ceil
import time


@dataclass(frozen=True)
class Job:
    xy: tuple
    lo: int
    hi: int
    shelf: int
    deadline: int
    weight: int = 1
    volume: int = 1


@dataclass(frozen=True)
class Instance:
    name: str
    H: int
    pads: tuple
    jobs: tuple
    groups: int = 1
    drones: int = 1
    batteries: int = 2
    charge: int = 2
    bottles: int = 2
    payload: int = 2
    endurance: int = 7
    stations: int = 1
    flight_fixed: int = 2
    air_cost: int = 1
    road_cost: int = 1

    def __post_init__(self):
        if any(type(x) is not int or x < 1 for x in (self.H,self.groups,self.drones,self.batteries,self.charge,self.bottles,self.payload,self.endurance,self.stations)):
            raise ValueError('positive integer clocks and capacities required')
        if len(set(self.pads)) != len(self.pads) or not self.pads or not self.jobs:
            raise ValueError('distinct pads and nonempty tasks required')
        if min(self.flight_fixed,self.air_cost,self.road_cost) < 0:
            raise ValueError('nonnegative operating costs required')
        for j in self.jobs:
            if not 0 <= j.lo <= j.hi < self.H or j.shelf < 0 or not 1 <= j.deadline <= self.H or min(j.volume,j.weight) < 1:
                raise ValueError('invalid job data')

    def xy(self,node):
        return self.pads[node] if node < len(self.pads) else self.jobs[node-len(self.pads)].xy

    def fly(self,a,b):
        return max(1,ceil(sum(abs(x-y) for x,y in zip(self.xy(a),self.xy(b)))/2))

    def road(self,p,q):
        return max(1,sum(abs(x-y) for x,y in zip(self.pads[p],self.pads[q])))


@dataclass(frozen=True)
class Route:
    g: int
    p: int
    launch: int
    q: int
    recovery: int
    visits: tuple

    @property
    def context(self):
        return self.g,self.p,self.launch,self.q,self.recovery

    def drone_slots(self):
        return tuple(range(self.launch,self.recovery))

    def battery_slots(self,inst):
        return tuple(range(self.launch,min(inst.H,self.recovery+inst.charge)))

    def window(self,inst,j,s):
        return self.recovery,min(inst.H-1,s+inst.jobs[j].shelf,inst.jobs[j].deadline-1)

    def cost(self,inst):
        return inst.flight_fixed+inst.air_cost*(self.recovery-self.launch)-sum(inst.jobs[j].weight*s for j,s in self.visits)


def case(name):
    jobs=(Job((1,0),1,9,6,12,1),Job((1,1),1,9,6,12,2),Job((2,0),1,9,6,12,1))
    kw=dict(name=name,H=12,pads=((0,0),(4,0)),jobs=jobs)
    if name == 'single_loose':
        return Instance(**kw)
    if name == 'single_tight':
        return Instance(**kw,batteries=1)
    if name == 'two_groups':
        return Instance(**{**kw,'jobs':jobs+(Job((3,1),2,8,5,12,2),)},groups=2,drones=1,batteries=1)
    if name == 'shelf_infeasible':
        return Instance(**{**kw,'jobs':(Job((2,2),2,6,0,12),)})
    raise ValueError(name)


def vehicle_arcs(inst):
    """(g,p,t,q,u,cost), including waiting arcs; no disconnected horizon nodes."""
    out=[]
    for g in range(inst.groups):
        for p in range(len(inst.pads)):
            for t in range(inst.H):
                for q in range(len(inst.pads)):
                    dt=1 if p == q else inst.road(p,q)
                    if t+dt <= inst.H:
                        out.append((g,p,t,q,t+dt,0 if p == q else dt*inst.road_cost))
    return tuple(out)


def contexts(inst):
    for g in range(inst.groups):
        for p in range(len(inst.pads)):
            for launch in range(inst.H):
                for q in range(len(inst.pads)):
                    for recovery in range(launch+1,min(inst.H-1,launch+inst.endurance)+1):
                        yield g,p,launch,q,recovery


def validate_route(inst,r):
    if r.context not in contexts(inst) or not r.visits:
        raise ValueError('invalid endpoints or empty flight')
    if len(r.visits)>inst.bottles or len({j for j,s in r.visits}) != len(r.visits):
        raise ValueError('bottle/elementarity violation')
    if sum(inst.jobs[j].volume for j,s in r.visits)>inst.payload:
        raise ValueError('payload violation')
    last=r.p;now=r.launch;m=len(inst.pads)
    for j,s in r.visits:
        job=inst.jobs[j]
        if not job.lo <= s <= job.hi or now+inst.fly(last,m+j)>s:
            raise ValueError('travel/window violation')
        a,b=r.window(inst,j,s)
        if a>b:
            raise ValueError('empty assay window')
        last=m+j;now=s+1
    if now+inst.fly(last,r.q)>r.recovery:
        raise ValueError('return violation')
    return True


def enumerate_routes(inst,deadline=None):
    """Independent exhaustive DFS: truth control only, never performance baseline."""
    out=[];m=len(inst.pads)
    for ctx in contexts(inst):
        if deadline is not None and time.perf_counter()>deadline:
            raise TimeoutError('route enumeration budget')
        g,p,l,q,e=ctx
        def dfs(last,now,visits,used,volume):
            if visits and now+inst.fly(last,q)<=e:
                out.append(Route(*ctx,visits))
            if len(visits)>=inst.bottles:
                return
            for j,job in enumerate(inst.jobs):
                if j in used or volume+job.volume>inst.payload:
                    continue
                lo=max(job.lo,now+inst.fly(last,m+j),e-job.shelf)
                hi=min(job.hi,e-1-inst.fly(m+j,q))
                if e>min(inst.H-1,job.deadline-1):
                    continue
                for s in range(lo,hi+1):
                    dfs(m+j,s+1,visits+((j,s),),used|{j},volume+job.volume)
        dfs(p,l,(),set(),0)
    return tuple(out)


def price(inst,common,visit,deadline=None):
    """All endpoint conditions, exact state dominance, visit-additive objective.

    common(ctx) includes endpoint/occupancy costs; visit(ctx,j,s) includes all
    cover/Hall/assay terms. No queue history; no assumption of nonnegative costs.
    """
    best=float('inf');winner=None;labels=0;conditions=0;m=len(inst.pads)
    for ctx in contexts(inst):
        if deadline is not None and time.perf_counter()>deadline:
            raise TimeoutError('exact pricing budget')
        conditions+=1;g,p,l,q,e=ctx
        # Layer by visited count guarantees every predecessor is finalized.
        layer={(0,p,l):(common(ctx),())}
        for count in range(inst.bottles+1):
            following={}
            for (mask,last,now),(value,path) in layer.items():
                labels+=1
                if mask and now+inst.fly(last,q)<=e and value<best:
                    best=value;winner=Route(*ctx,path)
                if count == inst.bottles:
                    continue
                volume=sum(job.volume for j,job in enumerate(inst.jobs) if mask>>j&1)
                for j,job in enumerate(inst.jobs):
                    if mask>>j&1 or volume+job.volume>inst.payload:
                        continue
                    low=max(job.lo,now+inst.fly(last,m+j),e-job.shelf)
                    high=min(job.hi,e-1-inst.fly(m+j,q))
                    if e>min(inst.H-1,job.deadline-1):
                        continue
                    for s in range(low,high+1):
                        state=mask|1<<j,m+j,s+1
                        val=value+visit(ctx,j,s)
                        if state not in following or val<following[state][0]:
                            following[state]=val,path+((j,s),)
            layer=following
    return dict(value=best,route=winner,labels=labels,conditions=conditions)
