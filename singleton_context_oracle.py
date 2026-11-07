"""Exact one-task-per-active-group J oracle, with constructive certificates.

No LP/MIP, group-state graph or numerical assignment solver is used. All costs
and assignment potentials are Fractions. Other responsibility counts are refused.
"""
from fractions import Fraction as F
import time


def clock(deadline):
    if time.perf_counter()>=deadline:raise TimeoutError('singleton context allowance')


def exact_assignment(costs):
    n=len(costs)
    if any(len(row)!=n for row in costs):raise ValueError('square assignment required')
    if any(v is not None and(type(v)is bool or not isinstance(v,(int,F)))for row in costs for v in row):
        raise ValueError('explicit rational assignment costs required')
    c=[[None if v is None else F(v)for v in row]for row in costs]
    u=[F(0)]*(n+1);v=[F(0)]*(n+1);p=[0]*(n+1);way=[0]*(n+1)
    for row in range(1,n+1):
        p[0]=row;j0=0;minimum=[None]*(n+1);used=[False]*(n+1)
        while True:
            used[j0]=True;i0=p[j0];delta=None;j1=None
            for j in range(1,n+1):
                if used[j]:continue
                value=c[i0-1][j-1]
                if value is not None:
                    rc=value-u[i0]-v[j]
                    if minimum[j]is None or rc<minimum[j]:minimum[j]=rc;way[j]=j0
                if minimum[j]is not None and(delta is None or minimum[j]<delta):
                    delta=minimum[j];j1=j
            if delta is None:
                rows=tuple(sorted({p[j]-1 for j in range(n+1)if used[j]}))
                neighbors=tuple(j for j in range(n)if any(c[i][j]is not None for i in rows))
                if len(neighbors)>=len(rows):raise ValueError('invalid Hall deficiency')
                return dict(status='infeasible',objective=None,hall_rows=rows,hall_neighbors=neighbors)
            for j in range(n+1):
                if used[j]:u[p[j]]+=delta;v[j]-=delta
                elif minimum[j]is not None:minimum[j]-=delta
            j0=j1
            if p[j0]==0:break
        while j0:
            j1=way[j0];p[j0]=p[j1];j0=j1
    assignment=[None]*n
    for j in range(1,n+1):assignment[p[j]-1]=j-1
    value=sum((c[i][j]for i,j in enumerate(assignment)),F(0))
    if value!=sum(u[1:])+sum(v[1:]):raise ValueError('assignment primal/dual mismatch')
    if any(c[i][j]is not None and u[i+1]+v[j+1]>c[i][j]for i in range(n)for j in range(n)):
        raise ValueError('assignment dual inequality violated')
    return dict(status='optimal',objective=value,row_to_column=tuple(assignment),
                row_dual=tuple(u[1:]),column_dual=tuple(v[1:]))


def road_data(inst):
    P=len(inst.pads);d=[[inst.road(a,b)for b in range(P)]for a in range(P)]
    paths=[[(a,)if a==b else(a,b)for b in range(P)]for a in range(P)]
    for k in range(P):
        for a in range(P):
            for b in range(P):
                candidate=d[a][k]+d[k][b]
                if candidate<d[a][b]:d[a][b]=candidate;paths[a][b]=paths[a][k]+paths[k][b][1:]
    return dict(distance=tuple(tuple(row)for row in d),paths=tuple(tuple(row)for row in paths))


def _trail(inst,roads,start,end,*,route=None,assay_pad=None,assay_start=None):
    trail=[];position=start;t=0
    def wait(until):
        nonlocal t
        if until<t:raise ValueError('negative carrier waiting interval')
        while t<until:trail.append((position,t,position,t+1));t+=1
    def move(target):
        nonlocal position,t
        path=roads['paths'][position][target]
        for a,b in zip(path,path[1:]):
            finish=t+inst.road(a,b);trail.append((a,t,b,finish));position=b;t=finish
    if route is not None:
        move(route.p);wait(route.launch);move(route.q);wait(route.recovery)
        if assay_pad==route.q:
            wait(max(route.recovery+inst.ground_duration,assay_start+inst.assay_duration))
        else:
            wait(route.recovery+inst.ground_duration);move(assay_pad)
            wait(assay_start+inst.assay_duration)
    move(end);wait(inst.H)
    return trail


def _single(inst,job,start,end,routes,roads,deadline):
    d=roads['distance'];Q=inst.assay_duration;U=inst.ground_duration;h=inst.handoff_delay
    best=None;witness=None;count=0;road_feasible=0;assay_candidates=0;j=inst.jobs[job]
    for route in routes:
        clock(deadline)
        if route.job!=job:continue
        count+=1;p,l,q,rho,s=route.p,route.launch,route.q,route.recovery,route.sample
        if d[start][p]>l or d[p][q]>rho-l or rho+U>inst.H:continue
        road_feasible+=1
        for c in range(len(inst.pads)):
            assay_candidates+=1
            k=rho+h if c==q else max(rho+h,rho+U+d[q][c])
            finish=(max(rho+U,k+Q)+d[q][end])if c==q else k+Q+d[c][end]
            if finish>inst.H or k>s+j.shelf or k+Q>min(j.deadline,inst.H):continue
            value=inst.flight_fixed+inst.air_cost*(rho-l)+j.weight*(k+Q-s)+inst.road_cost*(
                d[start][p]+d[p][q]+d[q][c]+d[c][end])
            key=(value,p,l,q,rho,s,c,k)
            if best is None or key<best:best=key;witness=(route,c,k)
    inventory=dict(routes_considered=count,road_feasible_routes=road_feasible,assay_candidates_checked=assay_candidates)
    if best is None:return dict(status='infeasible',objective=None,job=job,start=start,end=end,inventory=inventory)
    route,c,k=witness
    plan=dict(routes=[route],assays=[(job,route.sample,route.recovery,k)],
              trail=_trail(inst,roads,start,end,route=route,assay_pad=c,assay_start=k))
    return dict(status='optimal',objective=F(best[0]),job=job,start=start,end=end,
                assay_pad=c,plan=plan,inventory=inventory)


def single_task(inst,job,start,end,*,deadline=float('inf')):
    if type(job)is not int or not 0<=job<len(inst.jobs):raise ValueError('job outside domain')
    if any(type(p)is not int or not 0<=p<len(inst.pads)for p in (start,end)):raise ValueError('pad outside domain')
    from review88_physical.milp import independent_routes
    return _single(inst,job,start,end,independent_routes(inst,deadline),road_data(inst),deadline)


def solve_context(inst,start,end,consumed,*,deadline=float('inf')):
    began=time.perf_counter();start=tuple(start);end=tuple(end);consumed=tuple(consumed)
    if any(len(v)!=inst.groups for v in (start,end,consumed)):raise ValueError('complete group labels required')
    if any(type(v)is not int or v not in (0,1)for v in consumed):raise ValueError('singleton counts0/1 only')
    if sum(consumed)!=len(inst.jobs):raise ValueError('task count differs from fixed responsibilities')
    if any(type(p)is not int or not 0<=p<len(inst.pads)for p in start+end):raise ValueError('pad outside domain')
    from review88_physical.milp import independent_routes
    from service_separator25_joint.kernel import validate_boundary
    clock(deadline);roads=road_data(inst);routes=independent_routes(inst,deadline)
    active=tuple(g for g,v in enumerate(consumed)if v);idle=tuple(g for g,v in enumerate(consumed)if not v)
    idle_set=set(idle);by_job=[[]for _ in inst.jobs]
    for route in routes:by_job[route.job].append(route)
    for g in idle:
        if roads['distance'][start[g]][end[g]]>inst.H:
            return dict(status='infeasible',objective=None,idle_group=g,shortest_time=roads['distance'][start[g]][end[g]])
    cache={};costs=[]
    for j in range(len(inst.jobs)):
        row=[]
        for g in active:
            key=(j,start[g],end[g])
            if key not in cache:cache[key]=_single(inst,*key,by_job[j],roads,deadline)
            row.append(cache[key]['objective'])
        costs.append(row)
    clock(deadline);assignment=exact_assignment(costs)
    base=dict(context=(start,end,consumed),active_groups=active,assignment=assignment,
        cost_matrix=costs,cost_table=tuple(cache.values()),roads=roads,enumerated_routes=routes,
        native_optimizer_calls=0,resource_graph_builds=0)
    if assignment['status']=='infeasible':return dict(base,status='infeasible',objective=None)
    plans=[dict(routes=[],assays=[],trail=_trail(inst,roads,start[g],end[g]))if g in idle_set else None
           for g in range(inst.groups)]
    for j,col in enumerate(assignment['row_to_column']):
        g=active[col];plans[g]=cache[j,start[g],end[g]]['plan']
    idle_cost=sum((inst.road_cost*roads['distance'][start[g]][end[g]]for g in idle),F(0))
    objective=assignment['objective']+idle_cost;check=validate_boundary(inst,plans,start,end)
    if not check['valid']or check['objective']!=objective:raise ValueError(('singleton physical certificate mismatch',check))
    clock(deadline)
    return dict(base,status='optimal',optimal=True,objective=objective,lower_bound=objective,
        idle_cost=idle_cost,fleet=dict(group_plans=plans,check=check),seconds=time.perf_counter()-began)
