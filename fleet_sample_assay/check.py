"""Independent integer schedule checks; never accepts fractional solver output."""
from .physical import validate_route


def integer_schedule(inst,routes,result,network,tol=1e-7):
    if result['status']!=0 or network['status']!=0:
        raise ValueError('optimal solver results required for this control')
    if any(abs(x-round(x))>tol for x in list(result['z_values'])+list(result['route_values'])):
        raise ValueError('fractional operation')
    flights=[r for r,x in zip(routes,result['route_values']) if x>0.5]
    paths={};wait=set();move_cost=0
    for g in range(inst.groups):
        chosen=[a for a,x in zip(result['arcs'],result['z_values']) if x>0.5 and a[0]==g]
        node=(0,0);nodes={node};used=[]
        while node!=(0,inst.H):
            nxt=[a for a in chosen if (a[1],a[2])==node]
            if len(nxt)!=1:raise AssertionError('broken vehicle path')
            a=nxt[0];used.append(a);node=a[3],a[4];nodes.add(node);move_cost+=a[5]
            if a[1]==a[3]:wait.add((g,a[2]))
        if len(used)!=len(chosen):raise AssertionError('extra vehicle components')
        paths[g]=nodes
    jobs=[];supply={};flight_cost=0
    for r in flights:
        validate_route(inst,r)
        if (r.p,r.launch) not in paths[r.g] or (r.q,r.recovery) not in paths[r.g]:
            raise AssertionError('vehicle missed flight endpoint')
        flight_cost+=inst.flight_fixed+inst.air_cost*(r.recovery-r.launch)
        for j,s in r.visits:
            jobs.append(j);supply[r.g,j,s,r.recovery]=1
    if sorted(jobs)!=list(range(len(inst.jobs))):raise AssertionError('not exact task coverage')
    for g in range(inst.groups):
        for k in range(inst.H):
            if sum(r.g==g and r.launch<=k<r.recovery for r in flights)>inst.drones:
                raise AssertionError('aircraft capacity')
            if sum(r.g==g and r.launch<=k<r.recovery+inst.charge for r in flights)>inst.batteries:
                raise AssertionError('battery capacity')
    processed={key:0.0 for key in supply};occupancy={};turnaround=0
    for (g,j,s,e,k),x in network['assignments']:
        if abs(x-round(x))>tol:raise AssertionError('fractional integer assay allocation')
        if (g,j,s,e) not in supply:raise AssertionError('unexpected sample')
        if not e<=k<=min(s+inst.jobs[j].shelf,inst.jobs[j].deadline-1):raise AssertionError('shelf/result deadline')
        if (g,k) not in wait:raise AssertionError('assay while moving')
        processed[g,j,s,e]+=x;occupancy[g,k]=occupancy.get((g,k),0)+x
        turnaround+=inst.jobs[j].weight*(k+1-s)*x
    if any(abs(x-1)>tol for x in processed.values()) or any(x>inst.stations+tol for x in occupancy.values()):
        raise AssertionError('missing assay or insufficient stations')
    true_cost=move_cost+flight_cost+turnaround
    if abs(true_cost-result['objective'])>tol:raise AssertionError('objective mismatch')
    return dict(valid=True,cost=true_cost,flights=len(flights),vehicle_paths=len(paths))
