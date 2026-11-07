from review80_fairness.preprocess import shortest_road


def canonical(pairs):
    order=tuple(sorted(range(len(pairs)),key=lambda i:pairs[i]))
    return tuple(pairs[i] for i in order),order


def envelope(inst,j):
    job=inst.jobs[j];P=len(inst.pads)
    earliest=job.lo+inst.sampling_duration+min(inst.fly(P+j,p) for p in range(P))-inst.endurance
    latest_assay=min(inst.H-inst.assay_duration,job.hi+job.shelf,job.deadline-inst.assay_duration)
    end=max(job.deadline,latest_assay-inst.handoff_delay+inst.ground_duration)
    return max(0,earliest),min(inst.H,end)


def certify(inst,boundaries):
    if boundaries[0]!=0 or boundaries[-1]!=inst.H or any(a>=b for a,b in zip(boundaries,boundaries[1:])):
        raise ValueError('ordered complete horizon partition required')
    dist=shortest_road(inst);diameter=max(max(row) for row in dist);cuts=[]
    for c in boundaries[1:-1]:
        left=[j for j,x in enumerate(inst.jobs) if x.lo<c]
        right=[j for j,x in enumerate(inst.jobs) if x.lo>=c]
        a=max((envelope(inst,j)[1] for j in left),default=0)
        b=min((envelope(inst,j)[0] for j in right),default=inst.H)
        cuts.append(dict(epoch=c,left_clear=a,right_first=b,left_margin=c-a,right_margin=b-c,
                         valid=c-a>=diameter and b-c>=diameter))
    return dict(valid=all(c['valid'] for c in cuts),cuts=cuts,maximum_shortest_road=diameter)


def pack_limits(inst):
    """Partition mandatory aircraft intervals into cliques; no schedule search."""
    intervals=[];P=len(inst.pads)
    for j,job in enumerate(inst.jobs):
        left=job.hi-min(inst.fly(p,P+j) for p in range(P))
        right=job.lo+inst.sampling_duration+min(inst.fly(P+j,p) for p in range(P))+inst.ground_duration
        intervals.append((left,right))
    upper=0;clique=[];lo=hi=None
    for a,b in sorted(intervals):
        if a>=b:
            upper+=1;continue
        if clique and max(lo,a)>=min(hi,b):
            upper+=min(inst.drones,len(clique));clique=[]
        clique.append((a,b));lo=max(lo,a) if len(clique)>1 else a;hi=min(hi,b) if len(clique)>1 else b
    upper+=min(inst.drones,len(clique))
    n=len(inst.jobs)
    return max(0,n-(inst.groups-1)*upper),upper
