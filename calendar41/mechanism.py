"""Complete four-task instance enumeration with staggered verification."""
from itertools import permutations,product
from pathlib import Path
import argparse
import sys
import time
from review88_physical.types import Instance,Job,Route
from review88_physical.physical import validate_fleet
from review88_physical.evidence import encode,plain


def instance():
    jobs=tuple(Job(xy,t,t,2,deadline,1,name)for xy,t,deadline,name in[
        ((0,1),1,5,'old0'),((0,1),1,5,'old1'),
        ((5,1),16,20,'next0'),((5,1),16,20,'next1')])
    return Instance('staggered-two-group',30,((0,0),(5,0)),jobs,groups=2,drones=1,
        batteries=2,stations=1,assay_duration=2,ground_duration=3,handoff_delay=0,
        endurance=3,flight_fixed=2,air_cost=1,road_cost=1,air_speed=1)


def calendars():
    return (((5,6),(13,14),(24,25)),((5,6),(8,10),(24,25)))


def group_plan(start,old,new):
    trail=[(0,t,0,t+1)for t in range(start)]+[(0,start,1,start+5)]
    trail.extend((1,t,1,t+1)for t in range(start+5,25));trail.append((1,25,0,30))
    return {'routes':[Route(0,0,0,3,old,1),Route(1,15,1,18,new,16)],
            'assays':[(old,1,3,3),(new,16,18,18)],'trail':trail}


def fleet(first,second,old=(0,1),new=(2,3)):
    return [group_plan(first,old[0],new[0]),group_plan(second,old[1],new[1])]


def validate(inst,calendar,plans,promise=10):
    base=validate_fleet(inst,plans);errors=list(base['errors']);conflicts=[];releases={}
    if len(calendar)!=inst.groups:raise ValueError('one calendar per group required')
    for g,(reserved,p)in enumerate(zip(calendar,plans)):
        previous=0
        for a,b in reserved:
            if type(a)is not int or type(b)is not int or not 0<=previous<=a<b<=inst.H:
                raise ValueError('disjoint ordered verification slots required')
            previous=b
        waiting={t for p,t,q,u in p['trail']if p==q and u==t+1}
        for a,b in reserved:
            if not all(t in waiting for t in range(a,b)):
                errors.append('verification requires parked truck');conflicts.append(g)
            if any(max(start,a)<min(start+inst.assay_duration,b)for j,s,e,start in p['assays']):
                errors.append('assay overlaps verification');conflicts.append(g)
        for j,s,e,start in p['assays']:
            finish=start+inst.assay_duration;later=[b for a,b in reserved if a>=finish]
            if not later:errors.append('missing subsequent verification');continue
            release=min(later);releases[j]=release
            if release>min(inst.H,inst.jobs[j].lo+promise):errors.append('verified result promise')
    weighted_release=sum(inst.jobs[j].weight*r for j,r in releases.items())
    return {'valid':not errors,'errors':sorted(set(errors)),'base':base,
            'calendar_conflicts':sorted(set(conflicts)),'releases':releases,
            'transport':base['road_cost']+base['flight_cost'],'weighted_release':weighted_release}


def free_gap(calendar,a,b,cut):
    points=sorted({a,b,cut}|{max(a,min(b,x))for v in calendar for x in v})
    return max((y-x for x,y in zip(points,points[1:])
                if not any(max(x,l)<min(y,r)for l,r in calendar)),default=0)


def enumerate_instance():
    inst=instance();c=calendars();feasible=[];normalized={};candidates=0
    for old,new,first,second in product(permutations((0,1)),permutations((2,3)),range(6,11),range(6,11)):
        candidates+=1;p=fleet(first,second,old,new);checked=validate(inst,c,p)
        if not checked['valid']:continue
        q=fleet(8,10,old,new);normalized_check=validate(inst,c,q)
        entry={'old_assignment':old,'new_assignment':new,'road_starts':(first,second),
               'plans':p,'check':checked,'normalized_plans':q,'normalized_check':normalized_check}
        feasible.append(entry);normalized[(old,new)]=q
    swapped=fleet(10,8,old=(0,1),new=(3,2));swap_check=validate(inst,c,swapped)
    return {'instance':inst,'calendars':c,'cut':8,'candidates':candidates,'feasible':feasible,
            'normalized':list(normalized.values()),'individual_free_gaps':[free_gap(v,6,15,8)for v in c],
            'common_free_gap':free_gap(c[0]+c[1],6,15,8),'swap':{'plans':swapped,'check':swap_check},
            'scope':'Complete instance under the stated singleton, grid, calendar and deadline model. '
                    'Completeness uses the forced-event argument in the accompanying proof.',
            'new_optimizer_calls':0,'performance_claim':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    began=time.perf_counter();calls=[]
    def guard(frame,event,arg):
        name=frame.f_code.co_name if event=='call'else getattr(arg,'__name__','')if event=='c_call'else''
        if name in('optimize','solve_mip','solve_lp','milp','linprog'):
            calls.append(name);raise AssertionError('optimizer forbidden')
    previous=sys.getprofile();sys.setprofile(guard)
    try:r=enumerate_instance()
    finally:sys.setprofile(previous)
    if r['candidates']!=100 or len(r['feasible'])!=12 or len(r['normalized'])!=4:raise ValueError('predeclared population differs')
    if any(not x['normalized_check']['valid']or x['check']['base']['objective']!=56 for x in r['feasible']):
        raise ValueError('physical normalization differs')
    elapsed=time.perf_counter()-began
    if elapsed>60:raise TimeoutError('finite verification budget')
    r.update(status='PASS',seconds=elapsed,new_optimizer_calls=len(calls))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('xb')as f:f.write(encode(r))
    print({'status':r['status'],'candidates':100,'feasible':12,'normalized':4,'objective':56,
           'individual_gaps':r['individual_free_gaps'],'common_gap':r['common_free_gap'],'seconds':elapsed,'optimizer_calls':len(calls)})


if __name__=='__main__':main()
