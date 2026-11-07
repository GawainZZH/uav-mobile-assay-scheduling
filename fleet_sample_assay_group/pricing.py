"""Complete time-layer joint pricing with exact rational cost dominance."""
import math,time
from dataclasses import replace
from itertools import combinations
from fractions import Fraction as F
from fleet_sample_assay.physical import enumerate_routes,vehicle_arcs
from .exact import Scale,rational
from .column import decode


class Label:
    __slots__=('cost','parent','action')
    def __init__(self,cost,parent=None,action=None):self.cost=cost;self.parent=parent;self.action=action


class Limit(Exception):pass


class Engine:
    def __init__(self,inst,deadline):
        began=time.perf_counter();self.inst=replace(inst,groups=1)
        self.routes=enumerate_routes(self.inst,deadline);self.by_start={};self.masks=[];self.events=[]
        for i,r in enumerate(self.routes):
            if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('group physical pool build')
            self.by_start.setdefault((r.p,r.launch),[]).append(i)
            self.masks.append(sum(1<<j for j,s in r.visits))
            self.events.append((r.recovery,r.q,tuple(sorted((self.inst.jobs[j].weight,r.window(self.inst,j,s)[1]) for j,s in r.visits))))
        self.expired=[sum(1<<j for j,job in enumerate(self.inst.jobs) if job.hi<t) for t in range(inst.H+1)]
        self.arcs=vehicle_arcs(self.inst);self.road_by_start={}
        for a in self.arcs:
            if a[1]!=a[3]:self.road_by_start.setdefault((a[1],a[2]),[]).append(a)
        self.build_seconds=time.perf_counter()-began

    def price(self,prices,required=(),forbidden=(),sri=None,phase=1,seconds=30,label_limit=500000,compress=True):
        began=time.perf_counter();deadline=began+seconds;inst=self.inst;n=len(inst.jobs);stats=dict(processed_labels=0,transitions=0,
            peak_frontier=0,max_history_bits=0,compress=compress,physical_routes=len(self.routes))
        out=dict(status='time_limit',value=None,lower_bound=None,column=None)
        def clock():
            if time.perf_counter()>=deadline:raise TimeoutError('group pricing budget')
        try:
            clock()
            if len(prices)!=n or phase not in (0,1):raise ValueError('group price dimensions/phase')
            if any(type(j) is not int or not 0<=j<n for j in tuple(required)+tuple(forbidden)):raise ValueError('branch task domain')
            sri={tuple(sorted(S)):v for S,v in (sri or {}).items()}
            if any(len(S)!=3 or len(set(S))!=3 or min(S)<0 or max(S)>=n for S in sri):raise ValueError('SRI pricing domain')
            req=sum(1<<j for j in set(required));bad=sum(1<<j for j in set(forbidden))
            if req&bad:return {**out,**stats,'status':'infeasible','seconds':time.perf_counter()-began}
            if not all(math.isfinite(float(v)) for v in list(prices)+list(sri.values())):raise ValueError('finite pricing duals required')
            integral_physical=all(type(v) is int for v in [inst.flight_fixed,inst.air_cost,inst.road_cost]+[j.weight for j in inst.jobs])
            if not req and max(prices,default=0)<=0 and max(sri.values(),default=0)<=0 and (phase==0 or integral_physical):
                actions=[('wait',t,0,(),()) for t in range(inst.H)]
                column=decode(inst,self.routes,actions,required,forbidden);clock()
                return {**out,**stats,'status':'optimal','value':0.,'lower_bound':0.,'exact_value':('0','1'),
                    'column':column,'certificate_kind':'nonnegative_objective_empty_column','seconds':time.perf_counter()-began}
            values=list(prices)+list(sri.values())+[r.cost(inst) for r in self.routes]+[a[5] for a in self.arcs]
            weights=set(j.weight for j in inst.jobs)
            values.extend(w*(k+1) for w in weights for k in range(inst.H))
            scale=Scale(values);pi=[scale.integer(v) for v in prices];sigma={S:scale.integer(v) for S,v in sri.items()}
            assay={(w,k):scale.integer(w*(k+1))*phase for w in weights for k in range(inst.H)}
            route_cost=[]
            for i,r in enumerate(self.routes):
                if i%128==0:clock()
                jobs={j for j,s in r.visits}
                route_cost.append(scale.integer(r.cost(inst))*phase-sum(pi[j] for j in jobs)-sum(v for S,v in sigma.items() if len(jobs&set(S))>=2))
            roads={a:scale.integer(a[5])*phase for a in self.arcs if a[1]!=a[3]}
            layers=[{} for t in range(inst.H+1)];start=(0,0,(),(),())
            layers[0][start]=Label(0);terminal=None
            def offer(t,u,state,cost,parent,action):
                p,mask,flying,releases,pending=state
                expired=self.expired[u]^self.expired[t]
                if req&expired&~mask:return
                if compress:mask&=~self.expired[u]
                releases=tuple(x for x in releases if x>u)
                if any(d<u for w,d in pending):return
                key=p,mask,flying,releases,pending
                old=layers[u].get(key)
                if old is None or cost<old.cost:layers[u][key]=Label(cost,parent,action)
                stats['transitions']+=1
            for t in range(inst.H+1):
                clock();current=layers[t];stats['peak_frontier']=max(stats['peak_frontier'],sum(len(a) for a in layers[t:]))
                if t==inst.H:
                    for (p,mask,flying,releases,pending),label in current.items():
                        if p==0 and not flying and not pending and (terminal is None or label.cost<terminal.cost):terminal=label
                    clock()
                    break
                for (p,mask,flying,releases,pending),label in current.items():
                    if stats['processed_labels']>=label_limit:raise Limit('label_limit')
                    stats['processed_labels']+=1
                    if stats['processed_labels']%64==0:clock()
                    stats['max_history_bits']=max(stats['max_history_bits'],mask.bit_count())
                    arrivals=[f for f in flying if f[0]==t]
                    if any(q!=p for e,q,payload in arrivals):continue
                    active=tuple(f for f in flying if f[0]>t)
                    available=tuple(sorted(pending+tuple(token for e,q,payload in arrivals for token in payload)))
                    if any(d<t for w,d in available):continue
                    free=min(inst.drones-len(active),inst.batteries-len(releases))
                    candidates=[i for i in self.by_start.get((p,t),()) if not self.masks[i]&(mask|bad)]
                    k=min(inst.stations,len(available));choices=set(tuple(available[i] for i in ids) for ids in combinations(range(len(available)),k))
                    def launches(pos,remaining,used,selected):
                        yield selected,used
                        if remaining:
                            for z in range(pos,len(candidates)):
                                i=candidates[z]
                                if not self.masks[i]&used:
                                    yield from launches(z+1,remaining-1,used|self.masks[i],selected+(i,))
                    generated=0
                    for selected,newmask in launches(0,free,mask,()):
                        generated+=1
                        if generated%64==0:clock()
                        future=tuple(sorted(active+tuple(self.events[i] for i in selected)))
                        release=tuple(sorted(releases+tuple(self.routes[i].recovery+inst.charge for i in selected)))
                        rc=sum(route_cost[i] for i in selected)
                        # Original endpoints require a path node, not a wait arc:
                        # recover/launch first, then depart immediately if desired.
                        for a in self.road_by_start.get((p,t),()):
                            g,_,_,q,u,c=a
                            if q!=0 and u+inst.road(q,0)>inst.H:continue
                            if any(e<u for e,qf,payload in future):continue
                            offer(t,u,(q,newmask,future,release,available),label.cost+rc+roads[a],label,('move',t,p,q,u,selected))
                        for tested in choices:
                            rest=list(available)
                            for token in tested:rest.remove(token)
                            if any(d==t for w,d in rest):continue
                            value=label.cost+rc+sum(assay[w,t] for w,d in tested)
                            offer(t,t+1,(p,newmask,future,release,tuple(rest)),value,label,('wait',t,p,selected,tested))
                layers[t]={}
            if terminal is None:out['status']='infeasible';return {**out,**stats,'seconds':time.perf_counter()-began}
            actions=[];label=terminal
            while label.parent is not None:actions.append(label.action);label=label.parent
            actions.reverse();column=decode(inst,self.routes,actions,required,forbidden)
            actual=F(*map(int,column['exact_cost']))*phase-sum((rational(prices[j]) for j in column['covered']),F(0))
            for S,v in sri.items():
                count=sum(len({j for j,s in r.visits}&set(S))>=2 for r in column['routes']);actual-=rational(v)*count
            if actual!=F(terminal.cost,scale.denominator):raise ValueError('decoded exact reduced cost mismatch')
            clock();out.update(status='optimal',column=column,certificate_kind='complete_exact_layer_DP',**scale.result(terminal.cost))
        except TimeoutError:out['status']='time_limit'
        except Limit:out['status']='label_limit'
        out.update(stats,seconds=time.perf_counter()-began)
        return out
