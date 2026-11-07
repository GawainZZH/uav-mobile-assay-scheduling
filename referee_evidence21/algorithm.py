from itertools import product
from fractions import Fraction
import math
import time
import numpy as np
from service_separator25_joint.solve import signature,context,master,concatenate
from service_separator25_joint.interface import certify,canonical,pack_limits
from service_separator25_joint.cases import split
from service_separator25_joint.kernel import Kernel
from review88_physical.milp import dump_matrix
from .envelope import envelope

def merge_context(old,new):
    if old is None:return new.copy()
    r=new.copy()
    if r['status']=='infeasible':
        if old.get('objective')is not None:raise ValueError('infeasibility contradicts saved context fleet')
        return r
    r['lower_bound']=max(old.get('lower_bound',0),new.get('lower_bound',0))
    if old.get('objective')is not None and(new.get('objective')is None or old['objective']<new['objective']):
        for key in ['objective','fleet','vector','native_vector']:r[key]=old[key]
    if r.get('objective')is not None and r['lower_bound']>=r['objective']:
        if r['lower_bound']>r['objective']+1e-6:raise ValueError('context bound exceeds incumbent')
        r.update(status='optimal',optimal=True)
    return r

class BoundedKernel(Kernel):
    def __init__(self,inst,deadline,policy='J'):
        super().__init__(inst,deadline,policy,False)
        self.envelope=envelope(inst,deadline,self.m);self.native_calls+=1;self.attempts={}
    def lower(self,start,end,consumed):
        if self.envelope['status']=='infeasible_lp':return math.inf
        return max(super().lower(start,end,consumed),self.envelope['bound'])
    def solve(self,start,end,consumed,deadline):
        key=(tuple(start),tuple(end),tuple(consumed));old=self.solved.get(key)
        r=super().solve(start,end,consumed,deadline)
        out=merge_context(old,r);self.solved[key]=out
        self.attempts[key]=self.attempts.get(key,0)+1
        return out

def solve(inst,boundaries,seconds=90,*,seed=True,use_envelope=True,progress=None):
    began=time.perf_counter();deadline=began+seconds-2;objects={};events=[];trace=[];best=math.inf;bestpath=();bestfleet=None;lb=0
    cert=certify(inst,boundaries);cuts=(0,)+tuple(c['epoch']for c in cert['cuts']if c['valid'])+(inst.H,)
    stages=[];report=dict(status='unresolved',objective=None,lower_bound=0,optimal=False,fleet=None)
    def emit(kind):
        r=dict(seconds=time.perf_counter()-began,kind=kind,lower_bound=lb,objective=None if not math.isfinite(best)else best)
        trace.append(r)
        if progress:progress(r)
    def refresh(path=None):
        nonlocal best,bestpath,bestfleet
        if path is None:
            value,path,_=master(inst,stages,objects,upper=True,deadline=deadline)
        else:value=sum(objects[e['kernel']].solved[e['context']]['objective']for e in path)
        if value<best:
            fleet=concatenate(inst,stages,objects,path)
            if fleet['check']['objective']!=value:raise ValueError('complete fleet value')
            best=value;bestpath=path;bestfleet=fleet;emit('validated_fleet')
    def evaluate(edge,quantum):
        k=objects[edge['kernel']];r=k.solve(*edge['context'],min(deadline,time.perf_counter()+quantum))
        events.append(dict(kernel=edge['kernel'],**r));return r
    try:
        for local,ids in split(inst,cuts):
            key=signature(local)
            if key not in objects:
                if use_envelope:objects[key]=BoundedKernel(local,deadline)
                else:
                    k=BoundedKernel.__new__(BoundedKernel);Kernel.__init__(k,local,deadline,'J',False)
                    k.envelope=dict(status='disabled',bound=0,seconds=0);k.attempts={};objects[key]=k
            stages.append((local,ids,key))
        lb=sum(float(objects[key].envelope['bound'])for _,_,key in stages);emit('envelope')
        seed_deadline=min(deadline,began+.25*seconds);attempt_limit=max(12,4*len(stages));attempts=0
        def dfs(stage,state,path):
            nonlocal attempts
            if stage==len(stages):refresh(tuple(path));return True
            if time.perf_counter()>=seed_deadline or attempts>=attempt_limit:return False
            local,ids,key=stages[stage];lo,hi=pack_limits(local);n=len(ids)
            reserve=sum(pack_limits(future)[0]for future,_,_ in stages[stage+1:])
            vectors=[v for v in product(range(lo,hi+1),repeat=inst.groups)if sum(v)==n and all(c<=b-reserve for c,(p,b)in zip(v,state))]
            vectors.sort(key=lambda v:(sum((b-c)**2 for c,(p,b)in zip(v,state)),v))
            end=tuple(g%len(inst.pads)for g in range(inst.groups))if stage+1<len(stages)else(0,)*inst.groups
            for consumed in vectors:
                if time.perf_counter()>=seed_deadline or attempts>=attempt_limit:return False
                start=tuple(p for p,b in state);ctx,perm=context(start,end,consumed,True)
                nxt,reorder=canonical(tuple((p,b-c)for p,c,(_,b)in zip(end,consumed,state)))
                edge=dict(stage=stage,kernel=key,context=ctx,perm=perm,reorder=reorder,state=state,next_state=nxt,consumed=consumed)
                attempts+=1;r=evaluate(edge,min(2.,max(.02,seed_deadline-time.perf_counter())))
                if r.get('objective')is not None and dfs(stage+1,nxt,path+[edge]):return True
            return False
        if seed:dfs(0,tuple((0,inst.batteries)for _ in range(inst.groups)),[])
        while time.perf_counter()<deadline:
            low,path,size=master(inst,stages,objects,deadline=deadline)
            if math.isfinite(low):lb=max(lb,low)
            elif not path:report.update(status='infeasible',infeasible_proved=True);break
            report['interface']=size;refresh();emit('complete_interface')
            if math.isfinite(best)and lb>=best:report['status']='optimal';break
            if not path:break
            unknown=next((e for e in path if e['context']not in objects[e['kernel']].solved or not objects[e['kernel']].solved[e['context']]['optimal']),None)
            if unknown is None:refresh(path);report['status']='optimal';lb=best;break
            k=objects[unknown['kernel']];attempt=k.attempts.get(unknown['context'],0)
            evaluate(unknown,min(2.*2**min(attempt,7),max(.02,deadline-time.perf_counter())))
    except TimeoutError as exc:report['reason']=str(exc)
    except Exception as exc:report.update(status='correctness_error',reason=repr(exc))
    if bestfleet is not None and report['status']!='correctness_error':
        report.update(fleet=bestfleet,objective=best,status='optimal'if lb>=best else'feasible',optimal=lb>=best)
    report.update(lower_bound=lb,effective_boundaries=cuts,separators=cert,events=events,trace=trace,
        interface_path=bestpath,envelopes={key:k.envelope for key,k in objects.items()},kernel_count=len(objects),
        evaluated_contexts=sum(len(k.solved)for k in objects.values()),native_calls=sum(k.native_calls for k in objects.values()),
        kernel_models={key:dict(instance=k.inst,row_keys=k.m.M.keys,matrix=dump_matrix(k.m,k.cost,None))for key,k in objects.items()},
        configuration=dict(seed=seed,envelope=use_envelope,LP_prices=False,initial_quantum=2,geometric_quantum=True))
    report['seconds']=time.perf_counter()-began
    if report['seconds']>seconds:report.update(status='workflow_time_limit',optimal=False)
    return report
