"""Lazy original one-group kernels; no imported answer, bound or seed."""
import time
from dataclasses import replace
from fractions import Fraction as F
from service_separator25_joint.solve import signature
from group_subset_oracle.model import validate_subset


class Oracle:
    def __init__(self,domain,deadline,*,key_seconds=5,kernel_factory=None):
        if domain.inst is None:raise ValueError('original physical instance required')
        if key_seconds<=0:raise ValueError('positive local query budget required')
        self.d=domain;self.deadline=deadline;self.key_seconds=key_seconds;self.factory=kernel_factory
        self.kernels={};self.cache={};self.events=[];self.build_seconds=0.;self.query_seconds=0.

    def validate(self,i,plan):
        f=self.d.families[i];local,_=self.d.stages[f['stage']]
        v=validate_subset(local,plan,f['local_subset'],f['start'],f['end'],crew_capacity=self.d.inst.drones)
        if not v['valid']:raise ValueError(('local physical plan invalid',v['errors']))
        return F(str(v['objective']))

    def query(self,i):
        began=time.perf_counter();f=self.d.families[i];local,_=self.d.stages[f['stage']];template=signature(local)
        key=(template,tuple(f['local_subset']),f['start'],f['end'])
        if key in self.cache:
            r=self.cache[key];self.events.append(dict(family=i,cached=True,query_key=key,seconds=0.));return r
        if time.perf_counter()>=self.deadline:raise TimeoutError('local oracle global deadline')
        if template not in self.kernels:
            q=time.perf_counter();factory=self.factory
            if factory is None:
                from group_subset_oracle.kernel import SubsetKernel
                factory=SubsetKernel
            # Building is charged to the same global budget, not hidden in a per-key timer.
            self.kernels[template]=factory(replace(local,groups=1),self.deadline,crew_capacity=self.d.inst.drones)
            self.build_seconds+=time.perf_counter()-q
        r=self.kernels[template].query(f['local_subset'],f['start'],f['end'],min(self.deadline,time.perf_counter()+self.key_seconds))
        if r.get('objective')is not None and self.validate(i,r['plan'])!=F(str(r['objective'])):
            raise ValueError('G1 returned a wrong physical objective')
        self.cache[key]=r;elapsed=time.perf_counter()-began;self.query_seconds+=elapsed
        self.events.append(dict(family=i,cached=False,query_key=key,seconds=elapsed,record=r))
        return r

    def close(self):
        for k in self.kernels.values():k.close()

