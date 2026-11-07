"""A true incumbent does not require proving its local plans minimal."""
import time
from fractions import Fraction as F
from types import SimpleNamespace as NS
from review88_physical.physical import validate_fleet
from block_flow_solver.root import check
from block_price_tree.policy import Inspector


def verified_source(m,vector,fleet):
    x=[F(str(v))for v in vector]
    if len(x)!=len(m.cost)or any(not 0<=v<=u or flag and v.denominator!=1 for v,u,flag in zip(x,m.upper,m.integer)):
        raise ValueError('source domain differs')
    for r in m.rows:
        value=sum((a*x[j]for j,a in r['entries'].items()),F(0))-r['rhs']
        if value!=0 if r['eq']else value>0:raise ValueError('source exact row violation')
    actual=validate_fleet(m.domain.inst,fleet['group_plans']);cost=sum((c*v for c,v in zip(m.cost,x)),F(0))
    if not actual['valid']or cost!=F(str(actual['objective'])):raise ValueError('source physical cost mismatch')
    return NS(vector=x,solution=dict(objective=cost,fleet=dict(group_plans=fleet['group_plans'],check=actual),physical=actual,
        inherited_true_incumbent=True))


class WarmInspector(Inspector):
    def __init__(self,m,query,validate,reconstruct,deadline,source):
        super().__init__(m,query,validate,reconstruct,deadline)
        self.source=source;self.best=source.solution;self.point=list(source.vector);self.source_admissions=0

    def inspect(self,x):
        if time.perf_counter()<self.deadline:
            self.capture()
            if not self.pending and len(x)==len(self.source.vector)and max((abs(float(a)-float(b))for a,b in zip(x,self.source.vector)),default=0)<=1e-7:
                try:check(self.m,x,integer=True)
                except ValueError:pass
                else:
                    if abs(sum(float(c)*v for c,v in zip(self.m.cost,x))-float(self.source.solution['objective']))<=1e-6:
                        self.checks+=1;self.validated+=1;self.source_admissions+=1;return 'feasible'
        return super().inspect(x)
