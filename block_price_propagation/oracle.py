"""Root-LP adapter: each acquired price updates all matching families.

Use with the rebuilding root driver. This adapter is not a SCIP tree handler:
tree use additionally needs native constraint insertion and point rechecking.
"""
import time
from review88_physical.evidence import instance,plain
from service_separator25_joint.solve import signature
from .prices import PricePool
from .propagate import apply_prices


class SharedOracle:
    def __init__(self,m,base,*,model_for_query=None):
        self.m=m;self.base=base;self.provider=model_for_query or self._model
        self.pools={};self.seen=set();self.changes=[];self.events=[]
        self.compilation_seconds=0.;self.propagation_seconds=0.

    def _model(self,i):
        f=self.m.domain.families[i];local,_=self.m.domain.stages[f['stage']]
        k=self.base.kernels[signature(local)]
        return dict(instance=k.t.inst,crew_capacity=k.t.crew_capacity,cost_scale=k.t.cost_scale,matrix=k.matrix,row_keys=k.t.m.M.keys)

    def query(self,i):
        result=self.base.query(i);lp=result.get('LP')
        if not lp or not lp.get('price'):return result
        began=time.perf_counter();model=self.provider(i)
        key=(signature(instance(plain(model['instance']))),model['cost_scale'],model['crew_capacity'])
        if key not in self.pools:self.pools[key]=PricePool(model)
        else:
            expected=self.pools[key].model
            if expected['matrix']!=plain(model['matrix'])or expected['row_keys']!=plain(model['row_keys']):
                raise ValueError('a reused template model changed after price compilation')
        p=self.pools[key].compile(result);self.compilation_seconds+=time.perf_counter()-began
        if p.certificate_id not in self.seen:
            began=time.perf_counter();changes=apply_prices(self.m,[p]);self.propagation_seconds+=time.perf_counter()-began
            self.seen.add(p.certificate_id);self.changes.extend(changes)
            self.events.append(dict(source_family=i,certificate_id=p.certificate_id,updates=changes))
        return result
