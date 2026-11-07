"""Route-local rank and cost coefficients; no optimizer or new history state."""
import itertools,math,time
from fractions import Fraction as F


def down_float(value):
    """Convert an exact coefficient without accidentally increasing a <= row."""
    exact=value if isinstance(value,F) else F(value)
    result=float(exact)
    if not math.isfinite(result):raise ValueError('nonfinite completed coefficient')
    return math.nextafter(result,-math.inf) if F.from_float(result)>exact else result


def rank(masks,T,stations):
    best=0;unions=[0]*(1<<len(masks))
    for subset in range(1,1<<len(masks)):
        bit=subset&-subset
        union=unions[subset^bit]|masks[bit.bit_length()-1];unions[subset]=union
        best=max(best,subset.bit_count()-stations*(union&~T).bit_count())
    return best


def batch_bound(windows,weights,beta,stations):
    """K cheapest distinct slots per sample suffice for an exact local match.

    Integral prices give integer arithmetic. Other prices use downward rounded
    edge costs and sums, so this coefficient is conservative as well as tight
    to floating precision. Infeasible local columns return infinity.
    """
    K=len(windows)
    if len(weights)!=K or stations<1:raise ValueError('batch dimensions')
    if not K:return 0
    integral=all(math.isfinite(float(v)) and float(v).is_integer() for v in tuple(beta)+tuple(weights))
    choices=[]
    for W,w in zip(windows,weights):
        values=[]
        for k in W:
            value=int(w)*(k+1)-int(beta[k]) if integral else down_float(F(w)*(k+1)-F(float(beta[k])))
            values.append((value,k))
        choices.append(sorted(values)[:K])
    best=math.inf
    for assignment in itertools.product(*choices):
        slots=[v[1] for v in assignment]
        if any(slots.count(k)>stations for k in set(slots)):continue
        value=0
        for c,k in assignment:value=value+c if integral else math.nextafter(value+c,-math.inf)
        best=min(best,value)
    return best


class RouteKernel:
    """Complete-route coefficients; caches contain signatures, not LP variables."""
    def __init__(self,inst,routes,deadline=math.inf):
        self.inst=inst;self.routes=tuple(routes);self.windows=[];self.masks=[];self.signatures=[]
        self.local_bad=[];self.rank_cache={}
        for i,r in enumerate(routes):
            if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('route kernel build budget')
            W=tuple(tuple(range(e,d+1)) for e,d in (r.window(inst,j,s) for j,s in r.visits))
            masks=tuple(sum(1<<k for k in ww) for ww in W)
            signature=(r.g,tuple(sorted((ww[0],ww[-1],inst.jobs[j].weight) for ww,(j,s) in zip(W,r.visits))))
            self.windows.append(W);self.masks.append(masks);self.signatures.append(signature)
            if rank(masks,0,inst.stations):self.local_bad.append(i)
        self.bad=set(self.local_bad)

    def ranks(self,g,T,deadline=math.inf):
        key=(g,T)
        if key not in self.rank_cache:
            found=[]
            for i,r in enumerate(self.routes):
                if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('global rank completion budget')
                if r.g==g:
                    value=rank(self.masks[i],T,self.inst.stations)
                    if value:found.append((i,value))
            self.rank_cache[key]=tuple(found)
        return self.rank_cache[key]

    def costs(self,beta,deadline=math.inf):
        cache={};out=[]
        for i,(r,sig) in enumerate(zip(self.routes,self.signatures)):
            if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('global cost completion budget')
            if i in self.bad:out.append(0);continue
            if sig not in cache:
                W=self.windows[i];weights=tuple(self.inst.jobs[j].weight for j,s in r.visits)
                cache[sig]=batch_bound(W,weights,beta[r.g],self.inst.stations)
                if not math.isfinite(cache[sig]):raise ValueError('undetected locally infeasible route')
            out.append(cache[sig])
        return out

    def ordinary_costs(self,beta,deadline=math.inf):
        """Retain old price inequalities with conservative finite coefficients."""
        cache={};out=[]
        for i,r in enumerate(self.routes):
            if i%128==0 and time.perf_counter()>=deadline:raise TimeoutError('ordinary cost completion budget')
            value=F(0)
            for W,(j,s) in zip(self.windows[i],r.visits):
                key=r.g,j,s,r.recovery
                if key not in cache:
                    w=F(self.inst.jobs[j].weight)
                    cache[key]=min(w*(k+1)-F(float(beta[r.g][k])) for k in W)
                value+=cache[key]
            out.append(down_float(value))
        return out
