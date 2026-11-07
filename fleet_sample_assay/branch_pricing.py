"""Exact label pricing with complete context-specific adjacent-action prices."""
import time
from .physical import Route,contexts


def actions(route):
    chain=(('L',route.p,route.launch),)+tuple(('J',j,s) for j,s in route.visits)+(('R',route.q,route.recovery),)
    return tuple((route.context,a,b) for a,b in zip(chain,chain[1:]))


def price(inst,common,visit,edge,deadline=None):
    best=float('inf');winner=None;labels=0;conditions=0;m=len(inst.pads)
    for ctx in contexts(inst):
        if deadline is not None and time.perf_counter()>deadline:raise TimeoutError('action pricing budget')
        conditions+=1;g,p,l,q,e=ctx;layer={(0,p,l):(common(ctx),())}
        for count in range(inst.bottles+1):
            nxt={}
            for (mask,last,now),(value,path) in layer.items():
                labels+=1;previous=('J',last-m,now-1) if mask else ('L',p,l)
                if mask and now+inst.fly(last,q)<=e:
                    final=value+edge((ctx,previous,('R',q,e)))
                    if final<best:best=final;winner=Route(*ctx,path)
                if count==inst.bottles:continue
                volume=sum(job.volume for j,job in enumerate(inst.jobs) if mask>>j&1)
                for j,job in enumerate(inst.jobs):
                    if mask>>j&1 or volume+job.volume>inst.payload:continue
                    low=max(job.lo,now+inst.fly(last,m+j),e-job.shelf)
                    high=min(job.hi,e-1-inst.fly(m+j,q))
                    if e>min(inst.H-1,job.deadline-1):continue
                    for s in range(low,high+1):
                        state=mask|1<<j,m+j,s+1
                        val=value+visit(ctx,j,s)+edge((ctx,previous,('J',j,s)))
                        if state not in nxt or val<nxt[state][0]:nxt[state]=val,path+((j,s),)
            layer=nxt
    return dict(value=best,route=winner,labels=labels,conditions=conditions)
