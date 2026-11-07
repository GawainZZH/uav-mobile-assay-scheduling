"""Exact integral path decomposition with explicit consumable stock."""
from fractions import Fraction as F
from collections import defaultdict


def decompose(arcs,flows,*,groups,pads,budget,tasks):
    if type(groups)is not int or groups<1 or type(pads)is not int or pads<1 or type(budget)is not int or budget<0:
        raise ValueError('integer source domain required')
    tasks=tuple(tuple(s)for s in tasks);K=len(tasks);all_tasks=[j for s in tasks for j in s]
    if not K or len(set(all_tasks))!=len(all_tasks)or len(arcs)!=len(flows):raise ValueError('partition and flow dimensions')
    stock=[];balance=defaultdict(int);cover={j:0 for j in all_tasks};out=defaultdict(list);cost=F(0)
    for i,(e,v)in enumerate(zip(arcs,flows)):
        v=F(str(v))
        if v.denominator!=1 or v<0:raise ValueError('integral nonnegative flow required')
        v=int(v);stock.append(v);a=tuple(e['tail']);b=tuple(e['head']);S=tuple(e['subset'])
        if len(a)!=3 or len(b)!=3 or any(type(x)is not int for x in a+b):raise ValueError('integer node coordinates required')
        k,p,left=a;kk,q,right=b
        if not 0<=k<K or kk!=k+1 or not 0<=p<pads or not 0<=q<pads or not 0<=left<=budget or right!=left-len(S)or right<0:
            raise ValueError('stage, pad or consumable stock transition differs')
        if S!=tuple(sorted(set(S)))or not set(S)<=set(tasks[k]):raise ValueError('local task partition differs')
        if F(str(e['cost']))<0 or v>groups:raise ValueError('nonnegative cost and source capacity required')
        balance[a]+=v;balance[b]-=v;out[a].append(i);cost+=v*F(str(e['cost']))
        for j in S:cover[j]+=v
    if any(v!=1 for v in cover.values()):raise ValueError('global exact task cover differs')
    source=(0,0,budget)
    for (k,p,b),v in balance.items():
        expected=groups if(k,p,b)==source else None if k==K and p==0 else 0
        if expected is not None and v!=expected or expected is None and v>0:raise ValueError('source, conservation or terminal depot differs')
    if balance.get(source,0)!=groups:raise ValueError('source flow differs')
    paths=[];remaining=[]
    for g in range(groups):
        node=source;path=[]
        for k in range(K):
            i=next((i for i in out[node]if stock[i]),None)
            if i is None:raise ValueError('positive flow has no complete path')
            stock[i]-=1;path.append(i);node=tuple(arcs[i]['head'])
        if node[0]!=K or node[1]!=0:raise ValueError('path does not return to depot')
        paths.append(tuple(path));remaining.append(node[2])
    if any(stock):raise ValueError('unused flow component')
    return dict(paths=tuple(paths),remaining=tuple(remaining),objective=cost,coverage=cover)
