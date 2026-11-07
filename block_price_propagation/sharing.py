"""Equivalent price information in the same-problem strong formulation."""
from fractions import Fraction as F
from math import gcd,lcm
from block_flow_solver.matrix import row,add
from block_flow_solver.sharing import conditional_rows


def affine_row(price,cost,ownership,start,end):
    if set(ownership)!=set(range(len(price.cover))):raise ValueError('complete original task marginals required')
    e={j:-price.scale*price.cost_scale*F(str(c))for j,c in cost.items()if c};constant=price.constant
    used=list(ownership.values())
    for j,col in ownership.items():add(e,col,price.cover[j])
    for mapping,coefficients in((start,price.start),(end,price.end)):
        if type(mapping)is int:
            if not 0<=mapping<len(coefficients):raise ValueError('fixed endpoint outside domain')
            constant+=coefficients[mapping]
        else:
            if set(mapping)!=set(range(len(coefficients))):raise ValueError('complete endpoint marginals required')
            used.extend(mapping.values())
            for p,col in mapping.items():add(e,col,coefficients[p])
    if len(set(used))!=len(used)or set(used)&set(cost):raise ValueError('separate cost and feature coordinates required')
    r=row(e,-constant,False,'shared_G1_affine_price')
    # Exact row scaling reduces coefficients; it does not round a fractional RHS.
    den=lcm(*(v.denominator for v in list(r['entries'].values())+[r['rhs']]))
    nums=[int(v*den)for v in list(r['entries'].values())+[r['rhs']]];div=0
    for x in nums:div=gcd(div,abs(x))
    factor=F(den,div or 1)
    r['entries']={j:v*factor for j,v in r['entries'].items()};r['rhs']*=factor
    r['certificate_id']=price.certificate_id;r['validity']='unrounded_parametric_LP_price'
    return r


def rounded_rows(p,changes):
    rows=[];K=len(p.d.stages)
    for change in changes:
        k=change['stage']
        if type(k)is not int or not 0<=k<K:raise ValueError('stage outside original domain')
        local,ids=p.d.stages[k][:2];S=tuple(change['local_subset']);a=change['start'];b=change['end']
        if S!=tuple(sorted(set(S)))or any(type(j)is not int or not 0<=j<len(ids)for j in S):raise ValueError('original local task subset required')
        if tuple(ids[j]for j in S)!=tuple(change['subset']):raise ValueError('global/local labels differ')
        if any(type(x)is not int or not 0<=x<len(local.pads)for x in(a,b))or k==0 and a!=0 or k==K-1 and b!=0:
            raise ValueError('original endpoint domain required')
        for g in range(p.inst.groups):
            ends=(()if k==0 else(p.y[k,g,a],))+(()if k==K-1 else(p.y[k+1,g,b],))
            new=conditional_rows(p.costs[k,g],ends,{j:p.s[k,g,j]for j in range(len(ids))},S,lower=change['lower'])
            for r in new:r.update(certificate_id=change['certificate_id'],family=change['family'],stage=k,group=g)
            rows.extend(new)
    return rows
