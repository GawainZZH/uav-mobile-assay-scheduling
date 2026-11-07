"""Exact finite-box LP certificates compiled into subset/endpoint functions."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import ceil
from hashlib import sha256
from review88_physical.evidence import encode,instance,plain
from service_separator25_joint.solve import signature
from group_subset_oracle.trial import checked_price_matrix


def integer(v):
    x=F(str(v))
    if x.denominator!=1:raise ValueError('exact integer certificate coordinate required')
    return x.numerator


@dataclass(frozen=True)
class AffinePrice:
    constant: int
    cover: tuple
    start: tuple
    end: tuple
    scale: int
    cost_scale: int
    instance_signature: str
    crew_capacity: int
    certificate_id: str
    model_sha256: str

    def numerator(self,S,start,end):
        S=tuple(S)
        if S!=tuple(sorted(set(S)))or any(type(j)is not int or not 0<=j<len(self.cover)for j in S):raise ValueError('canonical local subset required')
        if any(type(p)is not int or not 0<=p<len(self.start)for p in(start,end)):raise ValueError('valid endpoint required')
        return self.constant+sum(self.cover[j]for j in S)+self.start[start]+self.end[end]

    def lower(self,S,start,end,*,rounded=True):
        value=F(self.numerator(S,start,end),self.scale)
        return F(ceil(value),self.cost_scale)if rounded else value/self.cost_scale


class PricePool:
    """The caller binds this model to the audited original G1 producer/source."""
    def __init__(self,model):
        model=plain(model)
        self.model=model;M=model['matrix'];self.keys=tuple(tuple(k)for k in model['row_keys'])
        self.inst=inst=instance(model['instance']);self.D=integer(model['cost_scale']);self.crew=integer(model['crew_capacity'])
        if inst.groups!=1 or self.D<1 or not 1<=self.crew<=inst.drones:raise ValueError('original G1 scale/crew domain required')
        nr,nc=M['shape']
        if not(len(self.keys)==len(M['rhs'])==len(M['equalities'])==nr and len(M['upper'])==len(M['cost'])==nc):raise ValueError('complete model arrays required')
        if len(set(self.keys))!=nr:raise ValueError('unique row roles required')
        if not len(M['row'])==len(M['col'])==len(M['data']):raise ValueError('complete sparse matrix required')
        self.cost=tuple(integer(x)for x in M['cost']);self.upper=tuple(integer(x)for x in M['upper']);self.rhs=tuple(integer(x)for x in M['rhs'])
        if min(self.upper,default=0)<0:raise ValueError('nonnegative finite column caps required')
        self.entries=[]
        for i,j,v in zip(M['row'],M['col'],M['data']):
            if type(i)is not int or type(j)is not int or not 0<=i<nr or not 0<=j<nc:raise ValueError('sparse coordinate outside model')
            self.entries.append((i,j,integer(v)))
        self.fingerprint=sha256(encode(model)).hexdigest();self.instance_signature=signature(inst);self.prices={}

    def compile(self,record):
        # This verifies every changed RHS role and the original query's scale.
        q=checked_price_matrix(self.model,record,tuple(record['subset']),record['start'],record['end'])
        lp=record.get('LP');p=None if not lp else lp.get('price')
        if p is None:raise ValueError('rational LP price required; no substitution of native integer bounds')
        w=tuple(integer(v)for v in p['weights']);s=integer(p['scale']);K=integer(p['constant'])
        if len(w)!=len(self.keys)or s<1:raise ValueError('price dimensions/scale')
        if any(v>0 for v,eq in zip(w,q['equalities'])if not eq):raise ValueError('wrong inequality price sign')
        rc=[v*s for v in self.cost]
        for i,j,a in self.entries:rc[j]-=w[i]*a
        constant=sum(u*min(0,v)for u,v in zip(self.upper,rc))
        if K!=constant:raise ValueError('box-repair constant differs from every-column reconstruction')
        original=F(K+sum(v*integer(b)for v,b in zip(w,q['rhs'])),s)
        if original!=F(str(lp['exact_dual'])):raise ValueError('original RHS price evaluation differs')
        if F(ceil(original),self.D)!=F(str(lp['bound'])):raise ValueError('cost-grid price bound differs')
        N=len(self.inst.jobs);P=len(self.inst.pads);cover={};starts={};ends={};pack=None;a0=K
        for key,b,v in zip(self.keys,self.rhs,w):
            if key[0]=='cover':cover[key[1]]=v
            elif key[0]=='pack':pack=v
            elif key[0]=='flow'and key[3]in(0,self.inst.H):
                (starts if key[3]==0 else ends)[key[2]]=v if key[3]==0 else-v
            else:a0+=b*v
        if set(cover)!=set(range(N))or set(starts)!=set(range(P))or set(ends)!=set(range(P))or pack is None:
            raise ValueError('complete cover/pack/endpoint row roles required')
        token=sha256(encode([self.fingerprint,w,K,s])).hexdigest()
        result=AffinePrice(a0,tuple(cover[j]+pack for j in range(N)),tuple(starts[p]for p in range(P)),
            tuple(ends[p]for p in range(P)),s,self.D,self.instance_signature,self.crew,token,self.fingerprint)
        if result.lower(record['subset'],record['start'],record['end'],rounded=False)!=original/self.D:
            raise ValueError('compiled affine function differs on original query')
        self.prices[token]=result
        return result
