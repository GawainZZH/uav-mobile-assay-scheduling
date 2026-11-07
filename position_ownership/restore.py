"""Complete real carrier paths before reconstructing old auxiliary coordinates."""
from copy import copy
from fractions import Fraction as F
from types import SimpleNamespace as NS
import numpy as np
from boundary_context_master.matrix import check_vector
from boundary_context_master.common import embed_source
from .model import lift


def strong_variant(p,*,keep_positions_integer=False):
    out=copy(p);out.matrix=dict(p.matrix);flags=[0]*len(p.matrix['cost'])
    flags[:p.physical_columns]=p.base['integer'][:p.physical_columns]
    if keep_positions_integer:
        for j in p.y.values():flags[j]=1
    out.matrix['integer']=flags
    out.adapter=NS(A=p.adapter.A,M=NS(cost=p.matrix['cost'],upper=p.matrix['upper'],rhs=p.matrix['rhs'],
        eq=p.matrix['equalities'],integer=flags))
    return out


def complete_real_fleet(p,vector,*,fixed_positions=None):
    from service_separator25_joint.model import complete
    from review88_physical.milp import restore
    fixed_positions={}if fixed_positions is None else dict(fixed_positions)
    # Every declared integer position is a compulsory branch condition.
    for c,g,a in p.y:
        if p.matrix['integer'][p.y[c,g,a]]:
            selected=[q for q in range(len(p.inst.pads))if vector[p.y[c,g,q]]>.5]
            if len(selected)!=1:raise ValueError('integer position not one-hot')
            actual=selected[0]
            if(c,g)in fixed_positions and fixed_positions[c,g]!=actual:raise ValueError('conflicting fixed boundary position')
            fixed_positions[c,g]=actual
    for(c,g),a in fixed_positions.items():
        if type(c)is not int or not 1<=c<len(p.d.stages)or type(g)is not int or not 0<=g<p.inst.groups or type(a)is not int or not 0<=a<len(p.inst.pads):
            raise ValueError('fixed boundary outside physical domain')
        if abs(float(vector[p.y[c,g,a]])-1)>1e-6:raise ValueError('fixed boundary disagrees with current physical flow')
    check_vector(p.adapter,vector)
    raw=np.asarray(vector[:p.base['shape'][1]],dtype=float).copy();m=p.m
    services=tuple(m.rids.values())+tuple(m.assignments.values())
    if np.max(abs(raw[list(services)]-np.rint(raw[list(services)])),initial=0)>1e-6:raise ValueError('integer service events required')
    raw[list(services)]=np.rint(raw[list(services)])
    # Recompute ownership marginals after numerical service normalization.
    normalized=lift(p,raw);check_vector(p.adapter,normalized)
    saved=list(m.M.upper)
    try:
        for(c,g),a in fixed_positions.items():
            cut=p.d.cuts[c]
            for(gg,i),col in m.zids.items():
                pp,t,q,u=m.arcs[i]
                if gg==g and t==cut and pp!=a:m.M.upper[col]=0
        completed=complete(m,raw,(0,)*p.inst.groups,(0,)*p.inst.groups)
        fleet,physical=restore(m,completed)
    finally:m.M.upper[:]=saved
    old,_,_=embed_source(p.d,p.a,fleet);actual=lift(p,old);check_vector(p.adapter,actual)
    for(c,g),a in fixed_positions.items():
        if abs(actual[p.y[c,g,a]]-1)>1e-9:raise ValueError('completion did not preserve fixed position')
    before=sum(float(c)*float(x)for c,x in zip(p.matrix['cost'],vector))
    if float(F(physical['objective']))>before+1e-5:raise ValueError('canonical carrier completion increased real objective')
    return fleet,actual
