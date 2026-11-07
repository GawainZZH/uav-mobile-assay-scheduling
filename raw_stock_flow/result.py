"""Exact currency conversion and instance-independent objective bounds."""
from dataclasses import replace
from fractions import Fraction as F
from math import lcm,gcd,ceil,isfinite


def currency(inst):
    values=[F(inst.road_cost),F(inst.air_cost),F(inst.flight_fixed)]+[F(j.weight)for j in inst.jobs]
    D=lcm(*(v.denominator for v in values));nums=[int(v*D)for v in values];g=0
    for v in nums:g=gcd(g,abs(v))
    work=replace(inst,road_cost=nums[0],air_cost=nums[1],flight_fixed=nums[2],
        jobs=tuple(replace(j,weight=w)for j,w in zip(inst.jobs,nums[3:])))
    return work,D,g or 1


def classify(bound,upper,floor,quantum):
    floor=F(str(floor));q=F(str(quantum))
    if q<=0 or floor<0:raise ValueError('valid objective floor/grid required')
    lower=q*ceil(floor/q)
    if bound is not None:
        if not isfinite(float(bound)):raise ValueError('finite native bound required')
        lower=max(lower,q*ceil(float(bound/q)-1e-6))
    if upper is not None:
        upper=F(str(upper))
        if upper<0 or(upper/q).denominator!=1 or lower>upper:raise ValueError('bound contradicts true objective upper/grid')
    return dict(lower=lower,upper=upper,closed=upper is not None and lower==upper,raw_native_bound=bound,
        inherited_floor=floor,quantum=q,proof_scope='analytic lower plus native numerical integer bound and independently validated physical upper')
