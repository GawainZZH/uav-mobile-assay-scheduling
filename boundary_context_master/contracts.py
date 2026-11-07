"""Exact algebra for boundary choices and symmetric context cost rows.

This module neither proves local schedule feasibility nor builds an optimization
model. Group identities remain in resource_path, even though cost patterns pool
homogeneous groups. Stored lower bounds and feasible upper bounds are distinct.
"""
from collections import Counter
from fractions import Fraction


def _integer(v):
    return type(v) is int


def _signature(pattern):
    pattern=tuple(pattern)
    if not pattern or any(not isinstance(t,tuple) or len(t)!=3 or
            any(not _integer(v) or v<0 for v in t) or not _integer(c) or c<=0
            for t,c in pattern):
        raise ValueError('nonempty integer type histogram required')
    if pattern!=tuple(sorted(pattern)) or len({t for t,c in pattern})!=len(pattern):
        raise ValueError('canonical unique histogram types required')
    return pattern


def histogram(actions):
    actions=tuple(tuple(a) for a in actions)
    if not actions or any(len(a)!=3 or any(not _integer(v) or v<0 for v in a) for a in actions):
        raise ValueError('nonnegative integer boundary actions required')
    return tuple(sorted(Counter(actions).items()))


def order_bits(pattern):
    """Sparse ones of the unique monotone threshold encoding."""
    return {(t,r):1 for t,c in _signature(pattern) for r in range(1,c+1)}


def deficit(pattern, values):
    pattern=_signature(pattern)
    return Fraction(sum(c for t,c in pattern))-sum(
        (Fraction(values.get((t,r),0)) for t,c in pattern for r in range(1,c+1)),
        Fraction(0))


def _rational(value):
    if not isinstance(value,(int,Fraction)) or isinstance(value,bool):
        raise ValueError('explicit integer or Fraction bound required')
    return Fraction(value)


def cost_row(pattern, bound, floor):
    """theta >= constant + coefficient * sum(support threshold bits)."""
    pattern=_signature(pattern);bound=_rational(bound);floor=_rational(floor)
    if floor<0 or bound<floor:
        raise ValueError('bound must dominate the nonnegative common floor')
    groups=sum(c for t,c in pattern);gap=bound-floor
    return dict(support=tuple(order_bits(pattern)),coefficient=gap,
                constant=bound-gap*groups,groups=groups,bound=bound,floor=floor)


def row_rhs(row, values):
    return row['constant']+row['coefficient']*sum(
        (Fraction(values.get(bit,0))for bit in row['support']),Fraction(0))


def resource_path(actions, tasks, limits, batteries, pads, depot=0):
    """Validate labeled choices, not local schedules or fractional LP points."""
    actions=tuple(tuple(tuple(a)for a in stage)for stage in actions)
    tasks=tuple(tasks);limits=tuple(limits);batteries=tuple(batteries)
    groups=len(batteries)
    if not actions or not groups or len(actions)!=len(tasks) or len(actions)!=len(limits):
        raise ValueError('complete nonempty stage/group domains required')
    if not _integer(pads) or pads<=0 or not _integer(depot) or not 0<=depot<pads or any(
            not _integer(b) or b<0 for b in batteries):
        raise ValueError('invalid pad/stock domain')
    positions=(depot,)*groups;remaining=batteries;contexts=[]
    for stage,n,(lo,hi)in zip(actions,tasks,limits):
        if not _integer(n) or n<0 or not _integer(lo) or not _integer(hi) or not 0<=lo<=hi:
            raise ValueError('invalid task or consumption interval')
        pattern=histogram(stage)
        if len(stage)!=groups or any(p>=pads or q>=pads or not lo<=v<=hi for p,q,v in stage):
            raise ValueError('action outside group/pad/consumption domain')
        if sum(v for p,q,v in stage)!=n:
            raise ValueError('stage task coverage differs')
        if tuple(p for p,q,v in stage)!=positions:
            raise ValueError('labeled group position continuity fails')
        remaining=tuple(b-v for b,(p,q,v)in zip(remaining,stage))
        if min(remaining)<0:
            raise ValueError('individual group stock exhausted')
        positions=tuple(q for p,q,v in stage);contexts.append(pattern)
    if positions!=(depot,)*groups:
        raise ValueError('terminal depot return fails')
    return dict(contexts=tuple(contexts),remaining=remaining,
                consumed=tuple(b-r for b,r in zip(batteries,remaining)))


def dimensions(groups,pads,limits):
    """Initial coordinate inventory, without constructing a matrix or graph."""
    limits=tuple(limits)
    if not _integer(groups) or groups<=0 or not _integer(pads) or pads<=0 or not limits:
        raise ValueError('positive groups/pads and nonempty stages required')
    types=[]
    for k,(lo,hi)in enumerate(limits):
        if not _integer(lo) or not _integer(hi) or not 0<=lo<=hi:
            raise ValueError('integer consumption interval required')
        types.append((1 if k==0 else pads)*(1 if k==len(limits)-1 else pads)*(hi-lo+1))
    binary=groups*sum(types)
    return dict(types_by_stage=tuple(types),x_binary=binary,threshold_binary=binary,
                cost_continuous=len(limits),total=2*binary+len(limits))
