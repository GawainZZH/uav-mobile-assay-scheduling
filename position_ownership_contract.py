"""Pure cost predicates for sufficient boundary-position/ownership integers.

No optimizer, matrix construction or experiment. Old auxiliary type variables
may remain continuous; they must not be the trigger of these conditions.
"""
from fractions import Fraction as F


def delta(endpoint_values,ownership_values,subset):
    S=set(subset)
    return len(endpoint_values)-sum(endpoint_values)+sum(1-v if j in S else v for j,v in enumerate(ownership_values))


def conditional_row(cost,endpoints,ownership,subset,bound):
    endpoints=tuple(endpoints);S=tuple(subset)
    if len(endpoints)>2 or len(set(endpoints))!=len(endpoints):raise ValueError('zero to two distinct free endpoint bits required')
    if S!=tuple(sorted(set(S)))or any(j not in ownership for j in S):raise ValueError('canonical task subset required')
    bits=tuple(ownership.values())
    if len(set(bits))!=len(bits)or set(endpoints)&set(bits)or set(cost)&(set(endpoints)|set(bits)):
        raise ValueError('physical cost and integer interface coordinates must be distinct')
    if any(F(str(c))<0 for c in cost.values()):raise ValueError('nonnegative group cost required')
    scale=F(1)if bound is None else F(str(bound))
    if scale<0:raise ValueError('nonnegative certified conditional lower bound required')
    row={}if bound is None else {j:-F(str(c))for j,c in cost.items()if c}
    for j in endpoints:row[j]=scale
    for j,col in ownership.items():row[col]=scale if j in S else-scale
    return dict(kind='infeasible'if bound is None else'cost',entries=row,
                rhs=scale*(len(endpoints)+len(S)-1),free_endpoints=len(endpoints))
