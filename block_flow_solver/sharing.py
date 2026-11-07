"""Transfer local lower/no-good rows and cost normalizations to strong MILP."""
from fractions import Fraction as F
from position_ownership_contract import conditional_row
from .matrix import row,add


def conditional_rows(cost,endpoints,ownership,subset,*,lower=None,upper=None,cap=None,infeasible=False):
    if infeasible:
        r=conditional_row(cost,endpoints,ownership,subset,None);r['validity']='local_infeasibility';return[r]
    result=[]
    if lower is not None:
        r=conditional_row(cost,endpoints,ownership,subset,lower);r['validity']='local_cost_lower_bound';result.append(r)
    if upper is not None:
        if cap is None or F(str(cap))<0:raise ValueError('original group-block finite cost cap required')
        U=F(str(upper));M=max(F(0),F(str(cap))-U)
        if U<0:raise ValueError('nonnegative feasible-plan upper required')
        # Validate the interface disjointness and canonical subset before forming the upper row.
        conditional_row(cost,endpoints,ownership,subset,0)
        e={j:F(str(c))for j,c in cost.items()if c};S=set(subset)
        for j in endpoints:add(e,j,M)
        for j,col in ownership.items():add(e,col,M if j in S else-M)
        r=row(e,U+M*(len(endpoints)+len(S)),False,'cost_normalization')
        r['validity']='optimality_preserving_normalization';result.append(r)
    return result
