"""Pure scheduling/evaluation contract; not an execution controller."""
from fractions import Fraction as F
from math import isfinite
from .cases import configurations


def arms():
    cases=[(c['name'],s)for c in configurations()for s in(11,29)]
    return[dict(case=f'RSC_{name}_s{seed}',configuration=name,calendar_seed=seed,method=m)
        for i,(name,seed)in enumerate(cases)for m in(('F','S')if i%2==0 else('S','F'))]


def assess(r,seconds):
    """Evaluate a previously physically checked result, without running a solver.

    Normal limits remain in PAR2. An impossible infeasible conclusion or any
    malformed result stops the batch. This does not replace physical replay.
    """
    fail=dict(eligible=False,**{'continue':False},par2=None,solved=False)
    if not isfinite(seconds)or seconds<0:return fail
    status=r.get('status')
    if status=='resource_limit':return dict(eligible=True,**{'continue':True},par2=420,solved=False)
    if status not in('optimal','unresolved'):return fail
    try:
        iv=r['interval'];lo=F(str(iv['lower']));hi=F(str(iv['upper']))
        if lo>hi or lo<0 or bool(iv['closed'])!=(lo==hi)or(status=='optimal')!=bool(iv['closed']):return fail
    except(KeyError,TypeError,ValueError,ZeroDivisionError):return fail
    solved=status=='optimal'and seconds<=210
    return dict(eligible=True,**{'continue':True},par2=seconds if solved else 420,solved=solved)
