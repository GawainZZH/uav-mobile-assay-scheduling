"""A conditional upper cap is never a primal-feasibility certificate."""
from fractions import Fraction as F


def exact(b):
    return bool(not b['infeasible'] and b['feasible'] and b['plan'] is not None and b['lower']==b['upper'])


def update(b,record,validate):
    n=dict(b)
    if record.get('status')=='infeasible':
        if not record.get('proof_scope'):raise ValueError('infeasibility requires its proof scope')
        if n['feasible']:raise ValueError('infeasibility contradicts a verified plan')
        n.update(infeasible=True,exact=False,proof_scope=record['proof_scope'])
    else:
        lb=record.get('lower_bound')
        if lb is not None:n['lower']=max(n['lower'],F(str(lb)))
        if record.get('objective')is not None:
            if record.get('plan')is None:raise ValueError('upper admission requires a plan')
            value=F(str(validate(record['plan'])))
            if value!=F(str(record['objective'])):raise ValueError('upper cost differs from physical validation')
            if not n['feasible']and value>n['upper']:raise ValueError('first admitted plan exceeds the declared conditional upper; no matching primal bound')
            if not n['feasible']or value<n['upper']:n.update(upper=min(n['upper'],value),plan=record['plan'],feasible=True)
            if record.get('status')=='optimal'and n['lower']!=value:raise ValueError('optimal record has no matching lower bound')
        if n['lower']<0 or n['lower']>n['upper']:raise ValueError('cost interval contradiction')
        n['exact']=exact(n)
        if record.get('closure_basis'):n['proof_scope']=record['closure_basis']
    b.clear();b.update(n)
    return b
