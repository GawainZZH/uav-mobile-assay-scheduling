"""Complete known-family information for a future fair strong baseline.

This is separate from the already frozen 239-bound delta replay.
"""
from fractions import Fraction as F
from math import isfinite,lcm
from block_flow_solver.sharing import conditional_rows
from block_flow_solver.matrix import activity
from vehicle_root_ablation.model import append_rows


def full_rows(p,bands):
    result=[];K=len(p.d.stages)
    for i,b in enumerate(bands):
        k=b['stage'];local,ids=p.d.stages[k][:2];S=tuple(ids.index(j)for j in b['subset']);a=b['start'];q=b['end']
        if not 0<=a<len(local.pads)or not 0<=q<len(local.pads)or k==0 and a!=0 or k==K-1 and q!=0:
            raise ValueError('original endpoint domain required')
        if b['feasible']and b.get('plan')is None:raise ValueError('upper normalization requires a previously validated plan')
        if b['feasible']and b['infeasible']:raise ValueError('contradictory family feasibility')
        cap=local.H*local.road_cost+len(ids)*(local.flight_fixed+local.air_cost*local.endurance)+local.H*sum((j.weight for j in local.jobs),F(0))
        for g in range(p.inst.groups):
            ends=(()if k==0 else(p.y[k,g,a],))+(()if k==K-1 else(p.y[k+1,g,q],))
            rows=conditional_rows(p.costs[k,g],ends,{j:p.s[k,g,j]for j in range(len(ids))},S,
                lower=None if b['infeasible']else F(str(b['lower'])),
                upper=F(str(b['upper']))if b['feasible']else None,cap=cap,infeasible=b['infeasible'])
            for r in rows:
                den=lcm(*(F(str(v)).denominator for v in list(r['entries'].values())+[r['rhs']]))
                r['entries']={j:F(str(v))*den for j,v in r['entries'].items()};r['rhs']=F(str(r['rhs']))*den
                r.update(family=i,stage=k,group=g,source_scope=b.get('proof_scope','retained_family_bound'))
            result.extend(rows)
    return result


def export(base,rows,source):
    matrix=append_rows(base,rows)
    if source is None:return dict(matrix=matrix,new_rows=rows,source_vector=None,warmstart_valid=False,source_violations=[])
    if len(source)!=base['shape'][1]or not all(isfinite(float(x))for x in source):raise ValueError('finite source vector of original dimension required')
    violations=[i for i,r in enumerate(rows)if activity(r,source)>1e-7]
    if any(rows[i]['kind']!='cost_normalization'for i in violations):raise ValueError('new lower/no-good contradicts original true source')
    return dict(matrix=matrix,new_rows=rows,source_vector=source if not violations else None,
        warmstart_valid=not violations,source_violations=violations,
        source_rejection_reason='old calendar violates an optimality-preserving cost normalization'if violations else None)
