"""Same raw-input/candidate information for both cold methods."""
import time
from fractions import Fraction as F
from unittest.mock import patch
from review88_physical.evidence import instance
from review88_physical.physical import validate_fleet
from stock_scale.initial import prepare,compact_rows
from stock_scale.contract import assess
from raw_stock_flow.result import classify
from raw_stock_flow.functional import report
from .projection import certify


def checkpoint(raw,fleet,lower=0,reason='verified common candidate'):
    actual=validate_fleet(instance(raw['instance']),fleet['group_plans'])
    if not actual['valid']:raise ValueError('checkpoint is not physically feasible')
    upper=F(str(actual['objective']));lower=F(str(lower))
    if not 0<=lower<=upper:raise ValueError('checkpoint bounds contradict')
    return dict(status='resource_limit',objective=upper,interval=dict(lower=lower,upper=upper,closed=lower==upper),
        fleet=dict(group_plans=fleet['group_plans'],check=actual),reason=reason,performance_claim=False)


def check_result(raw,r):
    if r.get('status')=='resource_limit'and r.get('interval',{}).get('upper')is None:
        iv=r.get('interval',{})
        if r.get('objective')is not None or iv.get('closed')or F(str(iv.get('lower',-1)))!=0:
            raise ValueError('unwitnessed resource result must have only the nonnegative floor')
        return dict(valid=True,objective=None,scope='no upper witness published within allowance')
    fleet=r.get('fleet')or(r.get('solution')or{}).get('fleet')
    if fleet is None:raise ValueError('every fixed case has a physical source; missing upper witness')
    actual=validate_fleet(instance(raw['instance']),fleet['group_plans'])
    iv=r['interval'];lo=F(str(iv['lower']));hi=F(str(iv['upper']))
    if not actual['valid']or hi!=F(str(actual['objective']))or hi!=F(str(r['objective']))or not 0<=lo<=hi:
        raise ValueError('result bound/physical objective differs')
    if bool(iv['closed'])!=(lo==hi):raise ValueError('closure flag contradicts interval')
    if not assess(r,0)['eligible']:raise ValueError('result status contradicts evaluation contract')
    return actual


def solve_record(record,method,deadline,save):
    began=time.perf_counter();raw=record['raw'];w=record['candidate']
    save(checkpoint(raw,w))
    a=prepare(raw,w,method,deadline)
    save(checkpoint(raw,a.source,a.floor/a.currency_scale,'same normalized source and analytic lower'))
    if method=='F':
        proof=certify(a)
        from raw_stock_flow.functional import flow
        # Reuse the unchanged single tree, pricing, learning and physical checks.
        # Only its initialization is supplied from this common raw preparation.
        with patch('raw_stock_flow.functional.initialize',return_value=a):r=flow(raw,deadline)
        r['compact_projection_support']=proof
    elif method=='S':r=strong(a,deadline)
    else:raise ValueError('fixed F/S methods')
    r.update(compact_fallback=a.compact_fallback,entry_seconds=time.perf_counter()-began)
    check_result(raw,r);return r


def strong(a,deadline):
    from raw_stock_flow.strong import build,append,complete,coordinates
    from cost_support.model import numeric_adapter
    from boundary_context_master.matrix import check_vector
    from boundary_context_master.native import solve_mip
    upper=a.source['check']['objective'];iv=classify(None,upper,a.floor,a.quantum)
    if iv['closed']:
        return report(a,dict(status='optimal',objective=upper,interval=iv,fleet=a.source,
            native_solves=0,matrix=None,reason='analytic closure'))
    q=build(a,deadline);rows=compact_rows(q,a);q.matrix=append(q.matrix,rows);q.adapter=numeric_adapter(q.matrix)
    check_vector(q.adapter,q.source_vector)
    native=solve_mip(q.adapter,deadline,initial=q.source_vector)
    if native['status']not in('optimal','timelimit','memlimit','userinterrupt','infeasible'):
        raise ValueError('unexpected physical native status '+native['status'])
    if native['status']=='infeasible':raise ValueError('native infeasible contradicts true common source')
    fleet=a.source;actual=q.source_vector
    if native['vector']is not None:
        candidate,vector=complete(q,native['vector'])
        if candidate['check']['objective']<=upper:fleet=candidate;actual=vector;upper=candidate['check']['objective']
    iv=classify(native['bound'],upper,a.floor,a.quantum)
    if native['status']=='optimal'and not iv['closed']:raise ValueError('native optimum lacks physical equality')
    return report(a,dict(status='optimal'if iv['closed']else'unresolved',objective=upper,interval=iv,
        fleet=fleet,native=native,native_solves=1,matrix=q.matrix,source_vector=q.source_vector,
        physical_vector=actual,analytic_rows=q.analytic_rows,compact_rows=rows,feature_coordinates=coordinates(q)))
