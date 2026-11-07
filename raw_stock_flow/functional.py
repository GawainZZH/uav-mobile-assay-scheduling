"""Fresh F/S entry points; every numerical input is derived from raw physics."""
import time
from fractions import Fraction as F
from review88_physical.evidence import plain
from review88_physical.physical import validate_fleet
from block_flow_solver.oracle import Oracle
from block_flow_solver.functional import dump
from block_flow_solver.restore import reconstruct
from block_price_propagation.oracle import SharedOracle
from price_tree_g5.warm import verified_source
from .initial import initialize,master
from .result import classify


def metadata(initial):
    return dict(currency_scale=initial.currency_scale,working_quantum=initial.quantum,working_floor=initial.floor,
        initializer_seconds=initial.seconds,historical_keys=0,old_answer_or_price_imported=False,
        initial_source_cost=None if initial.source is None else initial.source['check']['objective'],
        initial_family_count=initial.family_count,initial_stock_arc_count=initial.arc_count,
        initialization_records=initial.records,working_instance=plain(initial.domain.inst),raw_input_only=True)


def report(initial,r):
    D=initial.currency_scale;r['internal_objective']=r.get('objective');r['interval_internal']=r['interval']
    r['objective']=None if r['objective']is None else F(str(r['objective']))/D
    r['interval']=dict(r['interval'],lower=r['interval']['lower']/D,
        upper=None if r['interval']['upper']is None else r['interval']['upper']/D,
        quantum=r['interval']['quantum']/D,inherited_floor=r['interval']['inherited_floor']/D)
    fleet=r.get('solution',{}).get('fleet')if r.get('solution')else r.get('fleet')
    if fleet is not None:
        actual=validate_fleet(initial.original,fleet['group_plans'])
        if not actual['valid']or F(str(actual['objective']))!=r['objective']:raise ValueError('original currency/physical mismatch')
        r['physical_original']=actual
    r.update(metadata(initial));r['performance_claim']=False;return r


def flow(raw,deadline):
    began=time.perf_counter();initial=initialize(raw,deadline);m,x=master(initial);source=None
    if x is not None:
        source=verified_source(m,x,initial.source);source.solution.pop('inherited_true_incumbent',None);source.solution['raw_constructed_incumbent']=True
    base=Oracle(m.domain,deadline,key_seconds=5);shared=SharedOracle(m,base)
    try:
        from .engine import solve
        r=solve(m,shared.query,base.validate,reconstruct,deadline,floor=initial.floor,quantum=initial.quantum,source=source)
        from block_flow_solver.root import check
        feasible=False;error=None
        if r.get('master_vector')is not None:
            try:check(m,r['master_vector'],integer=True);feasible=True
            except ValueError as exc:error=repr(exc)
        r.update(matrix=dump(m),oracle_events=base.events,new_native_keys=len(base.cache),
            pricing_models={key:dict(instance=k.t.inst,crew_capacity=k.t.crew_capacity,cost_scale=k.t.cost_scale,
                matrix=k.matrix,row_keys=k.t.m.M.keys)for key,k in base.kernels.items()},
            shared_prices=[plain(p)for pool in shared.pools.values()for p in pool.prices.values()],
            propagation_events=shared.events,oracle_build_seconds=base.build_seconds,oracle_query_seconds=base.query_seconds,
            price_compilation_seconds=shared.compilation_seconds,price_propagation_seconds=shared.propagation_seconds,
            new_local_LP_calls=sum(k.lp_calls for k in base.kernels.values()),new_local_MIP_calls=sum(k.mip_calls for k in base.kernels.values()),
            final_master_point_feasible=feasible,final_master_point_error=error,total_seconds=time.perf_counter()-began)
        return report(initial,r)
    finally:base.close()


def strong(raw,deadline):
    began=time.perf_counter();initial=initialize(raw,deadline)
    upper=None if initial.source is None else initial.source['check']['objective'];interval=classify(None,upper,initial.floor,initial.quantum)
    if interval['closed']:
        return report(initial,dict(status='optimal',objective=upper,interval=interval,fleet=initial.source,
            matrix=None,native=None,native_solves=0,source_vector=None,physical_vector=None,
            analytic_rows=None,total_seconds=time.perf_counter()-began,reason='analytic bound equals independently validated source'))
    from .strong import build,complete,coordinates
    from boundary_context_master.native import solve_mip
    q=build(initial,deadline);native=solve_mip(q.adapter,deadline,initial=q.source_vector)
    if native['status']not in('optimal','timelimit','infeasible'):raise ValueError('unexpected fresh strong native status')
    fleet=None;actual=None
    if native['vector']is not None:fleet,actual=complete(q,native['vector']);upper=fleet['check']['objective']
    interval=classify(native['bound'],upper,initial.floor,initial.quantum)
    if native['status']=='infeasible'and upper is not None:raise ValueError('strong infeasible contradicts real source')
    if native['status']=='optimal'and not interval['closed']:raise ValueError('strong optimal lacks physical equality')
    return report(initial,dict(status='optimal'if interval['closed']else'infeasible'if native['status']=='infeasible'else'unresolved',
        objective=upper,interval=interval,fleet=fleet,native=native,native_solves=1,matrix=q.matrix,
        source_vector=q.source_vector,physical_vector=actual,analytic_rows=q.analytic_rows,feature_coordinates=coordinates(q),
        total_seconds=time.perf_counter()-began))
