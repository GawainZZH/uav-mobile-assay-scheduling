"""Functional entry points. Timing is recorded, never a fair performance result."""
from fractions import Fraction as F
from review88_physical.evidence import instance
from .domain import from_instance
from .matrix import build
from .oracle import Oracle
from .root import close_root
from .engine import solve
from .restore import reconstruct


def dump(m):
    return dict(cost=m.cost,upper=m.upper,integer=m.integer,rows=m.rows,flow_columns=m.flow_columns,
        family_columns=list(m.family_columns.items()),z=list(m.z.items()),theta=list(m.theta.items()),u=list(m.u.items()),
        families=m.bands,arcs=m.domain.arcs,tasks=m.domain.tasks)


def flow(raw,deadline):
    inst=instance(raw['instance']);d=from_instance(inst,raw['cuts']);m=build(d);o=Oracle(d,deadline,key_seconds=5)
    try:
        root=close_root(m,o.query,o.validate,deadline)
        result=dict(root=root)
        if root['status']not in('closed','numerical_support_closed'):
            result.update(status='unresolved',reason='root support not closed; no native tree started')
        else:
            result.update(solve(m,o.query,o.validate,reconstruct,deadline))
        result.update(matrix=dump(m),oracle_events=o.events,oracle_build_seconds=o.build_seconds,
            oracle_query_seconds=o.query_seconds,physical_kernels=len(o.kernels),new_native_keys=len(o.cache),
            pricing_models={key:dict(instance=k.t.inst,crew_capacity=k.t.crew_capacity,cost_scale=k.t.cost_scale,
                matrix=k.matrix,row_keys=k.t.m.M.keys)for key,k in o.kernels.items()},
            raw_input_only=True,old_answer_or_price_imported=False,
            strict_original_optimality_certificate=False,
            matrix_LP_price_does_not_upgrade_native_G1_integer_proofs=True,performance_claim=False)
        return result
    finally:o.close()


def strong_control(core,flow_result,deadline):
    """Retain the old strong matrix, add every acquired applicable local row."""
    from position_ownership.model import build as old_build
    from position_ownership.restore import strong_variant,complete_real_fleet
    from vehicle_root_ablation.model import append_rows,adapter
    from boundary_context_master.matrix import check_vector
    from boundary_context_master.native import solve_mip
    from .sharing import conditional_rows
    if flow_result['status']!='optimal':raise ValueError('functional control requires its verified flow source')
    inst=instance(core['instance']);cuts=tuple(core['cuts']);source=flow_result['solution']['fleet']
    p=old_build(inst,cuts,core['matrix'],core['physical_shape'],core['original_augmented_rows'],source,deadline=deadline)
    p=strong_variant(p,keep_positions_integer=True);rows=[];normalizations=0
    for f in flow_result['matrix']['families']:
        k=f['stage'];local,ids,_=p.d.stages[k];S=tuple(ids.index(j)for j in f['subset'])
        # Unqueried conditional caps are not feasible-plan normalizations.
        if not f['infeasible']and not f['feasible']and not F(str(f['lower'])):continue
        cap=local.road_cost*local.H+len(ids)*(local.flight_fixed+local.air_cost*local.endurance)+local.H*sum((j.weight for j in local.jobs),F(0))
        for g in range(inst.groups):
            ends=(()if k==0 else(p.y[k,g,f['start']],))+(()if k==len(p.d.stages)-1 else(p.y[k+1,g,f['end']],))
            new=conditional_rows(p.costs[k,g],ends,{j:p.s[k,g,j]for j in range(len(ids))},S,
                lower=None if f['infeasible']else f['lower'],upper=f['upper']if f['feasible']else None,
                cap=cap,infeasible=f['infeasible'])
            normalizations+=sum(r['kind']=='cost_normalization'for r in new);rows.extend(new)
    p.matrix=append_rows(p.matrix,rows);p.adapter=adapter(p.matrix)
    # The new flow source already uses the minimum certified plan of each selected
    # family; verify it again against the synchronized normalizations before solve.
    check_vector(p.adapter,p.source_vector)
    native=solve_mip(p.adapter,deadline,initial=p.source_vector)
    if native['vector']is None:raise ValueError('strong control lost its valid original source')
    fleet,actual=complete_real_fleet(p,native['vector']);cost=F(str(fleet['check']['objective']))
    closed=native['status']=='optimal'and abs(float(native['bound'])-float(cost))<=1e-6
    if not closed:raise ValueError('strong functional control did not close')
    return dict(status='optimal',objective=cost,lower_bound=cost,fleet=fleet,native=native,
        synchronized_rows=rows,synchronized_normalizations=normalizations,source_vector=p.source_vector,
        physical_vector=actual,matrix=p.matrix,source_is_new_flow_result=True,
        source_generation_is_not_free=True,all_old_core_rows_retained=True,
        proof_scope='native_numerical_branch_and_bound_with_physical_validation',performance_claim=False)
