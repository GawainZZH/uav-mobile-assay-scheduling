"""Raw-input single-tree run and complete applicable-information control."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
from dataclasses import replace
from review88_physical.evidence import instance,plain
from service_separator25_joint.solve import signature
from block_flow_solver.domain import from_instance
from block_flow_solver.matrix import build
from block_flow_solver.oracle import Oracle
from block_flow_solver.restore import reconstruct
from block_flow_solver.functional import dump
from block_price_propagation.oracle import SharedOracle
from block_price_propagation.sharing import affine_row
from block_price_propagation.prices import AffinePrice
from .engine import solve


def affine_rows(p,prices):
    rows=[];K=len(p.d.stages)
    for raw in prices:
        price=AffinePrice(**raw)if isinstance(raw,dict)else raw
        for k,stage in enumerate(p.d.stages):
            local,ids=stage[:2]
            if price.instance_signature!=signature(replace(local,groups=1))or price.crew_capacity!=p.inst.drones:continue
            for g in range(p.inst.groups):
                start=0 if k==0 else {a:p.y[k,g,a]for a in range(len(local.pads))}
                end=0 if k==K-1 else {a:p.y[k+1,g,a]for a in range(len(local.pads))}
                r=affine_row(price,p.costs[k,g],{j:p.s[k,g,j]for j in range(len(ids))},start,end)
                r.update(stage=k,group=g);rows.append(r)
    return rows


def flow(raw,deadline):
    inst=instance(raw['instance']);d=from_instance(inst,raw['cuts']);m=build(d);base=Oracle(d,deadline,key_seconds=5)
    shared=SharedOracle(m,base)
    try:
        result=solve(m,shared.query,base.validate,reconstruct,deadline)
        result.update(matrix=dump(m),oracle_events=base.events,oracle_build_seconds=base.build_seconds,
            oracle_query_seconds=base.query_seconds,physical_kernels=len(base.kernels),new_native_keys=len(base.cache),
            pricing_models={key:dict(instance=k.t.inst,crew_capacity=k.t.crew_capacity,cost_scale=k.t.cost_scale,
                matrix=k.matrix,row_keys=k.t.m.M.keys)for key,k in base.kernels.items()},
            shared_prices=[plain(p)for pool in shared.pools.values()for p in pool.prices.values()],
            propagation_events=shared.events,propagation_changes=shared.changes,
            compilation_seconds=shared.compilation_seconds,propagation_seconds=shared.propagation_seconds,
            raw_input_only=True,old_answer_or_price_imported=False,root_support_prerequisite=False,
            strict_original_optimality_certificate=False,performance_claim=False)
        return result
    finally:base.close()


def strong_control(core,flow_result,deadline):
    from position_ownership.model import build as old_build
    from position_ownership.restore import strong_variant,complete_real_fleet
    from vehicle_root_ablation.model import append_rows,adapter
    from boundary_context_master.matrix import check_vector
    from boundary_context_master.native import solve_mip
    from block_price_propagation.full_sync import full_rows
    if flow_result['status']!='optimal':raise ValueError('new physical source must close before functional control')
    inst=instance(core['instance']);source=flow_result['solution']['fleet']
    p=old_build(inst,tuple(core['cuts']),core['matrix'],core['physical_shape'],core['original_augmented_rows'],source,deadline=deadline)
    p=strong_variant(p,keep_positions_integer=True);base=p.matrix
    families=full_rows(p,flow_result['matrix']['families']);affine=affine_rows(p,flow_result['shared_prices'])
    p.matrix=append_rows(base,families+affine);p.adapter=adapter(p.matrix)
    check_vector(p.adapter,p.source_vector)
    native=solve_mip(p.adapter,deadline,initial=p.source_vector)
    if native['vector']is None:raise ValueError('control lost the verified new physical source')
    fleet,actual=complete_real_fleet(p,native['vector']);cost=F(str(fleet['check']['objective']))
    if native['status']!='optimal'or abs(float(native['bound'])-float(cost))>1e-6:
        raise ValueError('functional strong control did not close')
    return dict(status='optimal',objective=cost,lower_bound=cost,fleet=fleet,native=native,
        synchronized_rows=families+affine,family_rows=len(families),affine_rows=len(affine),
        source_vector=p.source_vector,physical_vector=actual,matrix=p.matrix,old_rows=base['shape'][0],
        source_is_new_flow_result=True,source_generation_is_not_free=True,all_old_core_rows_retained=True,
        proof_scope='native_numerical_branch_and_bound_with_physical_validation',performance_claim=False)
