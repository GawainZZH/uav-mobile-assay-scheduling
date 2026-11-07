"""One charged compact/monolithic functional workflow, with a common LP anchor."""
from fractions import Fraction as F
import math,time
import numpy as np

from referee_evidence21.algorithm import BoundedKernel
from service_separator25_joint.cases import split
from service_separator25_joint.solve import signature
from primal_transfer.native import FragmentProxy
from primal_transfer.incumbents import project_fleet,seed_kernel
from review88_physical.physical import validate_fleet
from review88_physical.milp import dump_matrix
from .common import describe,build_base,augment,embed_source,compact,complete_and_restore
from .matrix import encode
from .native import solve_lp,solve_mip
from .algorithm import solve_prepared


def run(inst,cuts,source,deadline,*,method,root_only=False,records=(),checkpoint=None):
    began=time.perf_counter();objects={};phases=[];partial_lower=0;common=None
    def measured(name,fn):
        t=time.perf_counter()
        if checkpoint:checkpoint(dict(stage=name,state='enter',elapsed=t-began))
        value=fn();row=dict(stage=name,state='complete',seconds=time.perf_counter()-t)
        if hasattr(value,'A'):row.update(shape=list(value.A.shape),nonzeros=value.A.nnz)
        phases.append(row)
        if checkpoint:checkpoint(row)
        return value
    try:
        checked=validate_fleet(inst,source['group_plans'])
        if not checked['valid']:raise ValueError('invalid common physical source')
        source=dict(source,check=checked);upper=F(checked['objective'])
        if partial_lower>upper:raise ValueError('invalid source objective')
        def kernels():
            for local,ids in split(inst,cuts):
                key=signature(local)
                if key in objects:continue
                k=BoundedKernel(local,deadline);p=k.mip
                k.mip=FragmentProxy(p.A,p.c,p.upper,p.lower,p.rhs,p.integer,reuse=True)
                objects[key]=k
        measured('local_matrices_and_stage_LP',kernels)
        fragments=project_fleet(inst,cuts,source)
        for part in fragments:seed_kernel(objects[part['kernel']],part)
        d=describe(inst,cuts,objects=objects,deadline=deadline)
        partial_lower=math.ceil(sum(d.floors))
        physical=measured('joint_strong_matrix',lambda:build_base(d,deadline))
        old_rows,old_columns=physical.A.shape
        old_lp=measured('joint_strong_LP',lambda:solve_lp(physical,deadline))
        if old_lp['status']!=0:raise TimeoutError('common original strong LP did not close')
        partial_lower=max(partial_lower,old_lp['integer_bound'])
        a=measured('joint_augmented_matrix',lambda:augment(d,physical,upper,records=records,deadline=deadline))
        source_vector,parts,actions=embed_source(d,a,source)
        new_lp=measured('joint_augmented_LP',lambda:solve_lp(a,deadline))
        if new_lp['status']!=0:raise TimeoutError('common augmented strong LP did not close')
        partial_lower=max(partial_lower,new_lp['integer_bound'])
        if new_lp['raw_objective']+1e-6<old_lp['raw_objective']:
            raise ValueError('augmentation weakened the common strong LP')
        anchor=max(old_lp['integer_bound'],new_lp['integer_bound'],math.ceil(sum(d.floors)))
        if anchor>upper:raise ValueError('common rational-price bound exceeds physical source')
        # The same rounded price bound is available in either representation.
        a.M.rhs[a.M.index['cb','anchor']]=-anchor
        a.anchor=F(anchor);a.A=a.M.sparse().tocsr();physical.A=a.A.tocsc()
        common=dict(complete=True,original_strong_LP=old_lp,augmented_strong_LP=new_lp,
            common_integer_anchor=anchor,source_objective=upper,
            original_shape=[old_rows,old_columns],augmented_shape=list(a.A.shape),
            original_rows_retained=True,source_vector=source_vector.tolist(),
            shared_cost_rows=tuple(records),source_is_only_primal=True)
        if method=='CB' and not root_only:
            answer=measured('compact_logical_search',lambda:solve_prepared(d,objects,source,anchor,deadline))
        elif method=='CB':
            master=measured('compact_root_matrix',lambda:compact(d,upper,records=records,anchor=anchor,deadline=deadline))
            hint=encode(master,actions,tuple(p['objective']for p in parts))
            native=measured('compact_root_integer_solve',lambda:solve_mip(master,deadline,initial=hint))
            bound=max(anchor,math.ceil(native['bound']-1e-6))if native['bound']is not None else anchor
            if bound>upper:raise ValueError('root master bound exceeds physical source')
            answer=dict(status='root_validated',optimal=bound==upper,objective=upper,lower_bound=bound,
                fleet=source,cuts=tuple(records),root_master=dict(native,
                    matrix=dump_matrix(master,np.asarray(master.M.cost),native['vector'])),
                root_shape=list(master.A.shape),root_nonzeros=master.A.nnz,local_calls=0)
        elif method=='PQ' and root_only:
            answer=dict(status='root_validated',optimal=anchor==upper,objective=upper,
                lower_bound=anchor,fleet=source,cuts=tuple(records),local_calls=0)
        elif method=='PQ':
            if anchor==upper:
                answer=dict(status='optimal',optimal=True,objective=upper,lower_bound=anchor,
                    fleet=source,cuts=tuple(records),native_integer_calls=0)
            else:
                native=measured('augmented_joint_integer_solve',lambda:solve_mip(a,deadline,initial=source_vector))
                bound=max(anchor,math.ceil(native['bound']-1e-6))if native['bound']is not None else anchor
                fleet=source;value=upper;completed=None
                if native['vector']is not None:
                    completed,restored=complete_and_restore(d,a,native['vector'])
                    if F(restored['check']['objective'])<value:fleet=restored;value=F(restored['check']['objective'])
                if bound>value:raise ValueError('monolithic numerical bound exceeds completed fleet')
                answer=dict(status='optimal'if bound==value else'feasible',optimal=bound==value,
                    objective=value,lower_bound=bound,fleet=fleet,cuts=tuple(records),
                    native_integer_calls=1,native=native,completed_vector=None if completed is None else completed.tolist(),
                    final_matrix=dump_matrix(a,np.asarray(a.M.cost),completed))
        else:raise ValueError('unknown compact prototype method')
        answer.update(method=method,root_only=root_only,common=common,phases=phases,
            total_workflow_seconds=time.perf_counter()-began,
            envelopes={key:k.envelope for key,k in objects.items()},
            local_models={key:dict(instance=k.inst,matrix=dump_matrix(k.m,k.cost,None))for key,k in objects.items()},
            local_native={key:dict(builds=k.mip.builds,resets=k.mip.context_resets,resumes=k.mip.resumes,
                optimize_seconds=k.mip.optimize_seconds)for key,k in objects.items()},
            full_resource_graph_builds=0,historical_bounds_imported=False,performance_claim=False)
        return answer
    except TimeoutError as exc:
        if partial_lower>upper:raise ValueError('partial common bound exceeds physical source')
        return dict(status='unresolved',optimal=False,reason=str(exc),objective=upper,
            lower_bound=partial_lower,fleet=source,method=method,root_only=root_only,
            cuts=tuple(records),common=common or dict(complete=False,common_integer_anchor=partial_lower),
            phases=phases,total_workflow_seconds=time.perf_counter()-began,
            full_resource_graph_builds=0,historical_bounds_imported=False,performance_claim=False)
    finally:
        for k in objects.values():k.mip.close();k.lp.clear()
