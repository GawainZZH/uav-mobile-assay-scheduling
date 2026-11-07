"""Reusable one-group pricing matrix, actual cover RHS and certified bounds."""
from fractions import Fraction as F
import math,time
import numpy as np

from fleet_sample_assay_certified.master import highs,pass_lp,native_run,check
from service_separator25_joint.dual import price
from primal_transfer.native import FragmentProxy
from review88_physical.milp import dump_matrix
from .model import build_template,context_rhs,elementary,restore_subset,residual


class SubsetKernel:
    def __init__(self,inst,deadline,*,crew_capacity=None):
        began=time.perf_counter();self.t=t=build_template(inst,deadline,crew_capacity=crew_capacity);self.seen=set()
        low=np.where(t.eq,t.base,-np.inf)
        self.lp=highs();pass_lp(self.lp,t.A,t.scaled_cost,t.upper,low,t.base)
        self.mip=FragmentProxy(t.A,t.scaled_cost,t.upper,low,t.base,t.integer,reuse=True)
        self.matrix=dump_matrix(t.m,t.scaled_cost,None);self.build_seconds=time.perf_counter()-began
        self.lp_calls=0;self.mip_calls=0

    def close(self):self.mip.close();self.lp.clear()

    def query(self,subset,start,end,deadline):
        began=time.perf_counter();t=self.t;subset=tuple(subset);key=(subset,start,end)
        if key in self.seen:raise ValueError('once-only subset query cannot repeat')
        self.seen.add(key);rhs=context_rhs(t,subset,start,end);lo=np.where(t.eq,rhs,-np.inf)
        lower=elementary(t,subset,start,end)
        if lower is None:return dict(status='infeasible',optimal=False,subset=subset,start=start,end=end,
            objective=None,lower_bound=None,plan=None,reason='no independent flight event')
        result=dict(status='unresolved',optimal=False,subset=subset,start=start,end=end,
            actual_rhs=rhs.tolist(),lower_bound=lower,elementary_bound=lower,
            objective=None,plan=None,physical=None,vector=None,native_vector=None,
            LP=None,native=None,cost_scale=t.cost_scale)
        ids=np.arange(len(rhs),dtype=np.int32)
        check(self.lp.changeRowsBounds(len(rhs),ids,lo,rhs))
        try:
            q=time.perf_counter();code,sol,info,msg=native_run(self.lp,deadline);self.lp_calls+=1
            result['LP']=dict(status=code,message=msg,seconds=time.perf_counter()-q)
            if code==2:
                if sol.value_valid:
                    try:residual(t,sol.col_value,rhs,relax=True)
                    except ValueError:pass
                    else:raise ValueError('LP infeasibility contradicts a feasible primal')
                result.update(status='infeasible',lower_bound=None,proof_scope='native_numerical_LP_infeasibility')
                return result
            if code==0 and sol.value_valid and sol.dual_valid:
                residual(t,sol.col_value,rhs,relax=True)
                p=price(t.A,t.scaled_cost,t.upper,t.eq,sol.row_dual);exact=p.evaluate(rhs)
                bound=F(math.ceil(exact),t.cost_scale);lower=max(lower,bound)
                result['LP'].update(vector=list(sol.col_value),raw_scaled_objective=float(t.scaled_cost@np.asarray(sol.col_value)),
                    price=dict(weights=p.weights.tolist(),constant=p.constant,scale=p.scale),
                    exact_dual=exact,bound=bound,scale=t.cost_scale)
            check(self.mip.changeRowsBounds(len(rhs),ids,lo,rhs))
            q=time.perf_counter();code,sol,info,msg=native_run(self.mip,deadline);self.mip_calls+=1
            native_bound=float(info.mip_dual_bound)
            result['native']=dict(status=code,message=msg,seconds=time.perf_counter()-q,
                scaled_bound=native_bound if math.isfinite(native_bound)and abs(native_bound)<1e18 else None,
                nodes=int(info.mip_node_count))
            if code==2:
                if sol.value_valid:raise ValueError('MIP infeasibility contradicts primal')
                result.update(status='infeasible',lower_bound=None,proof_scope='native_numerical_integer_infeasibility')
                return result
            if math.isfinite(native_bound)and abs(native_bound)<1e18:
                lower=max(lower,F(math.ceil(native_bound-1e-6),t.cost_scale))
            result['lower_bound']=lower
            if sol.value_valid:
                raw=np.asarray(sol.col_value);residual(t,raw,rhs)
                vector,plan,physical,note=restore_subset(t,raw,subset,start,end)
                upper=F(physical['objective'])
                if lower>upper:raise ValueError('subset bound exceeds physical plan')
                result.update(objective=upper,plan=plan,physical=physical,vector=vector.tolist(),
                    native_vector=raw.tolist(),normalization=note)
                if lower==upper:
                    lp_bound=result['LP'].get('bound',0)
                    basis=('elementary_bound_and_physical_plan'if result['elementary_bound']>=upper else
                           'rational_LP_price_and_physical_plan'if lp_bound>=upper else
                           'native_numerical_integer_bound_and_physical_plan')
                    result.update(status='optimal',optimal=True,closure_basis=
                        basis)
        except TimeoutError as exc:
            result.update(reason=str(exc),lower_bound=lower)
        finally:
            result['seconds']=time.perf_counter()-began
        return result
