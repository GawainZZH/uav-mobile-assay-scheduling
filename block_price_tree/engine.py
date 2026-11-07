"""One SCIP tree; every cross-family learned row is inserted globally."""
import math,time
from fractions import Fraction as F
from pyscipopt import Model,Conshdlr,SCIP_RESULT,quicksum
from .policy import Inspector


class Handler(Conshdlr):
    def __init__(self,policy,variables):self.policy=policy;self.variables=variables;self.static_delegations=0

    def inspect(self,solution):
        try:state=self.policy.inspect([self.model.getSolVal(solution,v)for v in self.variables])
        except Exception as exc:self.policy.error=repr(exc);state='error'
        if state in('unresolved','error'):self.model.interruptSolve()
        return state

    def insert(self,r,index):
        expr=quicksum(float(c)*self.variables[j]for j,c in r['entries'].items());rhs=float(r['rhs'])
        self.model.addCons(expr==rhs if r['eq']else expr<=rhs,name='learned_'+str(index),
            initial=True,local=False,removable=False)

    def enforce(self,is_lp):
        state=self.inspect(None);before=self.policy.inserted
        try:self.policy.flush(self.insert)
        except Exception as exc:
            self.policy.error=repr(exc);self.model.interruptSolve();return {'result':SCIP_RESULT.INFEASIBLE}
        if self.policy.inserted>before:return {'result':SCIP_RESULT.CONSADDED}
        if state=='linear_infeasible':
            self.static_delegations+=1
            return {'result':SCIP_RESULT.INFEASIBLE if is_lp else SCIP_RESULT.SOLVELP}
        if state in('unresolved','error','cuts'):return {'result':SCIP_RESULT.INFEASIBLE}
        return {'result':SCIP_RESULT.FEASIBLE}

    def consenfolp(self,constraints,nusefulconss,solinfeasible):return self.enforce(True)
    def consenfops(self,constraints,nusefulconss,solinfeasible,objinfeasible):return self.enforce(False)
    def conscheck(self,constraints,solution,checkintegrality,checklprows,printreason,completely,**kwargs):
        return {'result':SCIP_RESULT.FEASIBLE if self.inspect(solution)=='feasible'else SCIP_RESULT.INFEASIBLE}
    def conslock(self,constraint,locktype,nlockspos,nlocksneg):
        for v in self.variables:self.model.addVarLocksType(v,locktype,nlockspos+nlocksneg,nlockspos+nlocksneg)


def solve(m,query,validate,reconstruct,deadline):
    began=time.perf_counter();model=Model();model.hideOutput()
    try:
        for k,v in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
                    ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.),('misc/allowstrongdualreds',False)]:model.setParam(k,v)
        vs=[model.addVar(name='v'+str(j),lb=0,ub=float(u),obj=float(c),vtype='B'if flag and u<=1 else'I'if flag else'C')
            for j,(u,c,flag)in enumerate(zip(m.upper,m.cost,m.integer))]
        initial=len(m.rows)
        for r in m.rows:
            expr=quicksum(float(c)*vs[j]for j,c in r['entries'].items())
            model.addCons(expr==float(r['rhs'])if r['eq']else expr<=float(r['rhs']))
        policy=Inspector(m,query,validate,reconstruct,deadline);h=Handler(policy,vs)
        for v in vs:model.markDoNotAggrVar(v);model.markDoNotMultaggrVar(v);model.addVarLocks(v,1,1)
        model.includeConshdlr(h,'shared_stock_family_cost','global shared local scheduling costs',enfopriority=1000000,
            chckpriority=-10,needscons=False,sepafreq=-1,maxprerounds=0)
        left=deadline-time.perf_counter()
        if left<=0:raise TimeoutError('stock model construction exhausted allowance')
        model.setParam('limits/time',left);construction=time.perf_counter()-began;model.optimize()
        if policy.error:raise ValueError('shared stock handler: '+policy.error)
        status=str(model.getStatus());bound=float(model.getDualbound());bound=bound if math.isfinite(bound)and abs(bound)<1e18 else None
        value=None if policy.best is None else F(str(policy.best['objective']))
        closed=status=='optimal'and value is not None and bound is not None and abs(bound-float(value))<=1e-6 and not policy.reason
        if status=='optimal'and not closed:raise ValueError(('native optimal lacks same-cost physical incumbent',bound,value,len(policy.pending),policy.reason))
        return dict(status='optimal'if closed else'infeasible'if status=='infeasible'and not policy.reason else'unresolved',
            objective=value,lower_bound=value if closed else bound,solution=policy.best,master_vector=policy.point,
            native_status=status,native_bound=bound,native_solves=1,single_tree=True,nodes=int(model.getNNodes()),
            columns=len(vs),initial_rows=initial,final_rows=len(m.rows),inserted_rows=policy.inserted,
            checks=policy.checks,validated=policy.validated,static_delegations=h.static_delegations,
            events=policy.events,pending=list(policy.pending),reason=policy.reason,
            construction_seconds=construction,seconds=time.perf_counter()-began,
            closure_basis='valid_relaxation_bound_and_verified_original_solution'if closed else None,
            proof_scope='native_numerical_branch_and_bound_with_physical_validation',performance_claim=False)
    finally:model.freeProb()
