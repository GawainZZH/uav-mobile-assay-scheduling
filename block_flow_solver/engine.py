"""One native SCIP tree; physical incumbents require closed local families."""
import math,time
from fractions import Fraction as F
from pyscipopt import Model,Conshdlr,SCIP_RESULT,SCIP_STAGE,quicksum
from .root import check
from .matrix import band_rows,activity
from .bands import exact,update


class Handler(Conshdlr):
    def __init__(self,m,v,query,validate,reconstruct,deadline):
        self.m=m;self.v=v;self.query=query;self.validate=validate;self.reconstruct=reconstruct;self.deadline=deadline
        self.queried=set();self.pending=[];self.events=[];self.error=None;self.reason=None;self.best=None;self.point=None
        self.checks=0;self.validated=0;self.static_delegations=0;self.added=0

    def inspect(self,solution):
        self.checks+=1
        try:
            if time.perf_counter()>=self.deadline:
                self.reason='global deadline';self.model.interruptSolve();return 'unresolved'
            x=[self.model.getSolVal(solution,v)for v in self.v]
            try:check(self.m,x)
            except ValueError:return 'linear_infeasible'
            if any(abs(x[j]-round(x[j]))>1e-6 for j,v in enumerate(self.m.integer)if v):return 'fractional'
            selected=[i for i,js in self.m.family_columns.items()if sum(x[j]for j in js)>.5]
            learned=False;unresolved=False
            for i in selected:
                b=self.m.bands[i]
                if exact(b):continue
                if i in self.queried:unresolved=True;continue
                self.queried.add(i);old=dict(b);r=self.query(i);update(b,r,lambda p:self.validate(i,p))
                self.events.append(dict(family=i,record=r))
                if b['infeasible']or b['lower']!=old['lower']or b['upper']!=old['upper']:
                    rows=band_rows(self.m,i);self.m.rows.extend(rows);self.pending.extend(rows);learned=True
                if not exact(b):unresolved=True
            if learned:return 'cuts'
            if unresolved:
                self.reason='selected family cannot be certified within once-only query budget';self.model.interruptSolve();return 'unresolved'
            result=self.reconstruct(self.m,x);value=F(str(result['objective']))
            objective=sum(float(c)*v for c,v in zip(self.m.cost,x))
            if abs(float(value)-objective)>1e-6:raise ValueError('physical incumbent cost differs from master objective')
            self.validated+=1
            if self.best is None or value<F(str(self.best['objective'])):self.best=result;self.point=x
            return 'feasible'
        except TimeoutError as exc:
            self.reason=str(exc);self.model.interruptSolve();return 'unresolved'
        except Exception as exc:
            self.error=repr(exc);self.model.interruptSolve();return 'error'

    def enforce(self,is_lp):
        state=self.inspect(None);added=0
        for r in self.pending:
            self.model.addCons(quicksum(float(c)*self.v[j]for j,c in r['entries'].items())<=float(r['rhs']),
                initial=True,local=False,removable=False);added+=1
        self.pending.clear();self.added+=added
        if added:return {'result':SCIP_RESULT.CONSADDED}
        if state=='linear_infeasible':
            self.static_delegations+=1
            return {'result':SCIP_RESULT.INFEASIBLE if is_lp else SCIP_RESULT.SOLVELP}
        if state in('error','unresolved','cuts'):return {'result':SCIP_RESULT.INFEASIBLE}
        # Native integrality handler performs the actual flow branching.
        return {'result':SCIP_RESULT.FEASIBLE}

    def consenfolp(self,constraints,nusefulconss,solinfeasible):return self.enforce(True)
    def consenfops(self,constraints,nusefulconss,solinfeasible,objinfeasible):return self.enforce(False)
    def conscheck(self,constraints,solution,checkintegrality,checklprows,printreason,completely,**kwargs):
        return {'result':SCIP_RESULT.FEASIBLE if self.inspect(solution)=='feasible'else SCIP_RESULT.INFEASIBLE}
    def conslock(self,constraint,locktype,nlockspos,nlocksneg):
        for v in self.v:self.model.addVarLocksType(v,locktype,nlockspos+nlocksneg,nlockspos+nlocksneg)


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
        h=Handler(m,vs,query,validate,reconstruct,deadline)
        for v in vs:model.markDoNotAggrVar(v);model.markDoNotMultaggrVar(v);model.addVarLocks(v,1,1)
        model.includeConshdlr(h,'stock_family_cost','verified local scheduling costs',enfopriority=1000000,
            chckpriority=-10,needscons=False,sepafreq=-1,maxprerounds=0)
        left=deadline-time.perf_counter()
        if left<=0:raise TimeoutError('stock model construction exhausted allowance')
        model.setParam('limits/time',left);construction=time.perf_counter()-began;model.optimize()
        if h.error:raise ValueError('stock cost handler: '+h.error)
        status=str(model.getStatus());bound=float(model.getDualbound());bound=bound if math.isfinite(bound)and abs(bound)<1e18 else None
        value=None if h.best is None else F(str(h.best['objective']))
        # Pending global strengthening cannot invalidate a lower bound from the
        # current outer relaxation. A separately verified fleet may close it.
        closed=status=='optimal'and value is not None and bound is not None and abs(bound-float(value))<=1e-6 and not h.reason
        if status=='optimal'and not closed:raise ValueError(('native optimal without a verified same-cost original incumbent',bound,value,len(h.pending),h.reason,h.checks,h.validated))
        return dict(status='optimal'if closed else'infeasible'if status=='infeasible'and not h.reason else'unresolved',
            objective=value,lower_bound=value if closed else bound,solution=h.best,master_vector=h.point,
            native_status=status,native_bound=bound,native_solves=1,single_tree=True,nodes=int(model.getNNodes()),
            columns=len(vs),initial_rows=initial,final_rows=len(m.rows),inserted_rows=h.added,checks=h.checks,
            validated=h.validated,static_delegations=h.static_delegations,events=h.events,pending=h.pending,reason=h.reason,
            construction_seconds=construction,seconds=time.perf_counter()-began,
            closure_basis='valid_relaxation_bound_and_verified_original_solution'if closed else None,
            proof_scope='native_numerical_branch_and_bound_with_physical_validation',performance_claim=False)
    finally:model.freeProb()
