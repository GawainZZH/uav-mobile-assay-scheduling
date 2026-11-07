"""Generic single SCIP tree with optional physical incumbent and objective grid."""
import time,math
from pyscipopt import Model,quicksum
from block_price_tree.engine import Handler
from block_price_tree.policy import Inspector
from price_tree_g5.warm import WarmInspector
from .result import classify


def solve(m,query,validate,reconstruct,deadline,*,floor,quantum,source=None):
    began=time.perf_counter();upper=None if source is None else source.solution['objective']
    initial=classify(None,upper,floor,quantum)
    if initial['closed']:
        return dict(status='optimal',objective=upper,interval=initial,solution=source.solution,master_vector=source.vector,
            native_status='not_started_analytic_closure',native_bound=None,native_solves=0,single_tree=True,nodes=0,
            initial_rows=len(m.rows),final_rows=len(m.rows),inserted_rows=0,pending=[],events=[],reason=None,
            source_stored=False,source_admissions=0,construction_seconds=0.,seconds=time.perf_counter()-began)
    model=Model();model.hideOutput();calls=0
    try:
        for k,v in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
                    ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.),('misc/allowstrongdualreds',False)]:model.setParam(k,v)
        vs=[model.addVar(name='v'+str(j),lb=0,ub=float(u),obj=float(c),vtype='B'if flag and u<=1 else'I'if flag else'C')
            for j,(u,c,flag)in enumerate(zip(m.upper,m.cost,m.integer))]
        initial_rows=len(m.rows)
        for r in m.rows:
            expr=quicksum(float(c)*vs[j]for j,c in r['entries'].items())
            model.addCons(expr==float(r['rhs'])if r['eq']else expr<=float(r['rhs']))
        args=(m,query,validate,reconstruct,deadline)
        policy=Inspector(*args)if source is None else WarmInspector(*args,source)
        h=Handler(policy,vs)
        for v in vs:model.markDoNotAggrVar(v);model.markDoNotMultaggrVar(v);model.addVarLocks(v,1,1)
        model.includeConshdlr(h,'raw_stock_cost','global shared local scheduling costs',enfopriority=1000000,
            chckpriority=-10,needscons=False,sepafreq=-1,maxprerounds=0)
        stored=False
        if source is not None:
            sol=model.createOrigSol()
            for v,x in zip(vs,source.vector):model.setSolVal(sol,v,float(x))
            if not model.checkSol(sol,printreason=False,completely=True,original=True):
                model.freeSol(sol);raise ValueError('SCIP rejected raw constructed physical source')
            stored=bool(model.addSol(sol,free=True))
        left=deadline-time.perf_counter()
        if left<=0:raise TimeoutError('raw stock construction exhausted allowance')
        model.setParam('limits/time',left);construction=time.perf_counter()-began;calls+=1;model.optimize()
        if policy.error:raise ValueError('raw stock handler: '+policy.error)
        status=str(model.getStatus());bound=float(model.getDualbound());bound=bound if math.isfinite(bound)and abs(bound)<1e18 else None
        upper=None if policy.best is None else policy.best['objective'];interval=classify(bound,upper,floor,quantum)
        if status=='infeasible'and upper is not None:raise ValueError('infeasible contradicts real source')
        if status not in('optimal','timelimit','userinterrupt','infeasible'):raise ValueError('unexpected raw SCIP status '+status)
        if status=='optimal'and not interval['closed']:raise ValueError('native optimal has no matching physical incumbent')
        return dict(status='optimal'if interval['closed']else'infeasible'if status=='infeasible'and not policy.reason else'unresolved',
            objective=upper,interval=interval,solution=policy.best,master_vector=policy.point,native_status=status,
            native_bound=bound,native_solves=calls,single_tree=True,nodes=int(model.getNNodes()),
            initial_rows=initial_rows,final_rows=len(m.rows),inserted_rows=policy.inserted,pending=list(policy.pending),
            checks=policy.checks,validated=policy.validated,events=policy.events,reason=policy.reason,
            source_stored=stored,source_admissions=getattr(policy,'source_admissions',0),construction_seconds=construction,
            seconds=time.perf_counter()-began)
    finally:model.freeProb()
