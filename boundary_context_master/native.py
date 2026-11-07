"""Timed LP prices and single-thread integer solves for the finite prototype."""
import math,time
from fractions import Fraction
import numpy as np

from .matrix import check_vector


def solve_lp(data,deadline):
    from fleet_sample_assay_certified.master import highs,pass_lp,native_run
    from service_separator25_joint.dual import price
    from review88_physical.milp import dump_matrix
    began=time.perf_counter();M=data.M;c=np.asarray(M.cost,dtype=float)
    upper=np.asarray(M.upper,dtype=float);rhs=np.asarray(M.rhs,dtype=float);eq=np.asarray(M.eq,dtype=bool)
    h=highs()
    try:
        pass_lp(h,data.A,c,upper,np.where(eq,rhs,-np.inf),rhs)
        status,sol,info,msg=native_run(h,deadline)
        out=dict(status=status,message=msg,seconds=time.perf_counter()-began,
                 raw_objective=None,exact_bound=None,integer_bound=None)
        if status!=0 or not sol.value_valid or not sol.dual_valid:return out
        check_vector(data,sol.col_value,integer=False)
        p=price(data.A,c,upper,eq,sol.row_dual);bound=p.evaluate(rhs)
        out.update(raw_objective=float(c@np.asarray(sol.col_value)),exact_bound=bound,
            integer_bound=math.ceil(bound),matrix=dump_matrix(data,c,sol.col_value),
            price=dict(weights=p.weights.tolist(),constant=p.constant,scale=p.scale))
        out['seconds']=time.perf_counter()-began
        return out
    finally:
        h.clear()


def solve_mip(data,deadline,initial=None):
    from pyscipopt import Model,quicksum
    began=time.perf_counter();model=Model();model.hideOutput();M=data.M
    try:
        for k,v in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
                    ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.)]:model.setParam(k,v)
        variables=[model.addVar(name='c'+str(j),vtype='B'if M.integer[j]and M.upper[j]<=1 else'I'if M.integer[j]else'C',
            lb=0,ub=float(M.upper[j]),obj=float(M.cost[j]))for j in range(len(M.cost))]
        A=data.A.tocsr()
        for i in range(A.shape[0]):
            lo,hi=A.indptr[i:i+2]
            expr=quicksum(float(A.data[z])*variables[A.indices[z]]for z in range(lo,hi))
            model.addCons(expr==float(M.rhs[i])if M.eq[i]else expr<=float(M.rhs[i]))
        admitted=None
        if initial is not None:
            check_vector(data,initial);sol=model.createOrigSol()
            for v,x in zip(variables,initial):model.setSolVal(sol,v,float(x))
            if not model.checkSol(sol,printreason=False,completely=True,original=True):
                model.freeSol(sol);raise ValueError('SCIP rejected the original-matrix source')
            admitted=bool(model.addSol(sol,free=True))
        construction=time.perf_counter()-began;left=deadline-time.perf_counter()
        if left<=0:raise TimeoutError('native compact matrix construction exhausted allowance')
        model.setParam('limits/time',left);t=time.perf_counter();model.optimize();optimize=time.perf_counter()-t
        status=str(model.getStatus());bound=float(model.getDualbound());vector=None;value=None
        if model.getNSols():
            sol=model.getBestSol();vector=[model.getSolVal(sol,v)for v in variables]
            check_vector(data,vector);value=float(model.getSolObjVal(sol))
        return dict(status=status,bound=bound if math.isfinite(bound)and abs(bound)<1e18 else None,
            value=value,vector=vector,nodes=int(model.getNNodes()),seconds=time.perf_counter()-began,
            construction_seconds=construction,optimize_seconds=optimize,source_stored=admitted)
    finally:
        model.freeProb()
