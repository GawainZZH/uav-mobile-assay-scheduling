"""SCIP 10 default branch-and-cut under the same sparse matrix and wall budget."""
from types import SimpleNamespace
import math
import time
import highspy
import numpy as np
from pyscipopt import Model,quicksum


class Proxy:
    def __init__(self,A,c,upper,lower,rhs,integer):
        self.A=A.tocsr();self.c=np.asarray(c);self.upper=np.asarray(upper)
        self.lower=np.asarray(lower).copy();self.rhs=np.asarray(rhs).copy();self.integer=np.asarray(integer)
        self.runtime=0.;self.time_limit=0.;self.code=highspy.HighsModelStatus.kNotset;self.message='not run'
        self.solution=SimpleNamespace(value_valid=False,dual_valid=False,col_value=[])
        self.info=SimpleNamespace(mip_dual_bound=-math.inf,mip_node_count=0,objective_function_value=math.inf)

    def getRunTime(self):return self.runtime
    def setOptionValue(self,key,value):
        if key!='time_limit':raise ValueError('unexpected SCIP adapter option '+key)
        self.time_limit=value;return highspy.HighsStatus.kOk

    def changeRowsBounds(self,n,ids,lower,upper):
        self.lower[np.asarray(ids,dtype=int)]=lower;self.rhs[np.asarray(ids,dtype=int)]=upper
        return highspy.HighsStatus.kOk

    def run(self):
        began=time.perf_counter();deadline=began+self.time_limit-self.runtime
        model=Model();model.hideOutput()
        for k,v in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
                    ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.)]:model.setParam(k,v)
        variables=[model.addVar(name='x'+str(j),vtype='B' if self.integer[j] and self.upper[j]<=1 else 'I' if self.integer[j] else 'C',
            lb=0,ub=float(self.upper[j]),obj=float(self.c[j])) for j in range(len(self.c))]
        for i in range(self.A.shape[0]):
            lo,hi=self.A.indptr[i:i+2]
            expr=quicksum(float(self.A.data[k])*variables[self.A.indices[k]] for k in range(lo,hi))
            if self.lower[i]==self.rhs[i]:model.addCons(expr==float(self.rhs[i]),name='r'+str(i))
            else:
                model.addCons(expr<=float(self.rhs[i]),name='r'+str(i))
                if math.isfinite(self.lower[i]):model.addCons(expr>=float(self.lower[i]),name='l'+str(i))
        left=deadline-time.perf_counter()
        if left<=0:raise TimeoutError('SCIP construction exhausted budget')
        model.setParam('limits/time',left);model.optimize();status=str(model.getStatus())
        self.message='SCIP '+status
        self.code=(highspy.HighsModelStatus.kOptimal if status=='optimal' else
                   highspy.HighsModelStatus.kInfeasible if status=='infeasible' else highspy.HighsModelStatus.kTimeLimit)
        value=math.inf
        if model.getNSols():
            sol=model.getBestSol();x=[model.getSolVal(sol,v) for v in variables]
            value=float(model.getSolObjVal(sol));self.solution=SimpleNamespace(value_valid=True,dual_valid=False,col_value=x)
        else:self.solution=SimpleNamespace(value_valid=False,dual_valid=False,col_value=[])
        self.info=SimpleNamespace(mip_dual_bound=float(model.getDualbound()),mip_node_count=int(model.getNNodes()),
            objective_function_value=value)
        model.freeProb();self.runtime+=time.perf_counter()-began
        return highspy.HighsStatus.kOk

    def getModelStatus(self):return self.code
    def getSolution(self):return self.solution
    def getInfo(self):return self.info
    def modelStatusToString(self,status):return self.message
