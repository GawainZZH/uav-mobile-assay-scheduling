"""Retain original SCIP matrix; preserve the native tree only for unchanged RHS."""
from types import SimpleNamespace
import math
import time
import highspy
import numpy as np
from pyscipopt import Model,quicksum
from service_separator25_joint.scip_backend import Proxy as ColdProxy


class Proxy(ColdProxy):
    def __init__(self,*args,reuse=True,**kwargs):
        super().__init__(*args,**kwargs);self.reuse=reuse;self.model=None;self.variables=[];self.constraints=[]
        self.loaded=None;self.builds=0;self.context_resets=0;self.resumes=0;self.construction_seconds=0.;self.optimize_seconds=0.

    def close(self):
        if self.model is not None:self.model.freeProb();self.model=None

    def run(self):
        began=time.perf_counter();deadline=began+self.time_limit-self.runtime
        requested=(self.lower.tobytes(),self.rhs.tobytes())
        if self.model is not None and not self.reuse:self.close()
        if self.model is None:
            model=Model();model.hideOutput();self.model=model;self.builds+=1
            for k,v in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
                        ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.)]:model.setParam(k,v)
            self.variables=[model.addVar(name='x'+str(j),vtype='B'if self.integer[j]and self.upper[j]<=1 else'I'if self.integer[j]else'C',
                lb=0,ub=float(self.upper[j]),obj=float(self.c[j]))for j in range(len(self.c))]
            self.constraints=[]
            for i in range(self.A.shape[0]):
                a,b=self.A.indptr[i:i+2]
                expr=quicksum(float(self.A.data[k])*self.variables[self.A.indices[k]]for k in range(a,b))
                upper=model.addCons(expr<=float(self.rhs[i]),name='r'+str(i))
                if math.isfinite(self.lower[i]):model.chgLhs(upper,float(self.lower[i]))
                self.constraints.append(upper)
            self.loaded=requested
        elif self.loaded!=requested:
            self.model.freeTransform();self.context_resets+=1
            for i,cons in enumerate(self.constraints):
                self.model.chgRhs(cons,float(self.rhs[i]))
                self.model.chgLhs(cons,float(self.lower[i]))
            self.loaded=requested
        else:self.resumes+=1
        self.construction_seconds+=time.perf_counter()-began
        model=self.model;left=deadline-time.perf_counter()
        self.solution=SimpleNamespace(value_valid=False,dual_valid=False,col_value=[])
        if left>0:
            t=time.perf_counter();model.setParam('limits/time',model.getSolvingTime()+left);model.optimize()
            self.optimize_seconds+=time.perf_counter()-t;status=str(model.getStatus())
            if model.getNSols():
                sol=model.getBestSol();x=[model.getSolVal(sol,v)for v in self.variables]
                self.solution=SimpleNamespace(value_valid=True,dual_valid=False,col_value=x)
                value=float(model.getSolObjVal(sol))
            else:value=math.inf
            bound=float(model.getDualbound());nodes=int(model.getNNodes())
        else:status='local_construction_limit';value=math.inf;bound=-math.inf;nodes=0
        self.message='SCIP '+status
        self.code=highspy.HighsModelStatus.kOptimal if status=='optimal'else highspy.HighsModelStatus.kInfeasible if status=='infeasible'else highspy.HighsModelStatus.kTimeLimit
        self.info=SimpleNamespace(mip_dual_bound=bound,mip_node_count=nodes,objective_function_value=value)
        self.runtime+=time.perf_counter()-began
        return highspy.HighsStatus.kOk
