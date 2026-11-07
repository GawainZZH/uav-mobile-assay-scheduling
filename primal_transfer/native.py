"""Checked primal starts on the frozen reusable SCIP adapter's matrix domain.

The run loop is derived from revision22.native.Proxy so frozen source remains
unchanged. The additions are pre-search seed registration and checked admission
after the correct RHS is loaded. No source lower bound or tree is imported.
"""
from types import SimpleNamespace
import hashlib
import math
import time
import highspy
import numpy as np
from pyscipopt import Model,quicksum
from revision22.native import Proxy


class FragmentProxy(Proxy):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self._seeds={};self._submitted=None
        self.seeds_admitted=0;self.seeds_stored=0
        self.seed_seconds=0.;self.seed_registration_seconds=0.;self.seed_events=[]

    @staticmethod
    def _key(lower,rhs):
        return (np.asarray(lower,dtype=np.float64).tobytes(),
                np.asarray(rhs,dtype=np.float64).tobytes())

    @property
    def seed_count(self):return len(self._seeds)

    def register_seed(self,lower,rhs,vector,objective):
        """Register before the first run, after full original-matrix validation."""
        began=time.perf_counter()
        try:
            if self.builds:raise ValueError('register all source starts before native search')
            low=np.asarray(lower,dtype=float);right=np.asarray(rhs,dtype=float)
            x=np.asarray(vector,dtype=float);objective=float(objective)
            if x.shape!=self.c.shape or low.shape!=self.rhs.shape or right.shape!=self.rhs.shape:
                raise ValueError('seed dimensions do not match the fixed matrix')
            if not np.isfinite(x).all() or not np.isfinite(right).all() or np.isnan(low).any() or not math.isfinite(objective):
                raise ValueError('nonfinite primal seed')
            if np.any(low>right):raise ValueError('inconsistent row bounds')
            a=self.A@x;integer=np.asarray(self.integer,dtype=bool)
            error=max(np.max(low-a,initial=0),np.max(a-right,initial=0),
                      np.max(-x,initial=0),np.max(x-self.upper,initial=0),
                      np.max(abs(x[integer]-np.rint(x[integer])),initial=0),
                      abs(float(self.c@x)-objective))
            if error>1e-6:raise ValueError(('registered seed residual',error))
            key=self._key(low,right);old=self._seeds.get(key)
            if old is None or objective<old['objective']:
                saved=x.copy();saved.setflags(write=False)
                self._seeds[key]=dict(vector=saved,objective=objective)
        finally:self.seed_registration_seconds+=time.perf_counter()-began

    def _admit(self,key):
        seed=self._seeds.get(key)
        if seed is None or self._submitted==key:return
        began=time.perf_counter();model=self.model
        try:
            sol=model.createOrigSol()
            for var,value in zip(self.variables,seed['vector']):model.setSolVal(sol,var,float(value))
            if not model.checkSol(sol,printreason=False,completely=True,original=True):
                model.freeSol(sol)
                raise ValueError('native rejected a matrix-validated context seed')
            stored=bool(model.addSol(sol,free=True))
            self._submitted=key;self.seeds_admitted+=1;self.seeds_stored+=int(stored)
            self.seed_events.append(dict(context_hash=hashlib.sha256(key[0]+key[1]).hexdigest(),
                objective=seed['objective'],checked_original=True,stored=stored))
        finally:self.seed_seconds+=time.perf_counter()-began

    def run(self):
        began=time.perf_counter();deadline=began+self.time_limit-self.runtime
        requested=self._key(self.lower,self.rhs)
        if self.model is not None and not self.reuse:self.close()
        if self.model is None:
            model=Model();model.hideOutput();self.model=model;self.builds+=1;self._submitted=None
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
            self.model.freeTransform();self.context_resets+=1;self._submitted=None
            for i,cons in enumerate(self.constraints):
                self.model.chgRhs(cons,float(self.rhs[i]))
                self.model.chgLhs(cons,float(self.lower[i]))
            self.loaded=requested
        else:self.resumes+=1
        if time.perf_counter()<deadline:self._admit(requested)
        # Admission belongs to construction, and its time is deducted before solve.
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
