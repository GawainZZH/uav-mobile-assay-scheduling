from copy import deepcopy
from fractions import Fraction
from math import ceil,lcm
import time
import numpy as np
from service_separator25_joint.model import build
from service_separator25_joint.dual import price
from fleet_sample_assay_certified.master import highs,pass_lp,native_run
from review88_physical.milp import dump_matrix

def relaxed_boundary_model(inst,source=None,deadline=None):
    m=deepcopy(source)if source is not None else build(inst,deadline=deadline)
    M=m.M;m.base_columns=len(M.cost);m.free_start={};m.free_end={}
    for g in range(inst.groups):
        M.eq[M.index['pack',g]]=False;M.rhs[M.index['pack',g]]=len(inst.jobs)
        M.row(('free_start',g),1,True);M.row(('free_end',g),1,True)
        for p in range(len(inst.pads)):
            M.rhs[M.index['flow',g,p,0]]=0;M.rhs[M.index['flow',g,p,inst.H]]=0
            m.free_start[g,p]=M.column({('flow',g,p,0):-1,('free_start',g):1},0,upper=1)
            m.free_end[g,p]=M.column({('flow',g,p,inst.H):1,('free_end',g):1},0,upper=1)
    m.A=M.sparse().tocsc();return m

def embed_context(m,vector,context):
    start,end,consumed=context;x=np.zeros(len(m.M.cost));x[:m.base_columns]=vector
    for g in range(m.inst.groups):x[m.free_start[g,start[g]]]=1;x[m.free_end[g,end[g]]]=1
    return x

def envelope(inst,deadline,source=None):
    began=time.perf_counter();m=relaxed_boundary_model(inst,source,deadline)
    scale=lcm(*(Fraction(c).denominator for c in m.M.cost))
    c=np.array([int(Fraction(x)*scale)for x in m.M.cost],dtype=float)
    u=np.array(m.M.upper,dtype=float);rhs=np.array(m.M.rhs);eq=np.array(m.M.eq,dtype=bool)
    h=highs();pass_lp(h,m.A,c,u,np.where(eq,rhs,-np.inf),rhs)
    status,sol,info,msg=native_run(h,deadline)
    record=dict(status='unresolved',bound=0,seconds=time.perf_counter()-began,scale=scale)
    if status==0 and sol.dual_valid:
        p=price(m.A,c,u,eq,sol.row_dual);exact=p.evaluate(rhs)
        record.update(status='optimal_lp',bound=Fraction(ceil(exact),scale),
            exact_dual=exact,matrix=dump_matrix(m,c,sol.col_value),
            price=dict(weights=p.weights.tolist(),constant=p.constant,scale=p.scale),
            raw_lp_objective=float(c@np.asarray(sol.col_value)))
    elif status==2:record.update(status='infeasible_lp')
    record['seconds']=time.perf_counter()-began
    return record
