"""Same rows, append-only columns/cost cuts and retained native simplex basis."""
import time
from itertools import combinations
import highspy
import numpy as np
from scipy.sparse import coo_matrix,hstack,vstack,csc_matrix
from fleet_sample_assay.physical import vehicle_arcs
from fleet_sample_assay.integrated_master import route_rows
from fleet_sample_assay.master import Matrix,alpha
from fleet_sample_assay.guarded import _dump
from fleet_sample_assay import guarded


def highs():
    h=highspy.Highs()
    for key,value in (('output_flag',False),('threads',1),('random_seed',0),
                      ('primal_feasibility_tolerance',1e-8),('dual_feasibility_tolerance',1e-8),('mip_rel_gap',0.0)):
        check(h.setOptionValue(key,value))
    return h


def check(status):
    if status==highspy.HighsStatus.kError:raise ValueError('native HiGHS API error')


def pass_lp(h,A,c,upper,lower,rhs,integer=None):
    A=A.tocsc();lp=highspy.HighsLp();lp.num_col_=A.shape[1];lp.num_row_=A.shape[0]
    lp.col_cost_=c;lp.col_lower_=np.zeros(len(c));lp.col_upper_=upper
    lp.row_lower_=lower;lp.row_upper_=rhs
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_=A.indptr;lp.a_matrix_.index_=A.indices;lp.a_matrix_.value_=A.data
    if integer is not None:lp.integrality_=[highspy.HighsVarType.kInteger if x else highspy.HighsVarType.kContinuous for x in integer]
    check(h.passModel(lp))


def native_run(h,deadline):
    left=deadline-time.perf_counter()
    if left<=0:raise TimeoutError('native model build exhausted total budget')
    # Native run time is cumulative across reoptimizations on this object.
    check(h.setOptionValue('time_limit',h.getRunTime()+left))
    check(h.run());status=h.getModelStatus()
    code=0 if status==highspy.HighsModelStatus.kOptimal else (2 if status==highspy.HighsModelStatus.kInfeasible else 1)
    return code,h.getSolution(),h.getInfo(),h.modelStatusToString(status)


class IncrementalMaster:
    def __init__(self,inst,branches=(),use_sri=True):
        self.inst=inst;self.branches=tuple(branches);self.use_sri=use_sri
        self.arcs=vehicle_arcs(inst);self.routes=[];self.cuts=[];self.phase=1
        self.builds=1;self.solve_calls=0;self.basis_reuses=0;self.h=highs()
        M=Matrix();self.M=M
        for g in range(inst.groups):
            for p in range(len(inst.pads)):
                for t in range(inst.H+1):M.row(('flow',g,p,t),int(p==0 and t==0)-int(p==0 and t==inst.H),True)
            for p in range(len(inst.pads)):
                for t in range(inst.H):M.row(('launch',g,p,t));M.row(('recover',g,p,t))
            for k in range(inst.H):M.row(('drone',g,k),inst.drones);M.row(('battery',g,k),inst.batteries)
            for a in range(inst.H):
                for b in range(a,inst.H):M.row(('hall',g,a,b))
        for j in range(len(inst.jobs)):M.row(('cover',j),1,True)
        for i,(key,value) in enumerate(branches):M.row(('branch',i),value,True)
        if use_sri:
            for S in combinations(range(len(inst.jobs)),3):M.row(('sri',)+S,1)
        self.zids=[];self.rids=[]
        for ai,(g,p,t,q,u,cost) in enumerate(self.arcs):
            col={('flow',g,p,t):1,('flow',g,q,u):-1,('launch',g,p,t):-inst.drones,('recover',g,p,t):-inst.drones}
            if p==q:
                for a in range(t+1):
                    for b in range(t,inst.H):col['hall',g,a,b]=-inst.stations
            for i,(key,value) in enumerate(branches):
                if key==('vehicle',ai):col['branch',i]=1.0
            self.zids.append(M.column(col,upper=1,integer=True))
        self.theta=M.column({},upper=0)
        self.cover_aids=[M.column({('cover',j):1},1) for j in range(len(inst.jobs))]
        self.slack_aids=[]
        for key,eq in zip(M.keys,M.eq):
            if eq:self.slack_aids.append(M.column({key:1},1))
            self.slack_aids.append(M.column({key:-1},1))
        self.aids=self.cover_aids+self.slack_aids
        self.A=M.sparse().tocsc();self._pass()

    def _pass(self):
        M=self.M
        pass_lp(self.h,self.A,np.array(M.cost),np.array(M.upper),np.where(M.eq,M.rhs,-np.inf),np.array(M.rhs))

    def _append_routes(self,pool):
        if self.routes!=list(pool[:len(self.routes)]):raise ValueError('route pool must be append-only')
        M=self.M;inst=self.inst;old=len(M.cost)
        for r in pool[len(self.routes):]:
            self.rids.append(M.column(route_rows(inst,r,self.cuts,False,self.branches,self.use_sri),
                                      r.cost(inst) if self.phase==2 else 0,integer=True))
            self.routes.append(r)
        count=len(M.cost)-old
        if not count:return
        rows=[];cols=[];data=[]
        for j,col in enumerate(M.columns[old:]):
            for key,value in col.items():
                if value:rows.append(M.index[key]);cols.append(j);data.append(value)
        a=coo_matrix((data,(rows,cols)),shape=(len(M.keys),count)).tocsc()
        check(self.h.addCols(count,M.cost[old:],np.zeros(count),M.upper[old:],a.nnz,a.indptr,a.indices,a.data))
        self.A=hstack((self.A,a),format='csc')

    def _append_cuts(self,cuts):
        if self.cuts!=list(cuts[:len(self.cuts)]):raise ValueError('cut pool must be append-only')
        M=self.M;inst=self.inst;old=len(M.keys);rows=[];cols=[];data=[]
        for beta in cuts[len(self.cuts):]:
            n=len(self.cuts);key=('cost',n);M.row(key);self.cuts.append(beta)
            entries={self.theta:-1}
            for ai,a in enumerate(self.arcs):
                if a[1]==a[3]:entries[self.zids[ai]]=inst.stations*beta[a[0]][a[2]]
            for idx,r in zip(self.rids,self.routes):entries[idx]=sum(alpha(inst,beta,r.g,j,s,r.recovery) for j,s in r.visits)
            for idx,value in entries.items():
                if value:rows.append(M.index[key]-old);cols.append(idx);data.append(value);M.columns[idx][key]=value
        count=len(M.keys)-old
        if not count:return
        a=coo_matrix((data,(rows,cols)),shape=(count,len(M.cost))).tocsr()
        check(self.h.addRows(count,np.full(count,-np.inf),np.zeros(count),a.nnz,a.indptr,a.indices,a.data))
        self.A=vstack((self.A,a),format='csc')

    def _phase_two(self):
        if self.phase==2:return
        self.phase=2;M=self.M
        for idx,a in zip(self.zids,self.arcs):M.cost[idx]=a[5]
        for idx,r in zip(self.rids,self.routes):M.cost[idx]=r.cost(self.inst)
        M.cost[self.theta]=1;M.upper[self.theta]=np.inf
        for idx in self.aids:M.cost[idx]=0;M.upper[idx]=0
        ids=np.arange(len(M.cost),dtype=np.int32)
        check(self.h.changeColsCost(len(ids),ids,np.array(M.cost)))
        check(self.h.changeColsBounds(len(ids),ids,np.zeros(len(ids)),np.array(M.upper)))

    def solve(self,pool,cuts,phase,deadline,integer=False):
        began=time.perf_counter();M=self.M
        if phase==1 and self.phase!=1:raise ValueError('phase reversal requires a fresh node')
        self._append_routes(pool)
        if phase==2:self._phase_two();self._append_cuts(cuts)
        self.solve_calls+=1
        if self.solve_calls>1 and not integer and self.h.getBasis().valid:self.basis_reuses+=1
        h=self.h
        if integer:
            if phase!=2:raise ValueError('integer phase one not supported')
            h=highs();pass_lp(h,self.A,np.array(M.cost),np.array(M.upper),np.where(M.eq,M.rhs,-np.inf),np.array(M.rhs),M.integer)
        solve_start=time.perf_counter();status,sol,info,message=native_run(h,deadline)
        r=dict(status=status,message=message,objective=None,seconds=time.perf_counter()-began,
               solve_seconds=time.perf_counter()-solve_start,rows=self.A.shape[0],columns=self.A.shape[1],nonzeros=self.A.nnz,
               _matrix=self.A,_rhs=np.array(M.rhs),_eq=np.array(M.eq,dtype=bool),_cost=np.array(M.cost),_upper=np.array(M.upper),
               _keys=tuple(M.keys),_rids=list(self.rids),_zids=list(self.zids),_theta=self.theta,_aids=list(self.aids),_x=None)
        if not sol.value_valid or (not integer and status!=0):return r
        x=np.array(sol.col_value);cost=np.array(M.cost);eq=r['_eq'];rhs=r['_rhs'];ub=r['_upper']
        r.update(objective=float(cost@x),_x=x,route_values=x[self.rids],z_values=x[self.zids],theta=float(x[self.theta]),
                 artificial=x[self.aids],arcs=self.arcs)
        activity=self.A@x
        err=max(np.max(np.abs(activity[eq]-rhs[eq]),initial=0),np.max(np.maximum(0,activity[~eq]-rhs[~eq]),initial=0),
                np.max(np.maximum(0,-x),initial=0),np.max(np.maximum(0,x-ub),initial=0),abs(cost@x-info.objective_function_value))
        if not np.all(np.isfinite(x)) or err>1e-7:raise ValueError(('incremental primal residual',err))
        q=np.zeros((self.inst.groups,self.inst.H))
        for value,a in zip(x[self.zids],self.arcs):
            if a[1]==a[3]:q[a[0],a[2]]+=self.inst.stations*value
        r['q']=q
        if not integer:
            if not sol.dual_valid:raise ValueError('missing LP dual')
            pi=np.array(sol.row_dual)
            if not all(np.all(np.isfinite(v)) for v in (pi,cost,rhs,self.A.data)):
                raise ValueError('nonfinite finite-box certificate data')
            pi[~eq]=np.minimum(0,pi[~eq]);rc=cost-self.A.T@pi
            box=ub.copy();C0=2*self.inst.groups+len(self.inst.jobs)+sum(v for k,v in self.branches)
            box[~np.isfinite(box)]=C0 if phase==1 else 0
            box[self.rids]=1+C0 if phase==1 else 1
            box[self.theta]=0 if phase==1 else self.inst.H*sum(j.weight for j in self.inst.jobs)
            if not np.all(np.isfinite(box)) or np.any(box<0):raise ValueError('invalid finite box')
            if np.max(x-box,initial=0)>1e-7:raise ValueError('point violates finite box')
            lower=float(rhs@pi+np.minimum(0,rc)@box)
            # Weak duality needs a valid finite box and sign-safe pi, not a tight pi.
            # A positive slack is conservatism; only a crossed bound is contradictory.
            slack=r['objective']-lower
            if not np.isfinite(lower) or not np.all(np.isfinite(rc)) or slack < -1e-7:
                raise ValueError(('invalid finite-box lower bound',slack))
            r.update(dual={k:float(v) for k,v in zip(M.keys,pi)},verified_dual_bound=lower,
                     dual_certificate_slack=slack,route_box_upper=1+C0 if phase==1 else 1)
        r['seconds_with_validation']=time.perf_counter()-began
        if guarded.RAW_DIR is not None:
            coo=self.A.tocoo()
            _dump('global_node',dict(instance=self.inst,routes=pool,cuts=cuts,branches=self.branches,phase=phase,integer=integer,
                                 keys=M.keys,cost=M.cost,rhs=M.rhs,eq=M.eq,upper=M.upper,
                                 matrix_triplets=(coo.row,coo.col,coo.data,self.A.shape),
                                 point=x,dual=r.get('dual'),verified_bound=r.get('verified_dual_bound')))
        return r


def solve_node(inst,pool,*,cuts=(),branches=(),phase=2,integer=False,deadline=None,use_sri=True):
    return IncrementalMaster(inst,branches,use_sri).solve(pool,cuts,phase,time.perf_counter()+60 if deadline is None else deadline,integer)
