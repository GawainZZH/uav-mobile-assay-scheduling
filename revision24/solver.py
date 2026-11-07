"""Instrumented common MILP and unchanged SK; exact QC-cost face selection."""
from contextlib import contextmanager
from fractions import Fraction
from math import lcm
import importlib
import math
import time
import numpy as np
from pyscipopt import Model,quicksum,Eventhdlr,SCIP_EVENTTYPE
from review88_physical import milp as base
from service_separator25_joint.model import complete
from service_separator25_joint.monolithic import strengthening
from service_separator25_joint.scip_backend import Proxy
from revision22.clocks import pipeline
from review88_physical import physical
from revision22.strengthening import boundary_rule
from .candidate import admit,validate_candidate
from .quality_calendar import apply_model, release_epoch, validate_quality


def service_cost(m):
    scale=lcm(*(j.weight.denominator for j in m.inst.jobs))
    weights=[int(j.weight*scale)for j in m.inst.jobs]
    big=1+sum(w*j.deadline**2 for w,j in zip(weights,m.inst.jobs))
    c=np.zeros(len(m.M.cost))
    for (g,j,s,e,k),col in m.assignments.items():
        C=k+m.inst.assay_duration;c[col]=weights[j]*(big*C+C*C)
    return c,big,scale


def add_cost_face(m,face):
    key=('QC_combined_cost_optimal_face',);m.M.row(key,float(face),True)
    for col,cost in enumerate(m.M.cost):
        if cost:m.M.columns[col][key]=float(cost)


def checked_vector(m,raw):
    raw=np.asarray(raw);rhs=np.asarray(m.M.rhs);eq=np.asarray(m.M.eq,dtype=bool)
    x=complete(m,raw,(0,)*m.inst.groups,(0,)*m.inst.groups)
    a=m.A@x
    error=max(np.max(abs(a[eq]-rhs[eq]),initial=0),np.max(a[~eq]-rhs[~eq],initial=0),
        np.max(-x,initial=0),np.max(x-np.asarray(m.M.upper),initial=0),
        np.max(abs(x-np.rint(x)),initial=0))
    if error>1e-6:raise ValueError(('matrix',error))
    fleet,p=base.restore(m,x)
    return x,fleet,p


def transport_cost(m,quality_blocks=(),result_allowance=None):
    """Transportation first, then physical or quality-released completion."""
    inst=m.inst
    scale=lcm(*(Fraction(v).denominator for v in m.transport))
    W=lcm(*(j.weight.denominator for j in inst.jobs))
    maxima=[min(inst.H,j.lo+result_allowance) if quality_blocks else j.deadline for j in inst.jobs]
    cap=1+sum(int(j.weight*W)*maximum for j,maximum in zip(inst.jobs,maxima))
    c=np.asarray(m.transport,dtype=float)*scale*cap
    for (g,j,s,e,k),col in m.assignments.items():
        completion=k+inst.assay_duration
        release=release_epoch(completion,quality_blocks) if quality_blocks else completion
        c[col]=int(inst.jobs[j].weight*W)*(release if release is not None else 0)
    return c,dict(transport_scale=scale,service_scale=W,service_range=cap-1,scalar_priority=cap)


def _native_solve(inst,boundaries,seconds=60,policy='J',face=None,first_round=None,*,floors=False,responsibility=None,transport=False,qc_busy=(),progress=None,joint=True,initial_fleets=(),admission_sink=None,quality_blocks=(),result_allowance=None):
    began=time.perf_counter();deadline=began+seconds;reserve=min(1.,.1*seconds);metadata={};events=[];trace=[];best_seen=math.inf;last_trace=-1;last_bound=0
    r=dict(status='unresolved',optimal=False,objective=None,lower_bound=None,fleet=None,stages=[])
    with boundary_rule(joint),strengthening(inst,boundaries,metadata):
        m=base.build(inst,policy,deadline-reserve)
        if len(inst.jobs)==inst.groups*inst.batteries:
            for g in range(inst.groups):m.M.eq[m.M.index['pack',g]]=True
        if floors:
            from service_separator25_joint.cases import split
            from service_separator25_joint.solve import signature
            from referee_evidence21.envelope import envelope
            cert=metadata['certification'];cuts=(0,)+tuple(c['epoch']for c in cert['cuts']if c['valid'])+(inst.H,)
            env={};offset=0
            for stage,(local,ids)in enumerate(split(inst,cuts)):
                key=signature(local)
                if key not in env:env[key]=envelope(local,deadline-reserve)
                bound=env[key]['bound'];row=('shared_stage_cost_floor',stage);m.M.row(row,-float(bound))
                selected=set(ids)
                for (g,i),col in m.rids.items():
                    if m.routes[i].job in selected:m.M.columns[col][row]=-float(m.M.cost[col])
                for (g,j,s,e,k),col in m.assignments.items():
                    if j in selected:m.M.columns[col][row]=-float(m.M.cost[col])
                for (g,i),col in m.zids.items():
                    p,t,b,u=m.arcs[i]
                    if offset<=t<u<=offset+local.H and m.M.cost[col]:m.M.columns[col][row]=-float(m.M.cost[col])
                offset+=local.H
            r['envelopes']=env
        if responsibility is not None:
            for (g,i),col in m.rids.items():
                if g!=responsibility[m.routes[i].job]:m.M.upper[col]=0
        if qc_busy:
            for (g,i),col in m.zids.items():
                p,t,b,u=m.arcs[i]
                if p!=b and any(t<v+1 and v<u for a,z in qc_busy for v in range(a,z)):m.M.upper[col]=0
            for g in range(inst.groups):
                for a,z in qc_busy:
                    for t in range(a,z):m.M.rhs[m.M.index['slot',g,t]]=-1

        if first_round is not None:
            ids=[j for j,x in enumerate(inst.jobs)if x.lo<60]
            key=('fixed_myopic_first_responsibility',);m.M.row(key,len(ids),True)
            for i,route in enumerate(m.routes):
                if route.job in ids:m.M.columns[m.rids[0,i]][key]=1
        if face is not None:
            add_cost_face(m,face)
        if quality_blocks:
            r['quality_model']=apply_model(m,quality_blocks,result_allowance)
            r.update(quality_blocks=quality_blocks,result_allowance=result_allowance,
                     completion_semantics='successful_quality_release',
                     physical_combined_objective_semantics='transport_plus_weighted_assay_age')
        m.A=m.M.sparse().tocsc()
        c=np.asarray(m.M.cost,dtype=float)
        if face is not None:c,big,scale=service_cost(m);r.update(service_priority=big,weight_scale=scale)
        if transport:
            c,transport_metadata=transport_cost(m,quality_blocks,result_allowance)
            r.update(transport_metadata)

        def checked(raw):
            x,fleet,p=checked_vector(m,raw)
            if quality_blocks:
                quality=validate_quality(inst,fleet,quality_blocks,result_allowance)
                if not quality['valid']:raise ValueError(('quality replay',quality['errors']))
                fleet['quality']=quality
                p=dict(p,quality=quality)
            return x,fleet,p

        model=Model();model.hideOutput()
        for key,val in [('parallel/maxnthreads',1),('lp/threads',1),('randomization/randomseedshift',0),
            ('limits/gap',0.),('limits/absgap',0.),('limits/memory',2048.)]:model.setParam(key,val)
        variables=[model.addVar(name=f'x{j}',vtype='B'if m.M.integer[j]else'C',lb=0,
            ub=float(m.M.upper[j]),obj=float(c[j]))for j in range(len(c))]
        A=m.A.tocsr()
        for i in range(A.shape[0]):
            a,b=A.indptr[i:i+2];expr=quicksum(float(A.data[k])*variables[A.indices[k]]for k in range(a,b))
            model.addCons(expr==float(m.M.rhs[i])if m.M.eq[i]else expr<=float(m.M.rhs[i]))
        class Progress(Eventhdlr):
            def eventinit(self):self.model.catchEvent(SCIP_EVENTTYPE.BESTSOLFOUND|SCIP_EVENTTYPE.NODESOLVED,self)
            def eventexit(self):self.model.dropEvent(SCIP_EVENTTYPE.BESTSOLFOUND|SCIP_EVENTTYPE.NODESOLVED,self)
            def eventexec(self,event):
                nonlocal best_seen,last_trace,last_bound
                sol=self.model.getBestSol();now=time.perf_counter()-began;changed=False
                if sol is not None:
                    native=float(self.model.getSolObjVal(sol))
                    if native<best_seen-1e-6:
                        raw=[self.model.getSolVal(sol,v)for v in variables]
                        try:
                            (x,fleet,p),clocks=pipeline(physical,lambda:checked(raw),began)
                            clocks['detected']=now
                            value=float(p['road_cost']+p['flight_cost']if transport else p['objective'])
                            event_bound=float(self.model.getDualbound())
                            if not math.isfinite(event_bound)or abs(event_bound)>1e18:event_bound=0
                            if transport:event_bound=max(0,math.ceil((event_bound-r['service_range'])/r['scalar_priority']-1e-7)/r['transport_scale'])
                            events.append(dict(**clocks,objective=value,fleet=fleet,lower_bound=max(last_bound,event_bound)))
                            if quality_blocks:events[-1]['weighted_release']=p['quality']['weighted_release']
                            best_seen=native;changed=True
                        except ValueError:return
                bound=float(self.model.getDualbound())
                if not math.isfinite(bound)or abs(bound)>1e18:bound=0
                if transport:bound=max(0,math.ceil((bound-r['service_range'])/r['scalar_priority']-1e-7)/r['transport_scale'])
                last_bound=max(last_bound,bound)
                if changed or now-last_trace>=.2:
                    item=dict(seconds=time.perf_counter()-began,lower_bound=last_bound,objective=events[-1]['objective']if events else None,kind='native_progress')
                    trace.append(item);last_trace=now
                    if progress:progress(item)
        handler=Progress();model.includeEventhdlr(handler,'validated_progress','physical incumbents and global bounds')
        r['candidate_admissions']=[]if admission_sink is None else admission_sink
        for fleet in initial_fleets:
            try:admission=admit(model,variables,m,fleet)
            except Exception:
                model.freeProb()
                raise
            admission['seconds']=time.perf_counter()-began
            r['candidate_admissions'].append(admission)
        r['construction_seconds']=time.perf_counter()-began
        r['matrix_shape']=list(m.A.shape);r['matrix_nonzeros']=m.A.nnz
        left=deadline-reserve-time.perf_counter()
        if left<=0:
            model.freeProb()
            raise TimeoutError('build exhausted budget')
        model.setParam('limits/time',left);model.optimize();r['native_optimize_calls']=1;status=str(model.getStatus())
        bound=float(model.getDualbound());r['native_status']=status
        r['native_objective_bound']=bound if math.isfinite(bound)and abs(bound)<1e18 else None
        r['lower_bound']=(max(0,math.ceil((bound-r['service_range'])/r['scalar_priority']-1e-7)/r['transport_scale'])if transport else bound)if math.isfinite(bound)and abs(bound)<1e18 else None
        if model.getNSols():
            sol=model.getBestSol();raw=np.array([model.getSolVal(sol,v)for v in variables])
            (x,fleet,p),clocks=pipeline(physical,lambda:checked(raw),began);val=float(c@x)
            if face is not None and abs(float(p['objective'])-face)>1e-6:raise ValueError('cost-face completion')
            closed=status=='optimal'and abs(val-bound)<1e-5
            primary=float(p['road_cost']+p['flight_cost']if transport else p['objective'])
            r.update(status='optimal'if closed else'feasible',optimal=closed,objective=primary,physical_combined_objective=float(p['objective']),fleet=fleet,
                selection_objective=val if face is not None else None,
                stages=[dict(matrix=base.dump_matrix(m,c,x),native_vector=raw.tolist(),objective=val)])
            if transport and closed:r['lower_bound']=primary
            events.append(dict(**clocks,objective=primary,fleet=fleet,lower_bound=r['lower_bound'],kind='terminal_snapshot'))
            if quality_blocks:
                r.update(weighted_release=p['quality']['weighted_release'],quality=p['quality'])
                events[-1]['weighted_release']=p['quality']['weighted_release']
        elif status=='infeasible':r.update(status='infeasible')
        model.freeProb()
    r.update(seconds=time.perf_counter()-began,first_validated_fleet=events[0]if events else None,
        shared_strengthening=metadata,trace=trace,primal_events=events,objective_kind='transport'if transport else'combined')
    if r['seconds']>seconds:r.update(status='workflow_time_limit',optimal=False)
    return r


def solve(inst,boundaries,seconds=60,policy='J',face=None,first_round=None,*,floors=False,
          responsibility=None,transport=False,qc_busy=(),progress=None,joint=True,initial_fleets=(),
          quality_blocks=(),result_allowance=None):
    """Cold monolithic suffix; caller charges the shared source prefix separately.

    Candidate events and all returned clocks are relative to this worker's solver
    entry. A seed remains deliverable if subsequent matrix construction times out.
    The shared-candidate comparator is deliberately restricted to the common J
    combined-cost model, so no restricted policy can receive an unchecked seed.
    """
    if quality_blocks:
        if (not transport or initial_fleets or face is not None or floors or qc_busy
                or type(result_allowance) is not int or result_allowance<0):
            raise ValueError('quality mode requires transportation objective, integer allowance, and no seeds/face/floors/qc_busy')
    elif result_allowance is not None:
        raise ValueError('result allowance requires quality blocks')
    options=dict(policy=policy,face=face,first_round=first_round,floors=floors,
        responsibility=responsibility,transport=transport,qc_busy=qc_busy,joint=joint,
        quality_blocks=quality_blocks,result_allowance=result_allowance)
    if not initial_fleets:
        return _native_solve(inst,boundaries,seconds,progress=progress,**options)
    if policy!='J'or face is not None or first_round is not None or responsibility is not None or transport or qc_busy:
        raise ValueError('shared candidates require unrestricted J combined-cost model')
    began=time.perf_counter();seeds=[];events=[];admissions=[]
    for fleet in initial_fleets:
        (validated,check),clocks=pipeline(physical,lambda:validate_candidate(inst,fleet),began)
        seeds.append(validated)
        events.append(dict(**clocks,fleet=validated,objective=float(check['objective']),lower_bound=0,
            kind='shared_candidate',native_admitted=False))
        if progress:progress(dict(seconds=clocks['released'],objective=float(check['objective']),
            lower_bound=0,kind='shared_candidate'))
    offset=time.perf_counter()-began
    def forwarded(row):
        if progress:progress(dict(row,seconds=row['seconds']+offset))
    try:
        if offset>=seconds:raise TimeoutError('candidate replay exhausted suffix allowance')
        r=_native_solve(inst,boundaries,seconds-offset,progress=forwarded,initial_fleets=seeds,
            admission_sink=admissions,**options)
    except TimeoutError as exc:
        r=dict(status='unresolved',optimal=False,objective=None,lower_bound=0,fleet=None,
            reason=str(exc),trace=[],primal_events=[],candidate_admissions=admissions,native_optimize_calls=0,
            objective_kind='combined')
    if r.get('status')=='infeasible':raise ValueError('native infeasibility contradicts validated candidate')
    for row in r.get('trace',[]):row['seconds']+=offset
    for row in r.get('primal_events',[]):
        for phase in ['detected','reconstruction_started','reconstructed','validated','released','seconds']:
            row[phase]+=offset
    for row in r.get('candidate_admissions',[]):row['seconds']+=offset
    if 'construction_seconds'in r:r['construction_seconds']+=offset
    r['primal_events']=events+r.get('primal_events',[])
    eligible=[e for e in r['primal_events']if e['released']<=seconds]
    if eligible:
        best=min(eligible,key=lambda e:e['objective'])
        if r.get('lower_bound')is not None and r['lower_bound']>best['objective']+1e-6:
            raise ValueError('native bound exceeds validated candidate')
        if r.get('objective')is None or abs(best['objective']-r['objective'])>1e-6:
            r.update(status='feasible',optimal=False)
        r.update(objective=best['objective'],fleet=best['fleet'],physical_combined_objective=best['objective'])
    else:r.update(objective=None,fleet=None,optimal=False)
    r['first_validated_fleet']=eligible[0]if eligible else None
    r['candidate_validation_seconds']=offset
    r['seconds']=time.perf_counter()-began
    r.setdefault('native_optimize_calls',0)
    if r['seconds']>seconds:r.update(status='workflow_time_limit',optimal=False)
    return r


