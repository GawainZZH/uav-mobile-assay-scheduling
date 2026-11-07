"""One immutable original-G5 full-feedback experiment; no automatic retries.

Algorithm code remains the previous frozen version. This driver adds only the
once-only trial, budget accounting, process receipts and read-only analysis.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/revision39_compact_full_20261106'
BASE=ROOT/'results/revision39_compact_prototype_20261106'


def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,value):
    from review88_physical.evidence import encode
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('xb')as f:f.write(encode(value))


def design(prefix):
    if not 0<=prefix<=3:raise ValueError('original public source fee outside reserve')
    return [dict(index=i,case='G5',method='CB'if i==0 else'PQ',synchronized=i==2,
        root_only=False,root_source_fee=prefix,compute_internal=180,compute_external=210,
        internal=390 if i==2 else 180,external=420 if i==2 else 210)for i in range(3)]


def deadlines(arm,origin,charge,global_deadline):
    # Compute allowance is independent of the synchronized source's cost.
    fee=arm['root_source_fee']
    external=min(origin+arm['compute_external']-fee,
                 origin+arm['external']-charge,global_deadline-60)
    internal=min(origin+arm['compute_internal']-fee-2,
                 origin+arm['internal']-charge-2,external-2)
    if internal<=origin:raise TimeoutError('no remaining admitted workflow budget')
    return dict(internal=internal,external=external)


def source_fee(arm,prefix,cb):
    if not 0<=prefix<=3:raise ValueError('invalid original source fee')
    if not arm['synchronized']:return prefix
    if cb is None or cb['index']!=0 or cb['method']!='CB' or not 0<=cb['charged_seconds']<=210:
        raise ValueError('valid complete CB parent charge required')
    return prefix+cb['charged_seconds']


def admissible(result):
    if not (result.get('common')or{}).get('complete',False):return False
    lower=result['lower_bound'];upper=result['objective']
    if lower>upper+1e-6 or lower+1e-6<result['common']['common_integer_anchor']:
        raise ValueError('invalid retained global interval')
    if result.get('optimal')and abs(upper-lower)>1e-6:raise ValueError('false closure')
    return result['status']in ('optimal','feasible','unresolved')


def check_completion_costs(costs):
    if 'vector'not in costs:raise ValueError('physical completion vector missing')
    if 'native_vector'in costs and costs['vector']>costs['native_vector']+1e-6:
        raise ValueError('completion increased local cost')


def verify():
    m=read(OUT/'freeze.json')
    if sha(OUT/'freeze.json')!=(OUT/'freeze.sha256').read_text():raise ValueError('manifest drift')
    if len(m['sources'])<194:raise ValueError('incomplete frozen source inventory')
    for name,digest in m['sources'].items():
        if sha(ROOT/name)!=digest or sha(OUT/'source_snapshot'/name)!=digest:raise ValueError('source drift '+name)
    if sha(OUT/'input.json')!=m['input_sha256']:raise ValueError('input drift')
    return m


def prepare():
    if OUT.exists():raise FileExistsError('one preparation only; no overwrite')
    origin=time.perf_counter()
    from boundary_context_master.prototype import verify as original_verify
    old=original_verify()
    if read(BASE/'complete.json')['status']!='PASS':raise ValueError('predecessor did not pass')
    payload=read(BASE/old['cases']['G5']['path']);sources=dict(old['sources'])
    for name in ['boundary_context_trial.py','tests/test_compact_full_trial.py',
                 'docs/2026-11-06-revision39-compact-full-protocol.md']:
        sources[name]=sha(ROOT/name)
    OUT.mkdir()
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('source drift before freeze')
        q=OUT/'source_snapshot'/name;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes((ROOT/name).read_bytes())
    put(OUT/'input.json',payload)
    manifest=dict(schema='compact-full-trial-v1',preparation_origin=origin,sources=sources,
        input_sha256=sha(OUT/'input.json'),arms=design(payload['source_charge']),
        global_seconds=1200,wall_seconds=978,source_fee_reserve=219,sealing_reserve=3,
        threads=1,seed=0,memory_bytes=2*1024**3,retries=0,new_instances=0,
        predecessor_freeze_sha256=sha(BASE/'freeze.json'),performance_claim=False)
    put(OUT/'freeze.json',manifest);(OUT/'freeze.sha256').write_text(sha(OUT/'freeze.json'))
    print(json.dumps(dict(status='FROZEN_NO_OPTIMIZATION',sources=len(sources),arms=3)))


def command(index,origin,charge,global_deadline,probe):
    cmd=[sys.executable,'-B','-X','utf8','-m','boundary_context_trial','worker',
        '--index',str(index),'--origin',repr(origin),'--charge',repr(charge),'--deadline',repr(global_deadline)]
    return cmd+(['--probe']if probe else[])


def worker(index,origin,charge,global_deadline,probe=False):
    m=verify();arm=m['arms'][index];payload=read(OUT/'input.json')
    from review88_physical.evidence import instance
    from review88_physical.physical import validate_fleet
    inst=instance(payload['instance']);source=payload['fleet'];cuts=tuple(payload['cuts'])
    checked=validate_fleet(inst,source['group_plans'])
    if not checked['valid']:raise ValueError('invalid original public source')
    if probe:
        from boundary_context_master.common import describe,build_base,augment,embed_source
        limit=min(global_deadline-60,origin+28)
        d=describe(inst,cuts,deadline=limit);physical=build_base(d,limit)
        a=augment(d,physical,checked['objective'],deadline=limit);embed_source(d,a,source)
        put(OUT/'probes'/f'{index}.json',dict(status='PASS',optimizer_calls=0,
            case='G5',method=arm['method'],shape=list(a.A.shape),source_objective=checked['objective']))
        return
    limits=deadlines(arm,origin,charge,global_deadline);records=()
    if arm['synchronized']:
        previous=read(OUT/'arms/0.json')
        if not admissible(previous):raise ValueError('CB lacks completed common root')
        records=previous['cuts']
    from revision37.carrier_normalization import install
    notes=install()
    from boundary_context_master.workflow import run
    seq=0
    def checkpoint(row):
        nonlocal seq
        put(OUT/'phases'/str(index)/f'{seq:04d}.json',dict(row,charged_seconds=time.perf_counter()-origin+charge));seq+=1
    result=run(inst,cuts,source,limits['internal'],method=arm['method'],root_only=False,
               records=records,checkpoint=checkpoint)
    result.update(index=index,case='G5',synchronized=arm['synchronized'],source_charge=charge,
        normalization_notes=notes,charged_seconds=time.perf_counter()-origin+charge,
        workflow_elapsed_seconds=time.perf_counter()-origin)
    put(OUT/'arms'/f'{index}.json',result)
    if not admissible(result):raise ValueError('incomplete common root; stop trial')
    finished=time.perf_counter()-origin
    if finished+charge>arm['internal']:raise TimeoutError('internal complete charge exceeded')
    if finished+arm['root_source_fee']>arm['compute_internal']:
        raise TimeoutError('internal actual-compute allowance exceeded')
    print(json.dumps(dict(index=index,status=result['status'],lower=str(result['lower_bound']),upper=str(result['objective']))))


def audit_helper():
    path=ROOT/'docs/data/revision39_20261106/audit_completed_structural.py'
    spec=importlib.util.spec_from_file_location('compact_full_independent_saved_audit',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def analyze(m,receipts):
    from fractions import Fraction as F
    from scipy.sparse import coo_matrix
    from review88_physical.evidence import instance,encode
    from review88_physical.physical import validate_fleet
    from service_separator25_joint.kernel import validate_boundary
    audit=audit_helper();payload=read(OUT/'input.json');inst=instance(payload['instance'])
    rows=[];vectors=0;prices=0;local_fleets=0;maximum=0.
    def vector(matrix,x,**kwargs):
        nonlocal vectors,maximum
        v=audit.check_vector(matrix,x,**kwargs);vectors+=1;maximum=max(maximum,v['residual']);return v
    def sparse(a):return coo_matrix((a['data'],(a['row'],a['col'])),shape=a['shape']).tocsr()
    for arm,receipt in zip(m['arms'],receipts):
        raw=OUT/'arms'/f'{arm["index"]}.json'
        if sha(raw)!=receipt['raw_sha256']:raise ValueError('raw receipt drift')
        r=read(raw)
        if not admissible(r):raise ValueError('unusable root')
        check=validate_fleet(inst,r['fleet']['group_plans'])
        if not check['valid']or F(str(check['objective']))!=F(str(r['objective'])):raise ValueError('final fleet mismatch')
        if r['full_resource_graph_builds']!=0 or r['historical_bounds_imported']:raise ValueError('unexpected algorithm source')
        c=r['common'];old=c['original_strong_LP']['matrix'];new=c['augmented_strong_LP']['matrix']
        A=sparse(old);B=sparse(new);nr,nc=A.shape
        if (B[:nr,:nc]-A).nnz or B[:nr,nc:].nnz:raise ValueError('old strong matrix changed')
        for name in ['cost','rhs','equalities','upper','integer']:
            n=nc if name in ['cost','upper','integer']else nr
            if old[name]!=new[name][:n]:raise ValueError('old inventory changed '+name)
        vector(new,c['source_vector'])
        for lp in [c['original_strong_LP'],c['augmented_strong_LP']]:
            audit.check_price(lp['matrix'],dict(price=lp['price'],scale=1,exact_dual=lp['exact_bound'],bound=lp['integer_bound']))
            vector(lp['matrix'],lp['matrix']['x'],relax=True);prices+=1
        for lp in r.get('envelopes',{}).values():
            if lp['status']=='optimal_lp':
                audit.check_price(lp['matrix'],lp);vector(lp['matrix'],lp['matrix']['x'],relax=True);prices+=1
        for event in r.get('master_events',[]):
            if event['matrix']['x']is not None:vector(event['matrix'],event['matrix']['x'])
        for event in r.get('local_events',[]):
            rec=event['result'];model=r['local_models'][event['kernel']];local=instance(model['instance'])
            start,end,used=event['context'];matrix=model['matrix'];rhs=list(matrix['rhs'])
            P=len(local.pads);H=local.H;stride=P*(3*H+1)+2*H+1
            for g in range(local.groups):
                for p in range(P):
                    for t in [0,H]:
                        row=g*stride+p*(3*H+1)+t
                        expected=int(p==0 and t==0)-int(p==0 and t==H)
                        if rhs[row]!=expected or not matrix['equalities'][row]:raise ValueError('local flow row order differs')
                        rhs[row]=int(p==start[g])if t==0 else-int(p==end[g])
                row=g*stride+P*(3*H+1)+2*H
                if rhs[row]!=local.batteries:raise ValueError('local stock row order differs')
                rhs[row]=used[g]
            costs={}
            for field in ['vector','native_vector']:
                if rec.get(field)is not None:
                    costs[field]=vector(matrix,rec[field],rhs=rhs,completed=field=='vector')['objective']
            if rec.get('fleet')is not None:
                chk=validate_boundary(local,rec['fleet']['group_plans'],start,end);local_fleets+=1
                if not chk['valid']or F(str(chk['objective']))!=F(str(rec['objective'])):raise ValueError('local fleet mismatch')
                if [len(p['routes'])for p in rec['fleet']['group_plans']]!=list(used):raise ValueError('context stock mismatch')
                if abs(costs['vector']-float(chk['objective']))>1e-6:raise ValueError('context objective mismatch')
                check_completion_costs(costs)
        for delivery in r.get('deliveries',[]):
            chk=validate_fleet(inst,delivery['fleet']['group_plans'])
            if not chk['valid']or F(str(chk['objective']))!=F(str(delivery['objective'])):raise ValueError('delivery mismatch')
        if r.get('final_matrix')is not None:
            matrix=r['final_matrix']
            if r['native']['vector']is not None:vector(matrix,r['native']['vector'])
            if r['completed_vector']is not None:vector(matrix,r['completed_vector'],completed=True)
        if arm['synchronized']:
            if encode(read(OUT/'arms/0.json')['cuts'])!=encode(c['shared_cost_rows']):raise ValueError('cuts not synchronized')
        elif c['shared_cost_rows']:raise ValueError('unexpected initial context rows')
        rows.append(dict(index=arm['index'],method=arm['method'],synchronized=arm['synchronized'],
            status=r['status'],reason=r.get('reason'),lower=r['lower_bound'],upper=r['objective'],
            common_anchor=c['common_integer_anchor'],old_LP=c['original_strong_LP']['raw_objective'],
            augmented_LP=c['augmented_strong_LP']['raw_objective'],master_calls=r.get('master_calls',0),
            local_calls=r.get('local_calls',0),cost_rows=len(r['cuts']),charged_seconds=receipt['charged_seconds'],
            workflow_elapsed_seconds=receipt['seconds'],sampled_peak_tree_bytes=receipt['sampled_peak_tree_bytes']))
    cb,pq,sync=rows
    dominates=cb['lower']>=pq['lower']and cb['upper']<=pq['upper']and(cb['lower']>pq['lower']or cb['upper']<pq['upper'])
    return dict(status='PASS_ONCE_ONLY_FULL_FEEDBACK',rows=rows,saved_vectors_checked=vectors,
        saved_prices_recomputed=prices,local_fleets_replayed=local_fleets,max_residual=maximum,
        CB_interval_strictly_dominates_PQ=dominates,
        synchronized_anchor_strictly_higher=sync['common_anchor']>pq['common_anchor'],
        performance_claim=False,statistical_comparison=False,new_instances=0,retries=0)


def stage():
    from boundary_context_master import prototype as runner
    m=verify();origin=m['preparation_origin'];deadline=origin+m['wall_seconds'];receipts=[];fees=0.;stop=None
    # Reuse reviewed guard/receipt implementation with this driver's command and
    # compute-capped deadlines; no mutation to its source or algorithm code.
    runner.OUT=OUT;runner.command=command;runner.deadlines=deadlines
    put(OUT/'stage-start.json',dict(pid=os.getpid(),deadline=deadline))
    try:
        for arm in m['arms']:runner.child(arm,0,deadline,True)
        put(OUT/'preflight-complete.json',dict(status='PASS',probes=3,optimizer_calls=0))
        prefix=read(OUT/'input.json')['source_charge']
        for arm in m['arms']:
            fee=source_fee(arm,prefix,receipts[0]if receipts else None);fees+=fee
            if fees>m['source_fee_reserve']:raise ValueError('source reserve exceeded')
            receipt=runner.child(arm,fee,deadline)
            if receipt['seconds']+prefix>arm['compute_external']:raise TimeoutError('actual compute workflow exceeded')
            receipts.append(receipt)
        verify();analysis=analyze(m,receipts);put(OUT/'analysis.json',analysis)
    except BaseException as exc:stop=repr(exc)
    charged=time.perf_counter()-origin+fees
    status='PASS'if stop is None and len(receipts)==3 and charged<=m['global_seconds']else'STOP'
    put(OUT/'complete.json',dict(status=status,stop=stop,receipts=receipts,repeated_source_fees=fees,
        all_stage_charged_seconds=charged,performance_claim=False))
    print(json.dumps(dict(status=status,stop=stop,completed=len(receipts),charged_seconds=charged)))
    if status!='PASS':raise SystemExit(1)


def launch():
    from revision31.guard import guarded_command
    m=verify();deadline=m['preparation_origin']+m['wall_seconds']
    put(OUT/'launch-start.json',dict(pid=os.getpid(),deadline=deadline))
    cmd=[sys.executable,'-B','-X','utf8','-m','boundary_context_trial','stage']
    r=subprocess.run(guarded_command(cmd,deadline,os.getpid()),cwd=ROOT,
        timeout=max(.01,deadline-time.perf_counter())+2,creationflags=subprocess.CREATE_NO_WINDOW)
    put(OUT/'global-guard-complete.json',dict(returncode=r.returncode,finished=time.perf_counter(),deadline=deadline))
    raise SystemExit(r.returncode)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','launch','stage','worker'])
    p.add_argument('--index',type=int);p.add_argument('--origin',type=float);p.add_argument('--charge',type=float,default=0)
    p.add_argument('--deadline',type=float);p.add_argument('--probe',action='store_true');a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='launch':launch()
    elif a.action=='stage':stage()
    else:worker(a.index,a.origin,a.charge,a.deadline,a.probe)
