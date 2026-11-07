"""Frozen six-workflow correctness/root validation; no automatic retries."""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/revision39_compact_prototype_20261106'


def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verified_input(p,expected):
    if sha(p)!=expected:raise ValueError('input drift '+str(p))
    return read(p)


def root_completed(arm,result):
    if not result.get('common',{}).get('complete',False):return False
    if arm['root_only']and arm['method']=='CB':
        return (result.get('status')=='root_validated'and result.get('root_shape')is not None
                and result.get('root_master',{}).get('status')=='optimal')
    return result.get('status')in ('optimal','feasible','root_validated')


def put(p,value):
    from review88_physical.evidence import encode
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('xb')as f:f.write(encode(value))


def design():
    return [dict(index=i,case=c,method=m,root_only=c=='G5',
        internal=90 if c=='G5' else 60,external=120 if c=='G5' else 90)
        for i,(c,m)in enumerate([('H30','CB'),('H30','PQ'),('H42','CB'),('H42','PQ'),('G5','CB'),('G5','PQ')])]


def deadlines(arm,origin,charge,global_deadline):
    external=min(origin+arm['external']-charge,global_deadline-30)
    internal=min(origin+arm['internal']-charge-2,external-2)
    if internal<=origin:raise TimeoutError('source/global budget exhausted before admission')
    return dict(internal=internal,external=external)


def source_fee(arm,previous,root_fee):
    if arm['case']=='G5':
        if not 0<=root_fee<=3:raise ValueError('root source fee exceeds frozen reserve')
        return root_fee
    if arm['method']=='PQ':
        if previous is None or previous['case']!=arm['case']or previous['method']!='CB':
            raise ValueError('same-case compact cost-row source required')
        if not 0<=previous['seconds']<=90:raise ValueError('invalid complete compact-source charge')
        return previous['seconds']
    return 0.


def source_inventory(old):
    sources=dict(old)
    for p in (ROOT/'boundary_context_master').glob('*.py'):sources[p.relative_to(ROOT).as_posix()]=sha(p)
    for p in (ROOT/'tests').glob('test_boundary_*.py'):sources[p.relative_to(ROOT).as_posix()]=sha(p)
    for name in ['docs/2026-11-06-revision39-compact-prototype-protocol.md',
                 'docs/data/revision39_20261106/audit_completed_structural.py']:
        sources[name]=sha(ROOT/name)
    return sources


def verify():
    m=read(OUT/'freeze.json')
    if sha(OUT/'freeze.json')!=(OUT/'freeze.sha256').read_text():raise ValueError('manifest drift')
    for n,h in m['sources'].items():
        if sha(ROOT/n)!=h or sha(OUT/'source_snapshot'/n)!=h:raise ValueError('source drift '+n)
    for row in m['cases'].values():
        if sha(OUT/row['path'])!=row['sha256']:raise ValueError('input drift')
    return m


def prepare():
    if OUT.exists():raise FileExistsError('No repeat or overwrite of this prototype stage')
    origin=time.perf_counter();OUT.mkdir();old=read(ROOT/'results/revision39_structural_20261106/freeze.json')
    sources=source_inventory(old['sources'])
    for n,h in sources.items():
        if sha(ROOT/n)!=h:raise ValueError('old frozen source changed')
        q=OUT/'source_snapshot'/n;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes((ROOT/n).read_bytes())
    functional=ROOT/'results/revision39_primal_functional_20261106';fm=read(functional/'freeze.json');cases={}
    for short,name in [('H30','primal-H30'),('H42','primal-H42')]:
        row=fm['cases'][name];original=functional/row['path']
        payload=verified_input(original,row['sha256']);payload['source_charge']=0
        q=OUT/'inputs'/f'{short}.json';put(q,payload);cases[short]=dict(path=q.relative_to(OUT).as_posix(),sha256=sha(q))
    structural=ROOT/'results/revision39_structural_20261106';row=old['inputs']['Geo-groups5']
    prefix=read(structural/'prefixes.json')['Geo-groups5']
    from revision39.prefix import load_prefix
    public=load_prefix(structural/prefix['path'],prefix['sha256'],row['sha256'])
    if len(public['fleets'])!=1:raise ValueError('one common root source required')
    payload=dict(instance=verified_input(structural/row['path'],row['sha256']),cuts=row['metadata']['boundaries'],
        fleet=public['fleets'][0],source_charge=public['charge_seconds'],
        original_input_sha256=row['sha256'],original_prefix_sha256=prefix['sha256'])
    q=OUT/'inputs/G5.json';put(q,payload);cases['G5']=dict(path=q.relative_to(OUT).as_posix(),sha256=sha(q))
    m=dict(schema='compact-boundary-prototype-v1',preparation_origin=origin,sources=sources,cases=cases,
        arms=design(),global_seconds=720,source_fee_reserve=186,wall_seconds=534,
        memory_bytes=2*1024**3,threads=1,seed=0,retries=0,performance_claim=False,
        proof_expectations=dict(H30=18,H42=40),root_new_variables=1204)
    put(OUT/'freeze.json',m);(OUT/'freeze.sha256').write_text(sha(OUT/'freeze.json'))
    print(json.dumps(dict(status='FROZEN_NO_OPTIMIZATION',sources=len(sources),arms=6)))


def worker(index,origin,charge,global_deadline,probe=False):
    m=verify();arm=m['arms'][index];payload=read(OUT/m['cases'][arm['case']]['path'])
    from review88_physical.evidence import instance
    from boundary_context_master.common import describe,build_base,augment,embed_source
    from review88_physical.physical import validate_fleet
    inst=instance(payload['instance']);source=payload['fleet'];cuts=tuple(payload['cuts'])
    checked=validate_fleet(inst,source['group_plans'])
    if not checked['valid']:raise ValueError('invalid frozen source')
    if probe:
        # Matrix-only preflight. The elementary floor is valid without an LP.
        d=describe(inst,cuts,deadline=min(global_deadline-30,origin+28))
        physical=build_base(d,min(global_deadline-30,origin+28))
        a=augment(d,physical,checked['objective'],deadline=min(global_deadline-30,origin+28))
        vector,parts,actions=embed_source(d,a,source)
        put(OUT/'probes'/f'{index}.json',dict(status='PASS',optimizer_calls=0,
            case=arm['case'],method=arm['method'],shape=list(a.A.shape),source_objective=checked['objective']))
        return
    limits=deadlines(arm,origin,charge,global_deadline);records=()
    if arm['method']=='PQ'and not arm['root_only']:
        prior=read(OUT/'arms'/f'{index-1}.json')
        if prior['case']!=arm['case']or not prior['optimal']:raise ValueError('cost-row source is incomplete')
        records=prior['cuts']
    from revision37.carrier_normalization import install
    normalization_notes=install()
    from .workflow import run
    seq=0
    def checkpoint(row):
        nonlocal seq
        put(OUT/'phases'/str(index)/f'{seq:04d}.json',dict(row,charged_seconds=time.perf_counter()-origin+charge));seq+=1
    try:
        result=run(inst,cuts,source,limits['internal'],method=arm['method'],root_only=arm['root_only'],
                   records=records,checkpoint=checkpoint)
    except TimeoutError as exc:
        result=dict(status='unresolved',optimal=False,reason=str(exc),objective=checked['objective'],
            lower_bound=0,fleet=source,method=arm['method'],root_only=arm['root_only'])
    result.update(case=arm['case'],index=index,source_charge=charge,normalization_notes=normalization_notes,
        charged_seconds=time.perf_counter()-origin+charge)
    put(OUT/'arms'/f'{index}.json',result)
    if result['charged_seconds']>arm['internal']:raise TimeoutError('internal workflow overrun')
    if not arm['root_only']and (not result['optimal']or result['objective']!=m['proof_expectations'][arm['case']]):
        raise ValueError('tiny full-tree consistency did not pass; stop remaining workflows')
    print(json.dumps(dict(status=result['status'],case=arm['case'],method=arm['method'],
        lower_bound=str(result['lower_bound']),objective=str(result['objective']))))


def command(index,origin,charge,global_deadline,probe):
    cmd=[sys.executable,'-B','-X','utf8','-m','boundary_context_master.prototype','worker',
        '--index',str(index),'--origin',repr(origin),'--charge',repr(charge),'--deadline',repr(global_deadline)]
    return cmd+(['--probe']if probe else[])


def child(arm,charge,global_deadline,probe=False):
    import psutil
    from revision31.guard import guarded_command
    origin=time.perf_counter();limit=min(origin+30,global_deadline-30)if probe else deadlines(arm,origin,charge,global_deadline)['external']
    if limit<=origin:raise TimeoutError('global allowance cannot admit a worker')
    path=OUT/('probe_logs'if probe else'logs')/f'{arm["index"]}.txt';path.parent.mkdir(exist_ok=True)
    peak=0;samples=0;missing=0;cause=None
    with path.open('xb')as log:
        p=subprocess.Popen(guarded_command(command(arm['index'],origin,charge,global_deadline,probe),limit,os.getpid()),
            cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        while p.poll()is None:
            try:
                proc=psutil.Process(p.pid);rss=sum(q.memory_info().rss for q in [proc]+proc.children(recursive=True))
                peak=max(peak,rss);samples+=1
                if rss>2*1024**3:cause='memory_limit';p.kill()
            except(psutil.NoSuchProcess,psutil.AccessDenied):missing+=1
            if time.perf_counter()>limit+.3:cause='time_limit';p.kill()
            time.sleep(.02)
        code=p.wait()
    elapsed=time.perf_counter()-origin+.1
    receipt=dict(index=arm['index'],case=arm['case'],method=arm['method'],probe=probe,
        seconds=elapsed,source_charge=charge,charged_seconds=elapsed+charge,returncode=code,
        outcome=cause or('complete'if code==0 else'error_or_guard_deadline'),
        sampled_peak_tree_bytes=peak,samples=samples,missing_samples=missing)
    raw=OUT/('probes'if probe else'arms')/f'{arm["index"]}.json'
    receipt['raw_sha256']=sha(raw)if raw.exists()else None
    put(OUT/('probe_receipts'if probe else'receipts')/f'{arm["index"]}.json',receipt)
    if time.perf_counter()-origin>elapsed:raise ValueError('receipt sealing exceeded conservative charge')
    if code!=0 or cause or not raw.exists():raise ValueError(('worker stopped',arm['index'],receipt))
    if not probe and receipt['charged_seconds']>arm['external']:raise TimeoutError('external charged overrun')
    return receipt


def analyze(m,receipts):
    import numpy as np
    from scipy.sparse import coo_matrix
    from review88_physical.evidence import instance,encode
    from review88_physical.physical import validate_fleet
    helper=ROOT/'docs/data/revision39_20261106/audit_completed_structural.py'
    spec=importlib.util.spec_from_file_location('compact_saved_price_audit',helper)
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    rows=[];maxerr=0.;price_count=0;incomplete=[]
    for arm,receipt in zip(m['arms'],receipts):
        r=read(OUT/'arms'/f'{arm["index"]}.json');payload=read(OUT/m['cases'][arm['case']]['path'])
        physical=validate_fleet(instance(payload['instance']),r['fleet']['group_plans'])
        if not physical['valid']or str(physical['objective'])!=str(r['objective']):raise ValueError('independent fleet replay differs')
        common=r.get('common')
        if not root_completed(arm,r):
            incomplete.append(arm['index'])
            rows.append(dict(index=arm['index'],case=arm['case'],method=arm['method'],
                status='UNRESOLVED_COMMON_ROOT',lower=r['lower_bound'],upper=r['objective']))
            continue
        old=common['original_strong_LP']['matrix'];new=common['augmented_strong_LP']['matrix']
        for record in [common['original_strong_LP'],common['augmented_strong_LP']]:
            audit.check_price(record['matrix'],dict(price=record['price'],scale=1,
                exact_dual=record['exact_bound'],bound=record['integer_bound']));price_count+=1
        def sparse(v):return coo_matrix((v['data'],(v['row'],v['col'])),shape=v['shape']).tocsr()
        A=sparse(old);B=sparse(new);nr,nc=A.shape
        if (B[:nr,:nc]-A).nnz or B[:nr,nc:].nnz:raise ValueError('old strong matrix not retained')
        for field in ['cost','rhs','equalities','upper','integer']:
            n=nc if field in ['cost','upper','integer']else nr
            if old[field]!=new[field][:n]:raise ValueError('old strong inventory drift '+field)
        anchor=common['common_integer_anchor']
        if r['lower_bound']+1e-6<anchor or r['lower_bound']>r['objective']+1e-6:raise ValueError('root/global bound inconsistency')
        if not arm['root_only']and (not r['optimal']or r['objective']!=m['proof_expectations'][arm['case']]):
            raise ValueError('tiny optimum mismatch')
        if arm['method']=='PQ'and not arm['root_only']:
            previous=read(OUT/'arms'/f'{arm["index"]-1}.json')
            if encode(previous['cuts'])!=encode(common['shared_cost_rows']):raise ValueError('cost rows not synchronized')
        if arm['method']=='CB'and arm['root_only']:
            if r['root_shape'][1]!=1204 or r['root_master']['status']!='optimal':raise ValueError('G5 compact root incomplete')
        rows.append(dict(index=arm['index'],case=arm['case'],method=arm['method'],
            root_only=arm['root_only'],status=r['status'],lower=r['lower_bound'],upper=r['objective'],
            common_anchor=anchor,old_LP=common['original_strong_LP']['raw_objective'],
            augmented_LP=common['augmented_strong_LP']['raw_objective'],
            charged_seconds=receipt['charged_seconds'],sampled_peak_tree_bytes=receipt['sampled_peak_tree_bytes']))
    return dict(status='STOP_INCOMPLETE_ROOT'if incomplete else'PASS_FUNCTIONAL_AND_ROOT_NOT_PERFORMANCE',rows=rows,
        incomplete_workflows=incomplete,
        stage_prices_independently_recomputed=price_count,physical_fleets_replayed=6,
        original_rows_retained=True,performance_claim=False,new_instances=0,retries=0)


def stage():
    m=verify();origin=m['preparation_origin'];deadline=origin+m['wall_seconds'];receipts=[];fees=0;stop=None
    put(OUT/'stage-start.json',dict(pid=os.getpid(),deadline=deadline))
    try:
        for arm in m['arms']:child(arm,0,deadline,True)
        put(OUT/'preflight-complete.json',dict(status='PASS',probes=6,optimizer_calls=0,seconds=time.perf_counter()-origin))
        for arm in m['arms']:
            payload=read(OUT/m['cases'][arm['case']]['path'])
            prior=receipts[-1]if receipts else None
            fee=source_fee(arm,prior,payload['source_charge']);fees+=fee
            if fees>m['source_fee_reserve']:raise ValueError('source-fee reserve exceeded')
            receipts.append(child(arm,fee,deadline))
        verify();analysis=analyze(m,receipts);put(OUT/'analysis.json',analysis)
        if analysis['status']!='PASS_FUNCTIONAL_AND_ROOT_NOT_PERFORMANCE':
            raise TimeoutError('common root validation incomplete; no retries')
    except BaseException as exc:stop=repr(exc)
    charged=time.perf_counter()-origin+fees
    status='PASS'if stop is None and len(receipts)==6 and charged<=m['global_seconds']else'STOP'
    put(OUT/'complete.json',dict(status=status,stop=stop,receipts=receipts,
        all_stage_charged_seconds=charged,repeated_source_fees=fees,performance_claim=False))
    print(json.dumps(dict(status=status,stop=stop,completed=len(receipts),charged_seconds=charged)))
    if status!='PASS':raise SystemExit(1)


def launch():
    from revision31.guard import guarded_command
    m=verify();deadline=m['preparation_origin']+m['wall_seconds']
    put(OUT/'launch-start.json',dict(pid=os.getpid(),deadline=deadline))
    cmd=[sys.executable,'-B','-X','utf8','-m','boundary_context_master.prototype','stage']
    result=subprocess.run(guarded_command(cmd,deadline,os.getpid()),cwd=ROOT,
        timeout=max(.01,deadline-time.perf_counter())+2,creationflags=subprocess.CREATE_NO_WINDOW)
    put(OUT/'global-guard-complete.json',dict(returncode=result.returncode,
        finished=time.perf_counter(),deadline=deadline))
    raise SystemExit(result.returncode)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','launch','stage','worker'])
    p.add_argument('--index',type=int);p.add_argument('--origin',type=float);p.add_argument('--charge',type=float,default=0)
    p.add_argument('--deadline',type=float);p.add_argument('--probe',action='store_true');a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='launch':launch()
    elif a.action=='stage':stage()
    else:worker(a.index,a.origin,a.charge,a.deadline,a.probe)
