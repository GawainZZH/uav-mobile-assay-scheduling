"""One frozen original-context witness; no native optimization or retries."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/revision39_singleton_witness_20261106'
OLD=ROOT/'results/revision39_compact_full_20261106'


def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,value):
    from review88_physical.evidence import encode
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('xb')as f:f.write(encode(value))


def deadlines(manifest,origin):
    final=manifest['preparation_origin']+manifest['stage_seconds']-manifest['source_charge']-3
    worker=min(origin+manifest['compute_seconds'],final-120)
    if worker<=origin+2:raise TimeoutError('stage allowance cannot admit a new witness')
    return dict(worker=worker,final=final,audit_reserve=120.)


def verify():
    m=read(OUT/'freeze.json')
    if sha(OUT/'freeze.json')!=(OUT/'freeze.sha256').read_text():raise ValueError('manifest drift')
    if len(m['sources'])<197:raise ValueError('incomplete source inventory')
    for name,digest in m['sources'].items():
        if sha(ROOT/name)!=digest or sha(OUT/'source_snapshot'/name)!=digest:raise ValueError('frozen source drift')
    if sha(OUT/'input.json')!=m['input_sha256']:raise ValueError('input drift')
    return m


def prepare():
    if OUT.exists():raise FileExistsError('one numerical witness only')
    origin=time.perf_counter()
    from boundary_context_trial import verify as old_verify
    prior=old_verify();sources=dict(prior['sources']);raw=read(OLD/'arms/0.json');receipt=read(OLD/'receipts/0.json')
    if receipt['returncode']!=0 or sha(OLD/'arms/0.json')!=receipt['raw_sha256']:
        raise ValueError('original context source invalid')
    if len(raw['local_events'])!=9:raise ValueError('fixed source query history differs')
    event=raw['local_events'][8];ctx=event['context']
    if ctx!=[[0]*5,[0]*5,[0,1,1,1,1]]:raise ValueError('unexpected fixed context')
    payload=dict(instance=raw['local_models'][event['kernel']]['instance'],context=ctx,
        source_raw_sha256=receipt['raw_sha256'],source_query_index=8,
        source_query_seconds=event['seconds'],posthoc_original_effective_lower=raw['envelopes'][event['kernel']]['bound'])
    for name in ['singleton_context_oracle.py','singleton_witness.py','tests/test_singleton_context_oracle.py',
        'tests/test_singleton_witness_budget.py',
        'docs/2026-11-06-revision39-singleton-witness-protocol.md',
        'docs/2026-11-06-revision39-singleton-implementation.md',
        'docs/2026-11-06-revision39-singleton-context-proof.md',
        'docs/2026-11-06-revision39-singleton-context-review.md']:
        sources[name]=sha(ROOT/name)
    OUT.mkdir()
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('source changed before witness freeze')
        p=OUT/'source_snapshot'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((ROOT/name).read_bytes())
    put(OUT/'input.json',payload)
    m=dict(schema='singleton-context-witness-v1',preparation_origin=origin,sources=sources,
        input_sha256=sha(OUT/'input.json'),source_charge=receipt['charged_seconds'],
        original_source_receipt_sha256=sha(OLD/'receipts/0.json'),compute_seconds=60,stage_seconds=600,
        memory_bytes=2*1024**3,retries=0,contexts=1,native_optimizer_calls=0,performance_claim=False)
    put(OUT/'freeze.json',m);(OUT/'freeze.sha256').write_text(sha(OUT/'freeze.json'))
    print(json.dumps(dict(status='FROZEN_NO_CONTEXT_QUERY',sources=len(sources),contexts=1)))


def worker(origin,deadline):
    verify();payload=read(OUT/'input.json')
    from review88_physical.evidence import instance
    from singleton_context_oracle import solve_context
    # The historical lower bound is deliberately not an oracle argument.
    answer=solve_context(instance(payload['instance']),*payload['context'],deadline=deadline-2)
    answer.update(worker_seconds_before_output=time.perf_counter()-origin,
        fixed_original_context=True,historical_lower_bound_imported=False,performance_claim=False)
    put(OUT/'result.json',answer)
    if time.perf_counter()>deadline:raise TimeoutError('whole worker including output over limit')


def launch():
    import psutil
    from revision31.guard import guarded_command
    m=verify();origin=time.perf_counter();limits=deadlines(m,origin);deadline=limits['worker']
    put(OUT/'launch-start.json',dict(origin=origin,deadline=deadline,limits=limits,pid=os.getpid()))
    cmd=[sys.executable,'-B','-X','utf8','-m','singleton_witness','worker','--origin',repr(origin),'--deadline',repr(deadline)]
    peak=0;missing=0;samples=0;cause=None
    with (OUT/'worker.log').open('xb')as log:
        p=subprocess.Popen(guarded_command(cmd,deadline,os.getpid()),cwd=ROOT,stdin=subprocess.DEVNULL,
            stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        while p.poll()is None:
            try:
                proc=psutil.Process(p.pid);rss=sum(q.memory_info().rss for q in [proc]+proc.children(recursive=True))
                peak=max(peak,rss);samples+=1
                if rss>m['memory_bytes']:cause='memory_limit';p.kill()
            except(psutil.NoSuchProcess,psutil.AccessDenied):missing+=1
            if time.perf_counter()>deadline+.3:cause='time_limit';p.kill()
            time.sleep(.02)
        code=p.wait()
    raw=OUT/'result.json';digest=sha(raw)if raw.exists()else None
    # The expensive output hash is inside this charge, not after a 0.1s estimate.
    elapsed=time.perf_counter()-origin
    receipt=dict(returncode=code,outcome=cause or('complete'if code==0 else'error'),
        seconds_through_output_hash=elapsed,raw_sha256=digest,sampled_peak_tree_bytes=peak,
        samples=samples,missing_samples=missing)
    put(OUT/'receipt.json',receipt)
    whole_compute=time.perf_counter()-origin
    charged=time.perf_counter()-m['preparation_origin']+m['source_charge']
    passed=code==0 and cause is None and raw.exists()and whole_compute<=m['compute_seconds']and charged<=m['stage_seconds']
    put(OUT/'complete.json',dict(status='PASS_WITNESS_PENDING_INDEPENDENT_AUDIT'if passed else'STOP',
        charged_seconds=charged,source_charge=m['source_charge'],whole_compute_through_receipt_seconds=whole_compute,
        performance_claim=False))
    print(json.dumps(dict(status='PASS_WITNESS_PENDING_INDEPENDENT_AUDIT'if passed else'STOP',
                         charged_seconds=charged,compute_seconds=whole_compute)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','launch','worker'])
    p.add_argument('--origin',type=float);p.add_argument('--deadline',type=float);a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='launch':launch()
    else:worker(a.origin,a.deadline)
