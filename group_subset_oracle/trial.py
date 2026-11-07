"""Once-only original G5 subset-cost mechanism validation, with all fees."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/revision39_group_subset_20261106'
OLD=ROOT/'results/revision39_compact_full_20261106'


def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,value):
    from review88_physical.evidence import encode
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('xb')as f:f.write(encode(value))


def native_masks(njobs):return tuple(sorted((s for s in range(1<<njobs)if s.bit_count()>=2),key=lambda s:(s.bit_count(),s)))


def limits(m,origin,compute,reserve):
    final=m['preparation_origin']+m['stage_seconds']-m['source_charge']-3
    external=min(origin+compute,final-120);internal=external-reserve
    if internal<=origin:raise TimeoutError('stage budget cannot admit the worker')
    return dict(external=external,internal=internal,final=final)


def bind_source(prior,payload,receipt_path,raw_path):
    if sha(receipt_path)!=prior['original_source_receipt_sha256']:raise ValueError('sealed source receipt drift')
    if sha(raw_path)!=payload['source_raw_sha256']:raise ValueError('sealed source raw drift')
    receipt=read(receipt_path)
    if receipt['returncode']!=0 or receipt['raw_sha256']!=payload['source_raw_sha256']:
        raise ValueError('invalid source record')
    if receipt['charged_seconds']!=prior['source_charge']:raise ValueError('sealed source charge differs')
    return receipt


def bind_artifacts(result,receipt,out):
    prices=result.get('prices');artifacts=result.get('artifacts')
    if not isinstance(prices,list)or prices!=sorted(set(prices))or any(type(v)is not int or not 0<=v<16 for v in prices):
        raise ValueError('invalid price inventory')
    expected={'model.json'}|{f'prices/{v}.json'for v in prices}
    actual={'model.json'}|{p.relative_to(out).as_posix()for p in (out/'prices').glob('*.json')}
    if not isinstance(artifacts,dict)or set(artifacts)!=expected or actual!=expected:raise ValueError('artifact inventory differs')
    if sha(out/'result.json')!=receipt['raw_sha256']:raise ValueError('result receipt differs')
    for name,h in artifacts.items():
        if sha(out/name)!=h:raise ValueError('model or price artifact drift '+name)


def checked_price_matrix(model,record,subset,start=0,end=0):
    base=model['matrix'];keys=model['row_keys'];inst=model['instance'];H=inst['H'];P=len(inst['pads'])
    if len(keys)!=len(base['rhs'])or inst['groups']!=1:raise ValueError('G1 row inventory differs')
    rhs=list(base['rhs']);covers=set();pack=0;boundary=set()
    for i,key in enumerate(keys):
        if key[0]=='cover':
            j=key[1];covers.add(j);rhs[i]=int(j in subset)
            if not base['equalities'][i]:raise ValueError('cover must be equality')
        elif key[0]=='pack':
            if key[1]!=0 or not base['equalities'][i]:raise ValueError('G1 pack equality differs')
            rhs[i]=len(subset);pack+=1
        elif key[0]=='flow'and key[3]in (0,H):
            _,g,p,t=key
            if g!=0 or not base['equalities'][i]:raise ValueError('G1 boundary flow differs')
            boundary.add((p,t));rhs[i]=int(p==start)if t==0 else-int(p==end)
    if covers!=set(range(len(inst['jobs'])))or pack!=1 or boundary!={(p,t)for p in range(P)for t in (0,H)}:
        raise ValueError('required subset rows incomplete')
    if record['actual_rhs']!=rhs:raise ValueError('saved price is for a different cover/pack/endpoint RHS')
    if record['cost_scale']!=model['cost_scale']:raise ValueError('price cost scale differs')
    if record.get('LP')and record['LP'].get('price')is not None and record['LP']['scale']!=model['cost_scale']:
        raise ValueError('LP objective scale differs')
    return dict(base,rhs=rhs)


def verify():
    m=read(OUT/'freeze.json')
    if sha(OUT/'freeze.json')!=(OUT/'freeze.sha256').read_text():raise ValueError('manifest drift')
    if len(m['sources'])<205:raise ValueError('source inventory incomplete')
    for name,h in m['sources'].items():
        if sha(ROOT/name)!=h or sha(OUT/'source_snapshot'/name)!=h:raise ValueError('source drift '+name)
    for name in ['input','acceptance']:
        if sha(OUT/f'{name}.json')!=m[name+'_sha256']:raise ValueError('input drift '+name)
    return m


def prepare():
    if OUT.exists():raise FileExistsError('one frozen table only, no restart')
    origin=time.perf_counter()
    from singleton_witness import verify as old_verify
    old=old_verify();sources=dict(old['sources'])
    receipt=bind_source(old,read(ROOT/'results/revision39_singleton_witness_20261106/input.json'),OLD/'receipts/0.json',OLD/'arms/0.json')
    r=read(OLD/'arms/0.json')
    structural=ROOT/'results/revision39_structural_20261106'
    if sha(structural/'freeze.json')!=(structural/'freeze.sha256').read_text():raise ValueError('source metadata drift')
    metadata=read(structural/'freeze.json')['inputs']['Geo-groups5'];input_row=structural/metadata['path']
    if sha(input_row)!=metadata['sha256']or read(input_row)!=read(OLD/'input.json')['instance']:
        raise ValueError('original G5 identity differs')
    chosen=(0,2,7,1);key=r['local_events'][0]['kernel'];contexts=[];expected=[]
    for i in chosen:
        e=r['local_events'][i]
        if e['kernel']!=key or e['context'][:2]!=[[0]*5,[0]*5]:raise ValueError('fixed context differs')
        contexts.append(dict(source_query=i,context=e['context']))
        expected.append(dict(source_query=i,lower=e['result']['lower_bound'],upper=e['result']['objective']))
    fleet=r['local_events'][2]['result']['fleet'];plan=next(p for p in fleet['group_plans']if p['routes'])
    subset=sorted((q['job']if isinstance(q,dict)else q.job)for q in plan['routes'])
    payload=dict(instance=r['local_models'][key]['instance'],crew_capacity=metadata['metadata']['crew'],
        contexts=contexts,probe=dict(subset=subset,plan=plan,start=0,end=0),
        source_raw_sha256=receipt['raw_sha256'],source_metadata_sha256=sha(structural/'freeze.json'))
    for path in (ROOT/'group_subset_oracle').glob('*.py'):sources[path.relative_to(ROOT).as_posix()]=sha(path)
    for path in (ROOT/'tests').glob('test_group_subset_*.py'):sources[path.relative_to(ROOT).as_posix()]=sha(path)
    for name in ['docs/2026-11-06-revision39-group-subset-protocol.md',
                 'docs/2026-11-06-revision39-group-subset-factorization.md',
                 'docs/2026-11-06-revision39-group-subset-factorization-review.md']:
        sources[name]=sha(ROOT/name)
    OUT.mkdir()
    for name,h in sources.items():
        if sha(ROOT/name)!=h:raise ValueError('source changed before freeze')
        p=OUT/'source_snapshot'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((ROOT/name).read_bytes())
    put(OUT/'input.json',payload);put(OUT/'acceptance.json',dict(contexts=expected))
    m=dict(schema='group-subset-mechanism-v1',preparation_origin=origin,sources=sources,
        input_sha256=sha(OUT/'input.json'),acceptance_sha256=sha(OUT/'acceptance.json'),
        source_charge=receipt['charged_seconds'],stage_seconds=1200,worker_seconds=240,probe_seconds=30,
        per_query_seconds=15,query_output_reserve=1,worker_output_reserve=10,
        memory_bytes=2*1024**3,native_masks=native_masks(4),native_queries_maximum=11,
        threads=1,seed=0,retries=0,new_instances=0,contexts=4,performance_claim=False)
    put(OUT/'freeze.json',m);(OUT/'freeze.sha256').write_text(sha(OUT/'freeze.json'))
    print(json.dumps(dict(status='FROZEN_NO_OPTIMIZATION',sources=len(sources),native_queries_maximum=11)))


def probe(deadline):
    verify();data=read(OUT/'input.json')
    from review88_physical.evidence import instance
    from .model import build_template,embed_plan
    inst=instance(data['instance']);t=build_template(inst,deadline,crew_capacity=data['crew_capacity'])
    s=data['probe'];checks=[]
    if s['subset']!=list(range(4)):raise ValueError('full-task probe source differs')
    for mask in range(16):
        S=tuple(j for j in range(4)if mask>>j&1)
        plan=dict(trail=s['plan']['trail'],routes=[r for r in s['plan']['routes']if r['job']in S],
                  assays=[a for a in s['plan']['assays']if a[0]in S])
        x,physical=embed_plan(t,plan,S,s['start'],s['end'])
        checks.append(dict(mask=mask,objective=physical['objective'],valid=physical['valid']))
    put(OUT/'probe.json',dict(status='PASS',optimizer_calls=0,pricing_groups=t.m.inst.groups,
        shape=list(t.A.shape),subset_source_checks=checks,feasible_masks=list(range(16))))


def worker(deadline):
    m=verify();data=read(OUT/'input.json')
    from review88_physical.evidence import instance
    from service_separator25_joint.kernel import validate_boundary
    from singleton_context_oracle import single_task,_trail
    from .model import validate_subset,elementary
    from .kernel import SubsetKernel
    from .dp import combine_intervals
    origin=time.perf_counter();inst=instance(data['instance']);table={};phases=[];stop=None
    if len(inst.jobs)!=4:raise ValueError('fixed local task total differs')
    q=time.perf_counter();k=SubsetKernel(inst,deadline,crew_capacity=data['crew_capacity'])
    phases.append(dict(stage='G1_matrix_LP_and_proxy_setup',seconds=time.perf_counter()-q))
    try:
        t=k.t
        put(OUT/'model.json',dict(instance=t.m.inst,crew_capacity=t.crew_capacity,cost_scale=t.cost_scale,
            matrix=k.matrix,row_keys=t.m.M.keys,support=dict(routes=t.m.routes,arcs=t.m.arcs,
            rids=[(g,i,col)for(g,i),col in t.m.rids.items()],zids=[(g,i,col)for(g,i),col in t.m.zids.items()],
            assignments=list(t.m.assignments.items()))))
        plan=dict(routes=[],assays=[],trail=_trail(inst,t.roads,0,0))
        table[0]=dict(status='optimal',optimal=True,subset=(),start=0,end=0,
            lower_bound=F(0),objective=F(0),plan=plan,analytic=True,
            physical=validate_subset(inst,plan,(),0,0,crew_capacity=t.crew_capacity))
        put(OUT/'prices/0.json',table[0])
        q=time.perf_counter()
        for j in range(4):
            r=single_task(inst,j,0,0,deadline=deadline)
            if r['status']!='optimal':raise ValueError('unexpected infeasible singleton')
            r.update(subset=(j,),lower_bound=r['objective'],optimal=True,analytic=True,
                physical=validate_subset(inst,r['plan'],(j,),0,0,crew_capacity=t.crew_capacity),
                closure_basis='analytic_single_task_event_and_road_certificate')
            if not r['physical']['valid']:raise ValueError('analytic source physical error')
            table[1<<j]=r;put(OUT/'prices'/f'{1<<j}.json',r)
        phases.append(dict(stage='four_singleton_prices',seconds=time.perf_counter()-q))
        for mask in m['native_masks']:
            started=time.perf_counter()
            if started+1>=deadline:stop='worker_allowance';break
            subset=tuple(j for j in range(4)if mask>>j&1)
            r=k.query(subset,0,0,min(deadline,started+m['per_query_seconds']-m['query_output_reserve']))
            r['analytic']=False;table[mask]=r;put(OUT/'prices'/f'{mask}.json',r)
            if r['status']=='infeasible':raise ValueError('native infeasibility contradicts the validated subset source')
            if time.perf_counter()-started>m['per_query_seconds']:raise TimeoutError('query whole-output allowance exceeded')
            if not r['optimal']and r['status']!='infeasible':stop='first_unresolved_subset';break
        # Missing native costs remain valid intervals, never missing edges.
        complete=len(table)==16 and all(r.get('optimal')or r['status']=='infeasible'for r in table.values())
        defaults={}
        for mask in range(16):
            if mask in table:defaults[mask]=table[mask];continue
            subset=tuple(j for j in range(4)if mask>>j&1)
            lower=elementary(t,subset,0,0)
            defaults[mask]=dict(status='infeasible'if lower is None else'unresolved',
                lower_bound=0 if lower is None else lower,objective=None,plan=None)
        contexts=[];q=time.perf_counter()
        for row in data['contexts']:
            start,end,counts=map(tuple,row['context']);r=combine_intervals([defaults]*inst.groups,counts,4)
            fleet=None
            if r['objective']is not None:
                plans=[defaults[mask]['plan']for mask in r['upper_partition']]
                check=validate_boundary(inst,plans,start,end)
                if not check['valid']or F(check['objective'])!=r['objective']:raise ValueError('DP fleet mismatch')
                if [len(p['routes'])for p in plans]!=list(counts):raise ValueError('DP group labels differ')
                for mask,plan in zip(r['upper_partition'],plans):
                    S=tuple(j for j in range(4)if mask>>j&1)
                    if not validate_subset(inst,plan,S,0,0,crew_capacity=t.crew_capacity)['valid']:
                        raise ValueError('DP original crew/domain mismatch')
                fleet=dict(group_plans=plans,check=check)
            contexts.append(dict(r,source_query=row['source_query'],context=row['context'],fleet=fleet))
        phases.append(dict(stage='four_labeled_partition_DPs_and_physical_replay',seconds=time.perf_counter()-q))
        answer=dict(status='complete'if complete else'partial',reason=stop,contexts=contexts,
            prices=sorted(table),native_prices=sum(not r.get('analytic',False)for r in table.values()),
            lp_calls=k.lp_calls,mip_calls=k.mip_calls,pricing_groups=1,full_joint_G5_models=0,
            native_model_builds=k.mip.builds,context_resets=k.mip.context_resets,
            native_optimize_seconds=k.mip.optimize_seconds,phases=phases,
            seconds_before_output=time.perf_counter()-origin,old_source_costs_imported=False,performance_claim=False)
        names=['model.json']+[f'prices/{mask}.json'for mask in sorted(table)]
        answer['artifacts']={name:sha(OUT/name)for name in names}
        answer['seconds_before_output']=time.perf_counter()-origin
        put(OUT/'result.json',answer)
        if time.perf_counter()>deadline+10:raise TimeoutError('worker output exceeded allowance')
    finally:k.close()


def guarded(role,m):
    import psutil
    from revision31.guard import guarded_command
    origin=time.perf_counter();compute=m['probe_seconds']if role=='probe'else m['worker_seconds']
    bound=limits(m,origin,compute,2 if role=='probe'else m['worker_output_reserve'])
    cmd=[sys.executable,'-B','-X','utf8','-m','group_subset_oracle.trial',role,'--deadline',repr(bound['internal'])]
    put(OUT/f'{role}-start.json',dict(origin=origin,limits=bound,pid=os.getpid()))
    peak=0;samples=0;missing=0;cause=None
    with (OUT/f'{role}.log').open('xb')as log:
        p=subprocess.Popen(guarded_command(cmd,bound['external'],os.getpid()),cwd=ROOT,
            stdin=subprocess.DEVNULL,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        while p.poll()is None:
            try:
                proc=psutil.Process(p.pid);rss=sum(q.memory_info().rss for q in [proc]+proc.children(recursive=True))
                peak=max(peak,rss);samples+=1
                if rss>m['memory_bytes']:cause='memory_limit';p.kill()
            except(psutil.NoSuchProcess,psutil.AccessDenied):missing+=1
            if time.perf_counter()>bound['external']+.3:cause='time_limit';p.kill()
            time.sleep(.02)
        code=p.wait()
    path=OUT/('probe.json'if role=='probe'else'result.json');digest=sha(path)if path.exists()else None
    put(OUT/f'{role}-receipt.json',dict(returncode=code,cause=cause,raw_sha256=digest,
        seconds_through_hash=time.perf_counter()-origin,sampled_peak_tree_bytes=peak,samples=samples,missing_samples=missing))
    elapsed=time.perf_counter()-origin
    if code!=0 or cause or not path.exists()or elapsed>compute:raise ValueError(('worker failed',role,code,cause,elapsed))
    return elapsed


def analyze():
    from copy import deepcopy
    from review88_physical.evidence import instance
    from .model import validate_subset
    from .dp import combine_intervals
    from boundary_context_trial import audit_helper
    data=read(OUT/'input.json');expected=read(OUT/'acceptance.json')['contexts'];result=read(OUT/'result.json')
    bind_artifacts(result,read(OUT/'worker-receipt.json'),OUT)
    inst=instance(data['instance']);model=read(OUT/'model.json');base=model['matrix'];audit=audit_helper()
    prices={int(p.stem):read(p)for p in (OUT/'prices').glob('*.json')};vectors=0;proofs=0;maximum=0.
    for mask,r in prices.items():
        subset=tuple(j for j in range(4)if mask>>j&1)
        if tuple(r['subset'])!=subset:raise ValueError('subset record label differs')
        if r.get('plan')is not None:
            checked=validate_subset(inst,r['plan'],subset,0,0,crew_capacity=data['crew_capacity'])
            if not checked['valid']or F(str(checked['objective']))!=F(str(r['objective'])):raise ValueError('price physical mismatch')
        if r.get('analytic'):continue
        matrix=checked_price_matrix(model,r,subset)
        if r['LP']and r['LP'].get('price')is not None:
            audit.check_price(matrix,r['LP']);proofs+=1
            v=audit.check_vector(matrix,r['LP']['vector'],relax=True);vectors+=1;maximum=max(maximum,v['residual'])
        for key in ['native_vector','vector']:
            if r.get(key)is not None:
                v=audit.check_vector(matrix,r[key],completed=key=='vector');vectors+=1;maximum=max(maximum,v['residual'])
                if key=='vector'and abs(v['objective']/r['cost_scale']-float(F(str(r['objective']))))>1e-6:
                    raise ValueError('completed price objective differs')
    # Independent exhaustive labeled partitions, at most5^4 per saved context.
    low={};up={}
    for mask in range(16):
        if mask in prices:
            r=prices[mask];low[mask]=None if r['status']=='infeasible'else F(str(r['lower_bound']))
            up[mask]=None if r['objective']is None else F(str(r['objective']))
        else:
            # Recover the valid unknown lower used by the worker from independent
            # events, without constructing any matrix or solving a price.
            S=[j for j in range(4)if mask>>j&1];values=[]
            for j in S:
                values.append(min(inst.flight_fixed+inst.air_cost*(r['recovery']-r['launch'])+
                    inst.jobs[j].weight*(r['recovery']+inst.handoff_delay+inst.assay_duration-r['sample'])
                    for r in model['support']['routes']if r['job']==j))
            low[mask]=sum(values,F(0));up[mask]=None
    rows=[]
    for got,reference in zip(result['contexts'],expected):
        if got['source_query']!=reference['source_query']:raise ValueError('context order differs')
        counts=tuple(got['context'][2]);lower_values=[];upper_values=[]
        for allocation in product(range(inst.groups),repeat=4):
            if tuple(allocation.count(g)for g in range(inst.groups))!=counts:continue
            masks=tuple(sum(1<<j for j,g in enumerate(allocation)if g==i)for i in range(inst.groups))
            if all(low[s]is not None for s in masks):lower_values.append(sum((low[s]for s in masks),F(0)))
            if all(up[s]is not None for s in masks):upper_values.append(sum((up[s]for s in masks),F(0)))
        L=min(lower_values)if lower_values else None;U=min(upper_values)if upper_values else None
        if (L,U)!=(None if got['lower_bound']is None else F(str(got['lower_bound'])),
                   None if got['objective']is None else F(str(got['objective']))):raise ValueError('partition DP mismatch')
        if L is not None and L>reference['upper']or U is not None and U<reference['lower']:
            raise ValueError('new interval contradicts the saved source certificate')
        if got['optimal']and not reference['lower']<=U<=reference['upper']:raise ValueError('closed value mismatch')
        if got['fleet']is not None:
            from service_separator25_joint.kernel import validate_boundary
            start,end,_=map(tuple,got['context']);check=validate_boundary(inst,got['fleet']['group_plans'],start,end)
            if not check['valid']or F(check['objective'])!=U:raise ValueError('joint fleet replay mismatch')
        rows.append(dict(source_query=got['source_query'],lower=L,upper=U,optimal=got['optimal'],
            previous_lower=reference['lower'],previous_upper=reference['upper']))
    passed=result['status']=='complete'and all(x['optimal']for x in rows)
    return dict(status='PASS_FULL_DOMAIN_MECHANISM'if passed else'PARTIAL_VALID_INTERVALS_NO_EXTENSION',
        rows=rows,prices=len(prices),vectors_checked=vectors,LP_prices_recomputed=proofs,
        max_residual=maximum,original_joint_models_built=0,performance_claim=False)


def run():
    m=verify();stop=None;status='STOP';analysis=None
    try:
        guarded('probe',m)
        if read(OUT/'probe.json')['status']!='PASS':raise ValueError('preflight incomplete')
        guarded('worker',m);verify();analysis=analyze();put(OUT/'analysis.json',analysis);status=analysis['status']
    except BaseException as exc:stop=repr(exc)
    charge=time.perf_counter()-m['preparation_origin']+m['source_charge']
    if charge>m['stage_seconds']:status='STOP';stop=stop or'global charged allowance exceeded'
    put(OUT/'complete.json',dict(status=status,stop=stop,charged_seconds=charge,
        source_charge=m['source_charge'],performance_claim=False))
    print(json.dumps(dict(status=status,stop=stop,charged_seconds=charge)))
    if status=='STOP':raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','probe','worker'])
    p.add_argument('--deadline',type=float);a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='run':run()
    elif a.action=='probe':probe(a.deadline)
    else:worker(a.deadline)
