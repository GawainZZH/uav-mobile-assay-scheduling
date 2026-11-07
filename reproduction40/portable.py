"""Public reproduction adapter. Scientific algorithm files are not modified."""
import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import platform
import re
import statistics
import sys
import time


def unique(pairs):
    d={}
    for k,v in pairs:
        if k in d:raise ValueError('duplicate JSON key: '+k)
        d[k]=v
    return d


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'),object_pairs_hook=unique)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()


def verify(root):
    root=Path(root).resolve();m=load(root/'PUBLIC_MANIFEST.json');files=m['files']
    if not isinstance(files,dict)or not files:raise ValueError('nonempty file inventory required')
    targets=[]
    # All path/schema checks precede artifact reads and scientific imports.
    for n,h in files.items():
        if not isinstance(n,str)or'\\'in n:raise ValueError('canonical relative path required')
        p=PurePosixPath(n)
        if p.is_absolute()or p.as_posix()!=n or any(x in('..','.')for x in p.parts)or':'in n:
            raise ValueError('canonical relative path required')
        resolved=(root/n).resolve()
        if not resolved.is_relative_to(root):raise ValueError('relative path escapes root')
        if not isinstance(h,str)or not re.fullmatch('[0-9a-f]{64}',h):raise ValueError('invalid hash')
        targets.append((n,resolved,h))
    for n,p,h in targets:
        if not p.is_file()or digest(p)!=h:raise ValueError('hash drift: '+n)
    return {'status':'PASS','verified_files':len(files)}


def publish(path,value):
    from dataclasses import asdict,is_dataclass
    from fractions import Fraction
    def convert(x):
        if isinstance(x,Fraction):return x.numerator if x.denominator==1 else str(x)
        if is_dataclass(x):return convert(asdict(x))
        if isinstance(x,dict):
            out={}
            for k,v in x.items():
                key=str(k)
                if key in out:raise ValueError('JSON key conversion collision')
                out[key]=convert(v)
            return out
        if isinstance(x,(list,tuple)):return [convert(v)for v in x]
        if hasattr(x,'tolist'):return convert(x.tolist())
        if x is None or isinstance(x,(str,int,float,bool)):return x
        raise TypeError(type(x).__name__)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(convert(value),ensure_ascii=False,allow_nan=False).encode('utf-8')
    with path.open('xb')as f:f.write(encoded)


def arguments(case,method,seconds):
    if method not in('F','S'):raise ValueError('method must be F (SFBC) or S (MILP)')
    if not math.isfinite(seconds)or not 0<seconds<=3600:raise ValueError('finite positive budget required')
    if case!='DEMO'and not re.fullmatch(r'RSC_(A|G3|G6|P2|P4|K2|K6|N2|N6|B4|M16|M24)_s(11|29)',case):
        raise ValueError('unknown case')


def resource_reason(exc):
    if isinstance(exc,TimeoutError):return 'internal_time_limit'
    if isinstance(exc,MemoryError):return 'python_memory_limit'
    if isinstance(exc,ValueError)and str(exc)=='unexpected raw SCIP status memlimit':return 'native_memory_limit'
    if type(exc).__name__=='DomainLimit'and isinstance(exc.args[0],tuple)and exc.args[0][0]=='stock family/arc limit':
        return 'enumeration_limit'
    return None


def aggregate(rows):
    if len(rows)!=48:raise ValueError('all 48 receipts required')
    groups={'F':[],'S':[]};identities=set()
    for r in rows:
        a=r['arm'];method=a['method'];key=(a['case'],method)
        if method not in groups or key in identities:raise ValueError('duplicate or invalid arm')
        identities.add(key);v=r['assessment'];sec=r['seconds_through_hash'];par2=v['par2']
        if r.get('error')is not None or v.get('eligible')is not True:raise ValueError('ineligible receipt')
        if not math.isfinite(sec)or sec<0 or not math.isfinite(par2):raise ValueError('invalid duration')
        if v['solved']:
            if sec>210 or par2!=sec:raise ValueError('proof time inconsistent')
        elif par2!=420:raise ValueError('unproved arm must retain PAR2 penalty')
        groups[method].append(r)
    if any(len(g)!=24 for g in groups.values())or {r['arm']['case']for r in groups['F']}!={r['arm']['case']for r in groups['S']}:
        raise ValueError('24 matched cases required')
    methods={}
    for method,rs in groups.items():
        peaks=[r['sampled_peak_bytes']/1024**2 for r in rs]
        methods[method]={'cases':24,'proven_optima':sum(r['assessment']['solved']for r in rs),
                         'mean_PAR2':sum(r['assessment']['par2']for r in rs)/24,
                         'median_peak_MiB':statistics.median(peaks),'max_peak_MiB':max(peaks)}
    ratio=methods['F']['mean_PAR2']/methods['S']['mean_PAR2']
    return {'methods':methods,'PAR2_reduction_percent':100*(1-ratio),'new_optimizer_calls':0}


def summarize(root):
    root=Path(root);verify(root)
    rows=[load(root/f'data/receipts/{i}.json')for i in range(48)]
    from stock_scale.contract import arms,assess
    for i,arm in enumerate(arms()):
        summary=load(root/f'data/summaries/{i}.json');receipt=rows[i]
        if receipt['arm']!=arm or summary['arm']!=arm:raise ValueError('fixed arm identity differs')
        if assess(summary,receipt['seconds_through_hash'])!=receipt['assessment']:
            raise ValueError('receipt assessment differs from result')
        payload={'result':summary,'error':None,'full_result_sha256':receipt['raw_sha256']}
        encoded=(json.dumps(payload,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf8')
        if hashlib.sha256(encoded).hexdigest()!=receipt['retained_sha256']:
            raise ValueError('retained result binding differs')
    answer=aggregate(rows)
    expected=load(root/'data/analysis.json')
    for method in('F','S'):
        for field in('cases','proven_optima','mean_PAR2'):
            if not math.isclose(answer['methods'][method][field],expected['methods'][method][field],abs_tol=1e-9,rel_tol=1e-12):
                raise ValueError('published arithmetic differs')
    return answer


def replay(root):
    import zipfile
    root=Path(root);summarize(root)
    from stock_execution.entry import check_result
    from stock_scale.contract import arms
    expected=arms();checked=0
    with zipfile.ZipFile(root/'data/raw-results.zip')as z:
        names=[f'arms/{i}.json'for i in range(48)]
        if len(z.namelist())!=48 or set(z.namelist())!=set(names):raise ValueError('raw archive population differs')
        for i,arm in enumerate(expected):
            h=hashlib.sha256()
            with z.open(names[i])as f:
                for b in iter(lambda:f.read(1024**2),b''):h.update(b)
            if h.hexdigest()!=load(root/f'data/receipts/{i}.json')['raw_sha256']:
                raise ValueError('raw result binding differs')
            raw=load(root/f'data/inputs/{arm["case"]}.json');s=load(root/f'data/summaries/{i}.json')
            if not check_result(raw,s)['valid']:raise ValueError('physical result differs')
            checked+=1
    return {'status':'PASS','raw_results_verified':checked,'physical_results_replayed':checked,
            'new_optimizer_calls':0,'proof_scope':'Physical plans and archived-result binding, not a rational proof of the numerical integer tree.'}


def regenerate(root):
    root=Path(root);verify(root)
    from stock_scale.cases import configurations,make
    checked=0
    for c in configurations():
        for seed in(11,29):
            r=make(c,seed);name=r['raw']['instance']['name']
            if r['raw']!=load(root/f'data/inputs/{name}.json'):raise ValueError('raw regeneration differs')
            if r['candidate']['group_plans']!=load(root/f'data/candidates/{name}.json')['group_plans']:
                raise ValueError('candidate regeneration differs')
            checked+=1
    return {'status':'PASS','regenerated_cases':checked,'new_optimizer_calls':0}


def record(root,case):
    from stock_scale.cases import configurations,make
    if case=='DEMO':return make(dict(name='DEMO',G=2,P=2,K=2,n=2,B=2,extended=False),11)
    name,seed=case[4:].rsplit('_s',1)
    r=make(next(c for c in configurations()if c['name']==name),int(seed))
    if r['raw']!=load(Path(root)/f'data/inputs/{case}.json'):raise ValueError('cold generation differs')
    if r['candidate']['group_plans']!=load(Path(root)/f'data/candidates/{case}.json')['group_plans']:
        raise ValueError('cold candidate differs')
    return r


def run_one(root,case,method,seconds,output):
    began=time.perf_counter();arguments(case,method,seconds)
    if Path(output).exists():raise FileExistsError('output already exists')
    verify(root)
    from stock_execution.entry import solve_record,check_result
    r=record(root,case);latest=[]
    # Same 15-second reserve at the manuscript's 180-second internal allowance.
    reserve=min(15,seconds/4);deadline=began+seconds-reserve
    def save(x):
        latest[:]=[x]
    try:
        answer=solve_record(r,method,deadline,save)
    except Exception as exc:
        cause=resource_reason(exc)
        if cause is None:raise
        answer=dict(latest[-1],status='resource_limit',reason=cause)if latest else dict(
            status='resource_limit',objective=None,interval={'lower':0,'upper':None,'closed':False},fleet=None)
    checked=check_result(r['raw'],answer)
    result={'case':case,'method':method,'result':answer,'physical_validation':checked,
            'seconds_before_serialization':time.perf_counter()-began,'requested_internal_seconds':seconds,
            'native_reserve_seconds':reserve,'python':platform.python_version(),
            'platform':platform.platform(),'performance_claim':False,
            'clock_scope':'Portable single process; not the original guarded batch or its performance observations.'}
    publish(output,result)
    return {'status':answer['status'],'objective':answer['objective'],'physical_valid':checked['valid'],
            'elapsed_through_publication':time.perf_counter()-began,'output':str(output)}


def main(root=None):
    root=Path(root)if root is not None else Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for name in('verify','regenerate','summarize','replay'):sub.add_parser(name)
    run=sub.add_parser('run');run.add_argument('--case',required=True);run.add_argument('--method',choices=['F','S'],required=True)
    run.add_argument('--seconds',type=float,default=180);run.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if a.command=='run':answer=run_one(root,a.case,a.method,a.seconds,a.output)
    else:answer={'verify':verify,'regenerate':regenerate,'summarize':summarize,'replay':replay}[a.command](root)
    print(json.dumps(answer,default=str,ensure_ascii=False,allow_nan=False))


if __name__=='__main__':main()
