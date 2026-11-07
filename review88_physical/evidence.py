"""Immutable finite input/source manifest and exact portable evidence."""
from dataclasses import asdict,is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from .types import Instance,Job
from .cases import cases,tiny_cases,witness

ROOT=Path(__file__).resolve().parents[1]
MODULES=('__init__','types','routes','graph','physical','milp','pooled','cases','evidence','run','audit')


def plain(value):
    if is_dataclass(value):return plain(asdict(value))
    if isinstance(value,Fraction):return value.numerator if value.denominator==1 else str(value)
    if isinstance(value,dict):return {str(k):plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [plain(v) for v in value]
    if hasattr(value,'item'):return value.item()
    return value


def encode(value):
    return (json.dumps(plain(value),ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n').encode('utf8')


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(encode(value))


def read(path):return json.loads(Path(path).read_text(encoding='utf8'))


def instance(data):
    data=dict(data);data['jobs']=tuple(Job(**j) for j in data['jobs'])
    return Instance(**data)


def sources():
    from tools.package_referee75_revision import _local_import_closure
    members={f'review88_physical/{name}.py' for name in MODULES}
    members.update(('tests/test_review88_physical.py','docs/2026-10-19-review88-bridge-design.md',
                    'docs/plans/2026-10-19-review88-revision.md','tools/package_referee75_revision.py'))
    _local_import_closure(ROOT,members)
    return {name:sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(members)}


def freeze(directory,tiny=False):
    from .physical import validate_fleet
    import platform
    import highspy
    directory=Path(directory)
    if directory.exists():raise ValueError('freeze exists; replacement forbidden')
    inputs=tiny_cases() if tiny else cases();hashes={};arms=[]
    for index,inst in enumerate(inputs):
        path=directory/'inputs'/(inst.name+'.json');write(path,inst)
        hashes[path.name]=sha256(path.read_bytes()).hexdigest()
        p=witness(inst);v=validate_fleet(inst,p['group_plans'])
        if not v['valid']:raise ValueError(('pre-freeze witness',inst.name,v))
        write(directory/'witnesses'/(inst.name+'.json'),dict(fleet=p,physical=v))
        methods=('PF','PQ') if tiny or inst.name.startswith('B_') else ('J','TP','QC')
        if inst.name.startswith('B_') and index%2:methods=('PQ','PF')
        for method in methods:arms.append(dict(case=inst.name,method=method,
            policy=method if method in ('J','TP','QC') else 'J',
            seconds=10 if tiny else 45 if inst.name.startswith('B_') else 30,memory_GiB=2))
    source=sources()
    for name in source:
        path=directory/'source_snapshot'/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes((ROOT/name).read_bytes())
    manifest=dict(schema='review88-physical-v1',tiny=tiny,inputs=hashes,sources=source,arms=arms,
        serial=True,no_replacement=True,graph_state_limit=200000,graph_arc_limit=800000,
        global_wall_budget=120 if tiny else 1920,field_calibrated=False,
        protocol='ell1 with unchanged q includes conservative additional transfer reservation',
        machine=dict(platform=platform.platform(),python=platform.python_version(),highs=highspy.HIGHS_VERSION_MAJOR))
    write(directory/'freeze.json',manifest);return manifest


def verify(directory):
    directory=Path(directory);m=read(directory/'freeze.json')
    for name,h in m['sources'].items():
        if sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('source drift: '+name)
        if sha256((directory/'source_snapshot'/name).read_bytes()).hexdigest()!=h:
            raise ValueError('source snapshot drift: '+name)
    for name,h in m['inputs'].items():
        if sha256((directory/'inputs'/name).read_bytes()).hexdigest()!=h:raise ValueError('input drift: '+name)
    return m
