"""Standard-library access to content-addressed scientific evidence."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import urllib.request
import zipfile


def relative(s):
    if not isinstance(s,str)or not s or'\\'in s or':'in s:raise ValueError('canonical relative path required')
    p=PurePosixPath(s)
    if p.is_absolute()or p.as_posix()!=s or any(v in('..','.')for v in p.parts):
        raise ValueError('canonical relative path required')
    return s


def admissible(name,roles):
    relative(name);p=PurePosixPath(name);base=p.name.lower()
    if any(v in('.git','.codex','reviews')for v in p.parts):return False
    # Verified numerical rows/envelope replay report, despite its legacy name.
    if name=='docs/data/referee_evidence21_20261021/numerical/referee-report.json'and'scientific_summary'in roles:return True
    if re.search(r'reviewer[-_ ]?report|referee[-_ ]?report|writing[-_]guide|language[-_]scan|approval|authorization|permission',base):return False
    if base in('review_start_here.md','agents.md','profile.md','active.md'):return False
    if p.suffix.lower()in('.lic','.pem','.key','.env','.tex','.cls','.sty'):return False
    if name.startswith('manuscript/'):
        return p.suffix=='.pdf'and bool(set(roles)&{'active_figure','plot_geometry_or_composite_support'})
    if p.suffix.lower()=='.md':return base=='unsafe-prefixes.md'
    return p.suffix.lower()in('.json','.jsonl','.csv','.txt','.py','.npz','.npy','.pdf','.stdout','.stderr','.log','.sha256')


def redact(raw):
    for pat in [rb'gh[pousr]_[A-Za-z0-9]{24,}',rb'sk-[A-Za-z0-9]{24,}',
                rb'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----']:
        if re.search(pat,raw):raise ValueError('credential-like content requires exclusion before publication')
    pattern=rb'(?i)(?:[a-z]:[\\/]+users[\\/]+lenovo|[a-z]:[\\/]+pydir[\\/]+b&b&p)'
    return re.subn(pattern,b'ANONYMOUS_WORKSPACE',raw)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()


def validate_catalog(c):
    if not isinstance(c['files'],dict)or not c['files']:raise ValueError('nonempty evidence catalog required')
    def valid_hash(h):
        if not isinstance(h,str)or not re.fullmatch('[0-9a-f]{64}',h):raise ValueError('invalid SHA256')
    for n,v in c['parts'].items():
        relative(n)
        if len(PurePosixPath(n).parts)!=1 or not n.endswith('.zip'):raise ValueError('part filename required')
        valid_hash(v['sha256'])
    used=set()
    for n,v in c['files'].items():
        relative(n);valid_hash(v['sha256']);used.add(v['sha256'])
    for cohort,versions in c.get('source_versions',{}).items():
        for n,h in versions.items():
            relative(n)
            if not n.endswith('.py'):raise ValueError('only scientific Python source versions')
            valid_hash(h);used.add(h)
    for n,h in c.get('helpers',{}).items():
        relative(n);valid_hash(h);used.add(h)
    for alias,target in c.get('collection_aliases',{}).items():
        if not isinstance(alias,str)or not isinstance(target,str)or not alias or not target:raise ValueError('invalid collection alias')
        seen={alias}
        while target in c.get('collection_aliases',{}):
            if target in seen:raise ValueError('cyclic collection alias')
            seen.add(target);target=c['collection_aliases'][target]
        if not any(target in v['collections']for v in c['files'].values()):raise ValueError('empty collection alias target')
    if used!=set(c['objects']):raise ValueError('unreferenced or missing scientific object')
    for h,v in c['objects'].items():
        valid_hash(h)
        if v['part']not in c['parts']or type(v['bytes'])is not int or v['bytes']<0:raise ValueError('object descriptor differs')


class Store:
    def __init__(self,cache,catalog):
        validate_catalog(catalog);self.cache=Path(cache);self.c=catalog;self.opened={}

    def __enter__(self):return self

    def __exit__(self,*args):
        for z in self.opened.values():z.close()

    def selected(self,collection=None):
        if collection is None:return dict(self.c['files'])
        collection=self.canonical(collection)
        rows={n:v for n,v in self.c['files'].items()if collection in v['collections']}
        if not rows:raise ValueError('unknown or empty collection')
        return rows

    def canonical(self,collection):
        aliases=self.c.get('collection_aliases',{})
        while collection in aliases:collection=aliases[collection]
        return collection

    def needed(self,collection=None):
        collection=self.canonical(collection)
        hs={v['sha256']for v in self.selected(collection).values()}
        hs.update(self.c.get('helpers',{}).values())
        hs.update(v['sha256']for n,v in self.c['files'].items()if n.endswith('.py')and'code'in v['collections'])
        if collection is None:
            hs.update(self.c['objects'])
        else:hs.update(self.c.get('source_versions',{}).get(collection,{}).values())
        return hs

    def fetch(self,collection=None):
        self.cache.mkdir(parents=True,exist_ok=True);parts={self.c['objects'][h]['part']for h in self.needed(collection)}
        for n in sorted(parts):
            p=self.cache/n;v=self.c['parts'][n]
            if p.exists():
                if digest(p)!=v['sha256']:raise ValueError('cached part hash differs')
                continue
            url=v.get('url','')
            if not url.startswith('https://github.com/GawainZZH/uav-mobile-assay-scheduling/releases/download/'):
                raise ValueError('unexpected evidence download origin')
            with urllib.request.urlopen(url,timeout=120)as response,p.open('xb')as f:
                for b in iter(lambda:response.read(1024**2),b''):f.write(b)
            if digest(p)!=v['sha256']:raise ValueError('download hash differs')
        return {'downloaded_or_existing_parts':len(parts)}

    def archive(self,part):
        if part not in self.opened:
            p=self.cache/part;v=self.c['parts'][part]
            if not p.is_file():raise FileNotFoundError('fetch required evidence part: '+part)
            if digest(p)!=v['sha256']:raise ValueError('part hash differs')
            z=zipfile.ZipFile(p);expected={'objects/'+h for h,v in self.c['objects'].items()if v['part']==part}
            if len(z.namelist())!=len(expected)or set(z.namelist())!=expected:
                z.close();raise ValueError('archive member population differs')
            self.opened[part]=z
        return self.opened[part]

    def verify(self,collection=None):
        hs=self.needed(collection)
        for h in sorted(hs):
            v=self.c['objects'][h];z=self.archive(v['part']);actual=hashlib.sha256();size=0
            with z.open('objects/'+h)as f:
                for b in iter(lambda:f.read(1024**2),b''):actual.update(b);size+=len(b)
            if actual.hexdigest()!=h or size!=v['bytes']:raise ValueError('scientific object hash differs')
        return {'status':'PASS','logical_files':len(self.selected(collection)),'objects':len(hs),'optimizer_calls':0}

    def materialize(self,destination,collection):
        collection=self.canonical(collection)
        destination=Path(destination)
        if destination.exists():raise FileExistsError('restore destination must be new')
        self.verify(collection);selected=self.selected(collection)
        selected.update({n:v for n,v in self.c['files'].items()if n.endswith('.py')and'code'in v['collections']})
        for n,h in self.c.get('source_versions',{}).get(collection,{}).items():selected[n]={'sha256':h}
        for n,h in self.c.get('helpers',{}).items():selected[n]={'sha256':h}
        destination.mkdir(parents=True,exist_ok=False);root=destination.resolve()
        for n,v in selected.items():
            p=destination/relative(n)
            if not p.resolve().is_relative_to(root):raise ValueError('restore path escapes destination')
            p.parent.mkdir(parents=True,exist_ok=True);h=v['sha256'];z=self.archive(self.c['objects'][h]['part'])
            with z.open('objects/'+h)as source,p.open('xb')as output:
                for b in iter(lambda:source.read(1024**2),b''):output.write(b)
            if digest(p)!=h:raise ValueError('restored hash differs')
        return {'status':'PASS','files':len(selected),'output':str(destination),'optimizer_calls':0}


def main(root):
    from reproduction40.portable import verify,unique
    root=Path(root);verify(root)
    with gzip.open(root/'legacy/catalog.json.gz','rt',encoding='utf8')as f:c=json.load(f,object_pairs_hook=unique)
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['list','fetch','verify','materialize'])
    parser.add_argument('--collection');parser.add_argument('--cache',type=Path,default=root/'legacy-cache')
    parser.add_argument('--output',type=Path);a=parser.parse_args()
    with Store(a.cache,c)as s:
        if a.command=='list':r={'collections':c['collections'],'logical_files':len(c['files']),'objects':len(c['objects'])}
        elif a.command=='fetch':r=s.fetch(a.collection)
        elif a.command=='verify':r=s.verify(a.collection)
        else:
            if not a.collection or not a.output:parser.error('materialize requires --collection and --output')
            r=s.materialize(a.output,a.collection)
    print(json.dumps(r,ensure_ascii=False))
