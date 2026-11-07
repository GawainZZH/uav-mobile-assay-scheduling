import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import zipfile


class LegacyArchive(unittest.TestCase):
    def sample(self,root):
        raw=b'{"objective":56}';h=hashlib.sha256(raw).hexdigest()
        part=root/'part0.zip'
        with zipfile.ZipFile(part,'w')as z:z.writestr('objects/'+h,raw)
        c={'files':{'results/A/arms/0.json':{'sha256':h,'collections':['A']},
                    'results/B/arms/0.json':{'sha256':h,'collections':['B']}},
           'objects':{h:{'part':'part0.zip','bytes':len(raw)}},
           'parts':{'part0.zip':{'sha256':hashlib.sha256(part.read_bytes()).hexdigest(),'bytes':part.stat().st_size}},
           'source_versions':{}}
        return c,part

    def test_paths_reject_traversal_and_windows_drives(self):
        from legacy41.store import relative
        for s in ['../x','/tmp/x','C:/Users/x','a\\b','a/../b']:
            with self.assertRaises(ValueError):relative(s)
        self.assertEqual(relative('results/A/a.json'),'results/A/a.json')

    def test_review_and_credentials_are_not_scientific_members(self):
        from legacy41.store import admissible,redact
        self.assertFalse(admissible('docs/data/X/reviewer-report.md',[]))
        self.assertFalse(admissible('gurobi.lic',[]))
        self.assertTrue(admissible('results/X/arms/0.json',['native_evidence']))
        self.assertTrue(admissible('docs/data/referee_evidence21_20261021/numerical/referee-report.json',['scientific_summary']))
        with self.assertRaises(ValueError):redact(b'ghp_'+b'A'*32)

    def test_path_redaction_preserves_json_numbers_and_other_fields(self):
        from legacy41.store import redact
        a={'objective':56,'bound':2.5,'valid':True,'path':'C:/Users/Lenovo/work/run.json'}
        b,n=redact(json.dumps(a).encode())
        v=json.loads(b);self.assertEqual(n,1)
        self.assertEqual({k:v[k]for k in a if k!='path'},{k:a[k]for k in a if k!='path'})
        self.assertNotIn('Lenovo',v['path'])

    def test_duplicate_logical_data_stores_one_object_and_restores_requested_collection(self):
        from legacy41.store import Store
        with tempfile.TemporaryDirectory()as folder:
            r=Path(folder);c,p=self.sample(r)
            with Store(r,c)as s:
                self.assertEqual(s.verify()['objects'],1)
                out=r/'new';self.assertEqual(s.materialize(out,'A')['files'],1)
                self.assertEqual(json.loads((out/'results/A/arms/0.json').read_bytes())['objective'],56)
                self.assertFalse((out/'results/B/arms/0.json').exists())
                with self.assertRaises(FileExistsError):s.materialize(out,'A')

    def test_archive_drift_and_unreferenced_hidden_members_fail(self):
        from legacy41.store import Store
        with tempfile.TemporaryDirectory()as folder:
            r=Path(folder);c,p=self.sample(r)
            with zipfile.ZipFile(p,'a')as z:z.writestr('private-review.md',b'private')
            c['parts']['part0.zip']['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError,'member'):
                with Store(r,c)as s:s.verify()

    def test_alias_restores_matching_source_version_and_portable_helper(self):
        from legacy41.store import Store
        with tempfile.TemporaryDirectory()as folder:
            r=Path(folder);c,p=self.sample(r);h=next(iter(c['objects']))
            c['collection_aliases']={'alias-A':'A'}
            c['source_versions']={'A':{'scientific.py':h}}
            c['helpers']={'reproducibility/helper.py':h}
            with Store(r,c)as s:
                out=r/'new';s.materialize(out,'alias-A')
                self.assertTrue((out/'scientific.py').is_file())
                self.assertTrue((out/'reproducibility/helper.py').is_file())
                self.assertTrue((out/'results/A/arms/0.json').is_file())



if __name__=='__main__':unittest.main()
