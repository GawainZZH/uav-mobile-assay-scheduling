import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest


class PortableContract(unittest.TestCase):
    def test_manifest_rejects_parent_paths_before_reading_data(self):
        from reproduction40.portable import verify
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)
            (p/'PUBLIC_MANIFEST.json').write_text(json.dumps({'files':{'../outside':'0'*64}}))
            with self.assertRaisesRegex(ValueError,'relative'): verify(p)

    def test_duplicate_json_keys_are_rejected(self):
        from reproduction40.portable import load
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'data.json';p.write_text('{"a":1,"a":2}')
            with self.assertRaisesRegex(ValueError,'duplicate'):load(p)

    def test_hash_drift_cannot_pass(self):
        from reproduction40.portable import verify
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);(p/'a').write_bytes(b'original')
            (p/'PUBLIC_MANIFEST.json').write_text(json.dumps({'files':{'a':hashlib.sha256(b'original').hexdigest()}}))
            self.assertEqual(verify(p)['verified_files'],1)
            (p/'a').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'hash'):verify(p)

    def test_existing_output_is_not_overwritten(self):
        from reproduction40.portable import publish
        with tempfile.TemporaryDirectory()as folder:
            p=Path(folder)/'out.json';p.write_bytes(b'original')
            with self.assertRaises(FileExistsError):publish(p,{'new':True})
            self.assertEqual(p.read_bytes(),b'original')

    def test_exact_matrix_keys_and_rational_values_have_portable_json(self):
        from fractions import Fraction
        from reproduction40.portable import publish,load
        with tempfile.TemporaryDirectory()as folder:
            p=Path(folder)/'result.json';publish(p,{'prices':{(1,2):Fraction(1,3)},'upper':Fraction(7)})
            self.assertEqual(load(p),{'prices':{'(1, 2)':'1/3'},'upper':7})

    def test_all_48_receipts_are_required_and_proof_timing_is_consistent(self):
        from reproduction40.portable import aggregate
        rows=[]
        for i in range(48):
            method='F'if i%2==0 else'S'
            rows.append({'arm':{'case':'c'+str(i//2),'method':method},
                         'assessment':{'eligible':True,'solved':method=='F','par2':10 if method=='F' else 420},
                         'seconds_through_hash':10,'sampled_peak_bytes':1024**2,'error':None,'cause':None})
        a=aggregate(rows)
        self.assertEqual(a['methods']['F']['proven_optima'],24)
        self.assertEqual(a['methods']['S']['mean_PAR2'],420)
        with self.assertRaisesRegex(ValueError,'48'):aggregate(rows[:-1])
        rows[0]['assessment']['par2']=1
        with self.assertRaisesRegex(ValueError,'proof time'):aggregate(rows)

    def test_invalid_method_and_nonfinite_budget_fail_without_scientific_import(self):
        from reproduction40.portable import arguments
        with self.assertRaises(ValueError):arguments('c','unknown',20)
        with self.assertRaises(ValueError):arguments('c','F',float('nan'))

    def test_only_known_resource_exceptions_can_become_unresolved(self):
        from reproduction40.portable import resource_reason
        self.assertEqual(resource_reason(TimeoutError()),'internal_time_limit')
        self.assertEqual(resource_reason(ValueError('unexpected raw SCIP status memlimit')),'native_memory_limit')
        self.assertIsNone(resource_reason(ValueError('infeasible source')))

    def test_run_entry_retains_only_physically_valid_checkpoint_on_timeout(self):
        from unittest.mock import patch
        from reproduction40.portable import run_one,load
        from stock_scale.cases import make
        from stock_execution.entry import checkpoint
        r=make(dict(name='DEMO',G=2,P=2,K=2,n=2,B=2,extended=False),11)
        def limited(record,method,deadline,save):
            save(checkpoint(r['raw'],r['candidate']));raise TimeoutError()
        with tempfile.TemporaryDirectory()as folder:
            out=Path(folder)/'out.json'
            with patch('reproduction40.portable.verify'),patch('reproduction40.portable.record',return_value=r),patch('stock_execution.entry.solve_record',side_effect=limited):
                a=run_one(folder,'DEMO','F',20,out)
            self.assertEqual(a['status'],'resource_limit');self.assertTrue(a['physical_valid'])
            self.assertIsNotNone(load(out)['result']['interval']['upper'])

    def test_run_entry_does_not_publish_a_bug_as_valid_result(self):
        from unittest.mock import patch
        from reproduction40.portable import run_one
        with tempfile.TemporaryDirectory()as folder:
            out=Path(folder)/'out.json'
            with patch('reproduction40.portable.verify'),patch('reproduction40.portable.record',return_value={}),patch('stock_execution.entry.solve_record',side_effect=ValueError('real bug')):
                with self.assertRaisesRegex(ValueError,'real bug'):run_one(folder,'DEMO','F',20,out)
            self.assertFalse(out.exists())

    def test_shipped_receipt_chain_reproduces_original_summary(self):
        from reproduction40.portable import summarize
        root=Path(__file__).resolve().parents[1]
        if not(root/'PUBLIC_MANIFEST.json').exists():root=root/'public_release/uav-mobile-assay-stock-flow'
        if not(root/'PUBLIC_MANIFEST.json').exists():self.skipTest('public archive not exported')
        a=summarize(root)
        self.assertEqual(a['methods']['F']['proven_optima'],18)
        self.assertEqual(a['methods']['S']['proven_optima'],4)
        self.assertAlmostEqual(a['methods']['F']['mean_PAR2'],118.4740881875041)


if __name__=='__main__':unittest.main()
