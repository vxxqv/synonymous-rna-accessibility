"""Invented group data only. No biological file inputs."""
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
import analyze as A
import verify as V

def row(values):
    return {'assays':{a:{'groups':[{'bio':g+1,'samples':[a+str(g+1)],'values':[float(x)],'technical_finite':1,'median':float(x)}
                                        for g,x in enumerate(v)]} for a,v in values.items()}}

def pair(values, other=None):
    other=other or {a:[0]*4 for a in ('R','D','F')}
    return {'code':{'a':row(values),'b':row(other)},'geo':{'a':row(values),'b':row(other)}}

META={'a':'a','b':'b','site':1,'ROI':1,'both_SNV':True}

class Tests(unittest.TestCase):
    def reference(self,inputs):
        r=A.pair_result(META,inputs)
        source={m:list(rows.values()) for m,rows in inputs.items()}
        for m in source:
            for key,rw in zip(('a','b'),source[m]):rw['hgvs_nt']=key
        ref=V.independent_pair(META,source)
        self.assertEqual(r['category'],ref['category'])
        self.assertEqual(r['unstructured_log_radius_exclusive'],None if ref['radius'] is None else str(ref['radius']))
        return r

    def test_opposing(self):
        r=self.reference(pair({'R':[2]*4,'D':[0]*4,'F':[3]*4}))
        self.assertEqual(r['category'],'reversal_R_up_T_down');self.assertEqual(r['unstructured_log_radius_exclusive'],'1/2')
        self.assertEqual(sum(len(x['joint_cases']) for x in r['maps'].values()),128)

    def test_same_direction(self):
        self.assertEqual(self.reference(pair({'R':[2]*4,'D':[0]*4,'F':[-3]*4}))['category'],'same_direction_a_higher')

    def test_one_preparation(self):
        self.assertFalse(self.reference(pair({'R':[2,2,2,-1],'D':[0]*4,'F':[3]*4}))['category'].startswith('reversal_'))

    def test_one_map(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});x['geo']['a']['assays']['F']['groups'][0]['median']=-1.0
        x['geo']['a']['assays']['F']['groups'][0]['values']=[-1.0]
        self.assertFalse(self.reference(x)['category'].startswith('reversal_'))

    def test_tie(self):
        self.assertEqual(self.reference(pair({'R':[0,2,2,2],'D':[0]*4,'F':[3]*4}))['category'],'unresolved_with_ties')

    def test_missing(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});x['code']['a']['assays']['F']['groups'].pop()
        r=self.reference(x);self.assertEqual(r['category'],'unresolved_missing_or_misaligned');self.assertIsNone(r['unstructured_log_radius_exclusive'])
        self.assertEqual(len(r['maps']['code']['joint_cases']),64)

    def test_mask(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4})
        for name in ('a','b'):
            g=x['code'][name]['assays']['R']['groups'][0];g['samples']=['R1','Rrepeat'];g['values']=[g['median'],None];g['technical_finite']=1
        x['code']['b']['assays']['R']['groups'][0]['values']=[None,0.0]
        self.assertEqual(self.reference(x)['category'],'unresolved_missing_or_misaligned')

    def test_sample_mismatch(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});x['geo']['b']['assays']['D']['groups'][0]['samples']=['different']
        self.assertEqual(self.reference(x)['category'],'unresolved_missing_or_misaligned')

    def test_swap(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});y={m:{'a':v['b'],'b':v['a']} for m,v in x.items()}
        r=self.reference(y);self.assertEqual(r['category'],'reversal_R_down_T_up');self.assertEqual(r['unstructured_log_radius_exclusive'],'1/2')

    def test_common_offsets(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});baseline=self.reference(copy.deepcopy(x))
        for m in x:
            for key in ('a','b'):
                for a in ('R','D','F'):
                    for g in x[m][key]['assays'][a]['groups']:
                        g['median']+=5;g['values']=[v+5 for v in g['values']]
        y=self.reference(x);self.assertEqual(y['combined_bounds'],baseline['combined_bounds'])

    def test_shared_D_falsifier(self):
        r=self.reference(pair({'R':[0]*4,'D':[-1]*4,'F':[0]*4}))
        self.assertEqual(r['category'],'reversal_R_up_T_down');self.assertEqual(r['unstructured_log_radius_exclusive'],'1/4')

    def test_NC_no_cancellation(self):
        def ratio(Z):return (F(3)/Z-F(1,10))/(F(2)/Z-F(1,5))
        self.assertNotEqual(ratio(1),ratio(2));self.assertEqual((F(3)/1)/(F(2)/1),(F(3)/2)/(F(2)/2))

    def test_radius_boundary(self):
        # Four independent allele-assay perturbations move the R difference by4delta.
        margin=F(2);radius=margin/4
        for delta in (F(0),radius/2,radius):
            changed=margin-4*delta
            self.assertEqual(changed>0,delta<radius)
        self.assertEqual(margin-4*radius,0)

    def test_random_distinct_reference(self):
        rng=random.Random(98133)
        for _ in range(200):
            x=pair({a:[rng.randint(-5,5) for _ in range(4)] for a in ('R','D','F')})
            self.reference(x)

    def test_invalid(self):
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});x['code']['a']['assays']['R']['groups'][0]['median']=float('nan')
        with self.assertRaises(ValueError):self.reference(x)
        x=pair({'R':[2]*4,'D':[0]*4,'F':[3]*4});x['code']['a']['assays']['R']['groups'][0]['bio']=True
        with self.assertRaises(ValueError):self.reference(x)

    def test_optimization_refused(self):
        for mode in ('-O','-OO'):
            r=subprocess.run([sys.executable,mode,str(Path(A.__file__))],capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0);self.assertIn('Optimized Python unsupported',r.stderr)

    def test_fresh_directory_guard(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileExistsError):A.execute('nonexistent','nonexistent','nonexistent','nonexistent',temp)

    def test_hash_guard_before_decode(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp)
            for name in ('analyze.py','verify.py','test_synthetic.py'):
                (d/name).write_bytes((Path(__file__).parent/name).read_bytes())
            for name in ('proposed_protocol.md','synthetic_receipt.json'):(d/name).write_text('synthetic only')
            roles={n:d/(n+'.txt') for n in ('input','inventory','contract')}
            for path in roles.values():path.write_text('deliberately invalid JSON')
            hashes={n:A.sha(p) for n,p in roles.items()};hashes['input']='0'*64
            freeze={'status':'FROZEN_SYNTHETIC_PASS_ROOT_AUTHORIZED_EXPLORATORY_RUN','input_hashes':hashes,
                    'owned_files':{p.name:A.sha(p) for p in d.iterdir() if p.suffix in ('.py','.md','.json')}}
            A.dump(d/'freeze.json',freeze)
            with self.assertRaisesRegex(ValueError,'Input hash mismatch'):
                A.execute(d/'freeze.json',roles['input'],roles['inventory'],roles['contract'],d/'out')
            self.assertFalse((d/'out').exists())

if __name__=='__main__':
    if not __debug__:raise SystemExit('Optimized Python unsupported')
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    payload={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
             'independent_formulation_random_fields':200,'biological_inputs_read':False,
             'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(A.__file__),Path(V.__file__))}}
    out=Path(__file__).parent/'synthetic_receipt.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    if not result.wasSuccessful():raise SystemExit(1)
