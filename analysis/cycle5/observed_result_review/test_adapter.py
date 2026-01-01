"""Invented-data tests only. Never discovers or opens biological result paths."""
if not __debug__:raise SystemExit('Disabled assertions refused')
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile
sys.dont_write_bytecode=True
import verify_observed as adapter

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PRIMARY=BASE/'cambray_pilot_preflight/staged_analysis.py'
REFERENCE=BASE/'cambray_independent_review/reference.py'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

class AdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert adapter.sha(PRIMARY)=='5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5'
        assert adapter.sha(REFERENCE)==adapter.REFERENCE_SHA
        sys.path.insert(0,str(PRIMARY.parent));sys.path.insert(0,str(REFERENCE.parent))
        cls.primary=module('synthetic_primary_staged',PRIMARY)
        cls.fixture_module=module('existing_independent_fixture',REFERENCE.parent/'check_staged.py')
    def temp(self):return tempfile.TemporaryDirectory(prefix='adapter_invented_',dir=HERE)
    def stage(self,root,mutation=None,constant=False,overlap=False):
        auth,units,_=self.fixture_module.fixture(root,self.primary)
        if overlap:
            # Two genuine square faces share an edge and both original targets.
            # They remain descriptive, but each contains two union-mask IDs.
            split=json.loads((root/'split.json').read_text());old=split['units'][1]['coordinate_order_00_01_10_11'][:]
            first=split['units'][0]['coordinate_order_00_01_10_11'];second=[first[0],first[1],old[2],old[3]]
            split['units'][1]['coordinate_order_00_01_10_11']=second;split['units'][1]['withheld_ID']=first[1]
            split['units'][1]['diagonals']=[[second[0],second[3]],[second[1],second[2]]]
            split['units'][1]['component_index']=0
            split['units'][0]['isolated_observed_square']=split['units'][1]['isolated_observed_square']=False
            for unit in split['units'][2:]:unit['component_index']-=1
            components=[json.loads(line) for line in (root/'components.jsonl').read_text().splitlines()]
            components[0]['IDs']+=old[2:];components.pop(1)
            (root/'components.jsonl').write_text(''.join(json.dumps(c)+'\n' for c in components))
            existing=[line.split('\t') for line in (root/'edges.tsv').read_text().splitlines()[1:]]
            edges={tuple(sorted(e)) for e in existing if not (set(e)<=set(old))}
            for a,b in ((0,1),(0,2),(1,3),(2,3)):edges.add(tuple(sorted((second[a],second[b]))))
            (root/'edges.tsv').write_text('left_ID\tright_ID\n'+''.join(a+'\t'+b+'\n' for a,b in sorted(edges)))
            (root/'split.json').write_text(json.dumps(split));fixed={u['withheld_ID'] for u in split['units']}
            mask={'arms':{'development':{'masked_IDs':sorted(fixed)}},'units':[dict(u,prediction_eligible_from_metadata=len(set(u['coordinate_order_00_01_10_11'])&fixed)==1) for u in split['units']]}
            (root/'mask.json').write_text(json.dumps(mask))
            for role in ('split','mask','components','edges'):auth['inputs'][role]['sha256']=adapter.sha(Path(auth['inputs'][role]['path']))
            units=split['units']
        if mutation or constant:
            with zipfile.ZipFile(root/'toy.zip') as z:lines=z.read('invented.csv').splitlines(keepends=True)
            if constant:
                for n,line in enumerate(lines[1:-1],1):
                    cells=line.rstrip(b'\n').split(b',');cells[2:11]=[b'1']*9;lines[n]=b','.join(cells)+b'\n'
            if mutation:
                target=units[0]['withheld_ID'];known=next(i for i in units[0]['coordinate_order_00_01_10_11'] if i!=target)
                identifier=target if mutation in ('target','all_missing') else known
                for n,line in enumerate(lines[1:-1],1):
                    cells=line.rstrip(b'\n').split(b',')
                    if mutation=='all_missing' or cells[0]==b'"'+identifier.encode()+b'"':
                        cells[9 if mutation=='mean_hold' else 2]=b'9.5' if mutation=='mean_hold' else b'NA'
                        lines[n]=b','.join(cells)+b'\n'
            with zipfile.ZipFile(root/'toy.zip','w') as z:z.writestr('invented.csv',b''.join(lines))
            auth['inputs']['archive']['sha256']=adapter.sha(root/'toy.zip');adapter.write(root/'auth_changed.json',auth);authpath=root/'auth_changed.json'
        elif overlap:adapter.write(root/'auth_overlap.json',auth);authpath=root/'auth_overlap.json'
        else:authpath=root/'auth.json'
        self.primary.run('predict',authpath,root,'development',root/'prediction')
        auth['inputs']['prediction_seal']={'path':str(root/'prediction/seal.json'),'sha256':adapter.sha(root/'prediction/seal.json')}
        adapter.write(root/'evaluate_auth.json',auth)
        self.primary.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
        paths={'evaluation':root/'evaluation/evaluation.json','unsealed_inputs':root/'evaluation/unsealed_inputs.json',
            'predictions':root/'prediction/predictions.json','prediction_seal':root/'prediction/seal.json','known_inputs':root/'prediction/known_inputs.json',
            'split':root/'split.json','components':root/'components.jsonl','edges':root/'edges.tsv',
            'adapter':Path(adapter.__file__),'reference':REFERENCE}
        independent={'schema':'Cambray_independent_result_v1','execution_authorized':True,'data_scope':'synthetic','arm':'development',
            'inputs':{r:{'path':str(p.resolve()),'sha256':adapter.sha(p)} for r,p in paths.items()}}
        adapter.write(root/'independent_auth.json',independent)
        return independent,paths
    def compare(self,root):return adapter.run(root/'independent_auth.json',root,root/'comparison')
    def test_complete_all_core_quantities(self):
        with self.temp() as tmp:
            root=Path(tmp);self.stage(root);r=self.compare(root)
            self.assertEqual(r['status'],'verified_complete_arm');self.assertEqual(r['diagonals_verified'],24)
            self.assertEqual(r['preparation_scenarios_verified'],288);self.assertGreater(r['exact_critical_cells_verified'],0)
            self.assertEqual(r['prediction_cases'],12);self.assertFalse(r['primary_implementation_imported'])
    def test_missing_target_does_not_erase_sealed_prediction(self):
        with self.temp() as tmp:
            root=Path(tmp);_,paths=self.stage(root,'target');r=self.compare(root)
            self.assertEqual(r['status'],'verified_complete_arm');self.assertEqual(r['prediction_cases'],11)
            self.assertEqual(json.loads(paths['predictions'].read_text())['records'][0]['status'],'predicted')
    def test_missing_known_retained(self):
        with self.temp() as tmp:
            root=Path(tmp);self.stage(root,'known');r=self.compare(root)
            self.assertEqual(r['prediction_cases'],11)
    def test_all_missing_empty_metrics_and_coverage(self):
        with self.temp() as tmp:
            root=Path(tmp);self.stage(root,'all_missing');r=self.compare(root)
            self.assertEqual(r['prediction_cases'],0);self.assertEqual(r['diagonals_verified'],0)
            self.assertFalse(r['prediction_screen_passed'])
    def test_one_class_and_zero_baselines(self):
        with self.temp() as tmp:
            root=Path(tmp);self.stage(root,constant=True);r=self.compare(root)
            self.assertFalse(r['prediction_screen_passed'])
    def test_provenance_hold_only(self):
        with self.temp() as tmp:
            root=Path(tmp);self.stage(root,'mean_hold');r=self.compare(root)
            self.assertEqual(r['status'],'verified_source_hold');self.assertGreater(len(r['holds']),0)
            self.assertNotIn('primary_summary',r)
    def test_shared_mask_exclusions_do_not_drop_descriptions(self):
        with self.temp() as tmp:
            root=Path(tmp);_,paths=self.stage(root,overlap=True);r=self.compare(root)
            self.assertEqual(r['prediction_cases'],10);self.assertEqual(r['diagonals_verified'],24)
            records=json.loads(paths['predictions'].read_text())['records']
            self.assertEqual([x['status'] for x in records[:2]],['excluded_shared_corner_mask']*2)
    def tamper(self,root,transform):
        auth,paths=self.stage(root);evaluation=json.loads(paths['evaluation'].read_text());transform(evaluation)
        paths['evaluation'].write_text(json.dumps(evaluation));auth['inputs']['evaluation']['sha256']=adapter.sha(paths['evaluation'])
        (root/'independent_auth.json').write_text(json.dumps(auth))
        with self.assertRaises(ValueError):self.compare(root)
        self.assertTrue((root/'comparison/failure.json').exists())
    def test_tampered_label_rejected(self):
        with self.temp() as tmp:
            self.tamper(Path(tmp),lambda e:e['descriptions'][0]['diagonals'][0]['primary'].update(classification='wrong'))
    def test_tampered_exact_geometry_rejected(self):
        with self.temp() as tmp:
            self.tamper(Path(tmp),lambda e:e['descriptions'][0]['diagonals'][0]['threshold_geometry_tolerance_zero']['cells'][0].update(classification='wrong'))
    def test_tampered_group_metric_rejected(self):
        with self.temp() as tmp:
            self.tamper(Path(tmp),lambda e:e['prediction_evaluation']['summaries']['additive']['continuous']['protein'].update(group_MAE={'fake':0}))
    def test_wrong_arm_before_selected_exports(self):
        with self.temp() as tmp:
            root=Path(tmp);auth,_=self.stage(root);auth['arm']='reserved_validation';(root/'independent_auth.json').write_text(json.dumps(auth))
            original=adapter.load
            def guarded(path):
                if Path(path).name in ('evaluation.json','unsealed_inputs.json'):raise AssertionError('Wrong arm reached selected outcome loader')
                return original(path)
            with mock.patch.object(adapter,'load',side_effect=guarded):
                with self.assertRaises(ValueError):self.compare(root)
    def test_hash_failure_before_selected_outcomes_or_output(self):
        with self.temp() as tmp:
            root=Path(tmp);auth,_=self.stage(root);auth['inputs']['unsealed_inputs']['sha256']='0'*64;(root/'independent_auth.json').write_text(json.dumps(auth))
            original=adapter.load
            def guarded(path):
                if Path(path).name!='independent_auth.json':raise AssertionError('Hash failure reached selected outcome loader')
                return original(path)
            with mock.patch.object(adapter,'load',side_effect=guarded):
                with self.assertRaises(ValueError):self.compare(root)
            self.assertFalse((root/'comparison').exists())
    def test_sealed_payload_tamper_detected(self):
        with self.temp() as tmp:
            root=Path(tmp);auth,paths=self.stage(root);paths['predictions'].write_bytes(paths['predictions'].read_bytes()+b' ')
            auth['inputs']['predictions']['sha256']=adapter.sha(paths['predictions']);(root/'independent_auth.json').write_text(json.dumps(auth))
            with self.assertRaises(ValueError):self.compare(root)
    def test_extra_input_role_duplicate_json_and_overwrite(self):
        with self.temp() as tmp:
            root=Path(tmp);auth,_=self.stage(root);auth['inputs']['archive']={'path':'must_not_open','sha256':'0'*64};(root/'independent_auth.json').write_text(json.dumps(auth))
            with self.assertRaises(ValueError):self.compare(root)
            (root/'duplicate.json').write_text('{"x":1,"x":2}')
            with self.assertRaises(ValueError):adapter.load(root/'duplicate.json')
            (root/'exists').mkdir()
            with mock.patch.object(adapter,'authorize',side_effect=AssertionError('Authorization accessed')):
                with self.assertRaises(FileExistsError):adapter.run(root/'absent',root,root/'exists')
    def test_optimized_execution_refused_before_io(self):
        with self.temp() as tmp:
            root=Path(tmp)
            for flag in ('-O','-OO'):
                out=root/flag[1:]
                p=subprocess.run([sys.executable,flag,str(Path(adapter.__file__)), '--authorization',str(root/'absent'),'--root',str(root),'--out',str(out)],capture_output=True,text=True)
                self.assertNotEqual(p.returncode,0);self.assertIn('Disabled assertions',p.stderr);self.assertFalse(out.exists())

if __name__=='__main__':
    before={str(p):adapter.sha(p) for p in (PRIMARY,REFERENCE,Path(adapter.__file__),Path(__file__))}
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AdapterTests))
    unchanged=all(adapter.sha(p)==h for p,h in before.items())
    receipt={'status':'passed' if r.wasSuccessful() and unchanged else 'failed','tests':r.testsRun,'failures':len(r.failures),'errors':len(r.errors),
        'input_hashes_before_after_equal':before if unchanged else False,'real_biological_sources_opened':0,
        'fixture':'Twelve invented squares, eight groups, one partially measured component-only vertex; no archive/outcome path discovery.'}
    with (HERE/'synthetic_receipt.json').open('x',encoding='utf-8') as stream:json.dump(receipt,stream,indent=2);stream.write('\n')
    print(json.dumps(receipt));raise SystemExit(0 if receipt['status']=='passed' else 1)
