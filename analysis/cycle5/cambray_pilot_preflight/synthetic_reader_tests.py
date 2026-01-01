"""Synthetic-only staged-reader, sealing, scoring and refusal checks.

Run with the existing Python interpreter. Temporary fixtures contain invented
identifiers and values, never the biological source archive. All temporary files
are removed after each test. The final result is printed, not written over prior
preflight receipts.
"""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile
import staged_analysis as runner
import graph_rules as rules
import masking_v2 as masking

HERE=Path(__file__).resolve().parent


def fixture(root,poison_hidden=False,missing_target=False,missing_known=False,author_delta=0,constant=False):
    units=[];components=[];edges=[];records={};groups=[]
    for k,arm in enumerate(('development','development','reserved_validation')):
        seed=str(k+2);ids=[f'toy_{k}_{i}' for i in range(4)];groups.append([seed])
        target=sorted(ids)[rules.corner_index(sorted(ids))];m=ids.index(target)
        fields=[[1.,1.] for i in range(4)]
        if k==0 and not constant:fields[m]=[0.,2.];fields[m^3]=[2.,0.]
        units.append({'square_index':k,'seed':seed,'arm':arm,'coordinate_order_00_01_10_11':ids,'withheld_ID':target,
            'diagonals':[[ids[0],ids[3]],[ids[1],ids[2]]],'component_index':k,'component_size':4,'isolated_observed_square':True})
        components.append({'IDs':ids})
        for a in range(4):
            for b in (a^1,a^2):
                if a<b:edges.append((ids[a],ids[b]))
        for i,identifier in enumerate(ids):
            values=[fields[i][0]]*4+[fields[i][1]]*3
            means=[fields[i][0],fields[i][1]]
            if author_delta and k==1:means[0]+=author_delta
            if missing_target and k==0 and identifier==target:values[0]='NA'
            if missing_known and k==0 and identifier==ids[m^1]:values[0]='NA'
            if poison_hidden and (arm=='reserved_validation' or identifier==target):values=[b'\xffinvalid']*7
            if poison_hidden:means=[b'\xffauthor']*2
            records[identifier]=values+means
    split={'units':units,'seed_groups':groups}
    for name,obj in [('split.json',split),('mask.json',masking.union_mask_manifest(split))]:
        (root/name).write_text(json.dumps(obj),encoding='utf-8')
    (root/'components.jsonl').write_text('\n'.join(json.dumps(c) for c in components)+'\n')
    (root/'edges.tsv').write_text('left_ID\tright_ID\tposition_1based\n'+''.join(f'{a}\t{b}\t1\n' for a,b in edges))
    (root/'protocol.md').write_text('Frozen synthetic-only protocol fixture.\n')
    header=['','unused_before']+runner.REPLICATES+runner.AUTHOR_MEANS+['unused_after']
    lines=[','.join('"'+c+'"' for c in header).encode()+b'\n']
    for identifier,values in records.items():
        cells=[b'"'+identifier.encode()+b'"',b'\xffignored']+[v if isinstance(v,bytes) else str(v).encode() for v in values]+[b'\xffignored']
        lines.append(b','.join(cells)+b'\n')
    lines.append(b'"unrelated",\xff\xfe,"malformed unreachable remainder\n')
    with zipfile.ZipFile(root/'data.zip','w') as archive:archive.writestr('data.csv',b''.join(lines))
    paths={'protocol':root/'protocol.md','split':root/'split.json','mask':root/'mask.json','components':root/'components.jsonl',
        'edges':root/'edges.tsv','archive':root/'data.zip','graph_rules':Path(rules.__file__),'masking':Path(masking.__file__),'runner':Path(runner.__file__)}
    auth={'schema':'Cambray_staged_v1','execution_authorized':True,'allowed_phases':['predict','evaluate'],
        'allowed_arms':['development','reserved_validation'],'expected_records':13,
        'inputs':{role:{'path':str(p.resolve()),'sha256':runner.sha(p)} for role,p in paths.items()}}
    (root/'auth.json').write_text(json.dumps(auth))
    return units,auth


def add_seal(root,auth,directory='prediction'):
    auth=copy.deepcopy(auth);p=root/directory/'seal.json'
    auth['inputs']['prediction_seal']={'path':str(p),'sha256':runner.sha(p)}
    (root/'evaluate_auth.json').write_text(json.dumps(auth));return auth


class ReaderTests(unittest.TestCase):
    def temp(self):return tempfile.TemporaryDirectory(prefix='c5_synthetic_',dir=HERE)

    def test_hidden_other_arm_and_author_mean_poison_never_decoded(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root,poison_hidden=True)
            result=runner.run('predict',root/'auth.json',root,'development',root/'prediction')
            self.assertEqual(result['status'],'predictions_sealed_before_target_access')
            payload=json.loads((root/'prediction/known_inputs.json').read_text())
            fixed={u['withheld_ID'] for u in units if u['arm']=='development'}
            self.assertFalse(set(payload['rows'])&fixed)
            self.assertEqual(payload['reader']['decoded_fields'],runner.REPLICATES)
            self.assertEqual(len(payload['rows']),6)
            self.assertTrue(all(set(r)==set(runner.REPLICATES) for r in payload['rows'].values()))

    def test_end_to_end_controls_exact_predictions_and_descriptions(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction')
            prediction=json.loads((root/'prediction/predictions.json').read_text())
            target0=next(r for r in prediction['records'] if r['square_index']==0)
            self.assertEqual(target0['predictions']['additive'],[0,2])
            self.assertEqual(target0['predictions']['adjacent_mean'],[1,1])
            self.assertTrue(target0['classifications']['additive']['joint_only_obstruction'])
            self.assertFalse(target0['classifications']['adjacent_mean']['joint_only_obstruction'])
            add_seal(root,auth)
            runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
            evaluated=json.loads((root/'evaluation/evaluation.json').read_text())
            self.assertEqual(evaluated['primary_summary']['eligible_squares'],2)
            self.assertEqual(evaluated['primary_summary']['eligible_diagonals'],4)
            self.assertEqual(evaluated['primary_summary']['joint_only_obstructions'],1)
            self.assertEqual(evaluated['primary_summary']['preparation_consistent_diagonals'],1)
            self.assertEqual(evaluated['primary_summary']['range_certified_diagonals'],1)
            self.assertEqual(evaluated['summary_by_observed_component_stratum']['isolated_observed_square']['selected_squares'],2)
            self.assertEqual(evaluated['summary_by_observed_component_stratum']['larger_observed_component']['selected_squares'],0)
            self.assertFalse(evaluated['coverage_passed'])
            scores=evaluated['prediction_evaluation']
            self.assertTrue(scores['both_true_classes'])
            self.assertFalse(scores['discrimination_screen_passed'])
            for trait in runner.TRAITS:
                self.assertEqual(scores['summaries']['additive']['continuous'][trait]['MAE'],0)
                self.assertEqual(scores['summaries']['adjacent_mean']['continuous'][trait]['MAE'],.5)
            self.assertEqual(scores['summaries']['additive']['pooled_confusion'],{'TP':1,'TN':1,'FP':0,'FN':0})
            self.assertEqual(scores['summaries']['additive']['pooled_balanced_accuracy'],1)
            for square in evaluated['descriptions']:
                self.assertEqual(len(square['diagonals']),2)
                for diagonal in square['diagonals']:
                    self.assertEqual(len(diagonal['preparation_scenarios']),12)
                    self.assertEqual(diagonal['observed_component']['point']['classification'],diagonal['primary']['classification'])

    def test_target_missingness_only_discovered_after_seal(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root,missing_target=True)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction')
            first=json.loads((root/'prediction/predictions.json').read_text())['records'][0]
            self.assertEqual(first['status'],'predicted')
            self.assertNotIn('target_finite',first)
            add_seal(root,auth)
            runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
            result=json.loads((root/'evaluation/evaluation.json').read_text())
            self.assertEqual(len(result['coverage']),2)
            self.assertFalse(result['coverage'][0]['eligible'])
            self.assertEqual(result['prediction_evaluation']['records'][0]['status'],'unscorable_missing_target_or_square')

    def test_known_missingness_retained_without_target_access(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root,missing_known=True,poison_hidden=True)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction')
            result=json.loads((root/'prediction/predictions.json').read_text())
            self.assertEqual(result['records'][0]['status'],'missing_known_replicates')

    def test_mean_rounding_uses_max_tolerance_and_holds(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root,author_delta=1.005e-8)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction')
            add_seal(root,auth)
            result=runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
            self.assertEqual(result['status'],'held_for_author_mean_provenance')
            payload=json.loads((root/'evaluation/evaluation.json').read_text())
            self.assertNotIn('descriptions',payload)
            self.assertTrue(any(m['status']=='mean_discrepancy' for m in payload['mean_diagnostics']))

    def test_zero_baseline_and_one_class_fail_without_relative_division(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root,constant=True)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction');add_seal(root,auth)
            runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
            result=json.loads((root/'evaluation/evaluation.json').read_text())['prediction_evaluation']
            self.assertFalse(result['both_true_classes']);self.assertFalse(result['discrimination_screen_passed'])
            self.assertTrue(all(c['relative_gain'] is None and not c['passes_10_percent'] for c in result['MAE_comparisons']))
            self.assertIsNone(result['summaries']['additive']['pooled_balanced_accuracy'])

    def test_missing_author_mean_is_provenance_hold(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root)
            with zipfile.ZipFile(root/'data.zip') as archive:content=archive.read('data.csv')
            lines=content.splitlines(keepends=True)
            cells=lines[1].rstrip(b'\n').split(b',');cells[9]=b'NA';lines[1]=b','.join(cells)+b'\n'
            with zipfile.ZipFile(root/'data.zip','w') as archive:archive.writestr('data.csv',b''.join(lines))
            auth['inputs']['archive']['sha256']=runner.sha(root/'data.zip')
            (root/'auth.json').write_text(json.dumps(auth))
            runner.run('predict',root/'auth.json',root,'development',root/'prediction');add_seal(root,auth)
            result=runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'evaluation')
            self.assertEqual(result['status'],'held_for_author_mean_provenance')
            payload=json.loads((root/'evaluation/evaluation.json').read_text())
            self.assertTrue(any(m['status']=='missing_author_mean' for m in payload['mean_diagnostics']))
            self.assertNotIn('prediction_evaluation',payload)

    def test_partial_component_trait_preserves_known_failure(self):
        # Unknown protein cannot override the known growth failure on a chain.
        graph=[[1],[0,2],[1]];field=[[1,1],[None,0],[1,1]]
        result=runner.field_check(graph,field,(0,2),0)
        self.assertFalse(result['complete']);self.assertIsNone(result['point'])
        self.assertFalse(result['bounds']['joint_optimistic'])
        row={c:{'value':1,'missing_reason':None,'raw':'1'} for c in runner.REPLICATES}
        row[runner.PROTEIN[0]]={'value':None,'missing_reason':'NA','raw':'NA'}
        self.assertEqual(runner.point(row),[None,1])
        self.assertEqual(runner.ranges(row),[None,[1,1]])

    def test_seal_and_hash_refusal_precede_reader(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root)
            runner.run('predict',root/'auth.json',root,'development',root/'prediction');add_seal(root,auth)
            prediction=root/'prediction/predictions.json'
            prediction.write_bytes(prediction.read_bytes()+b' ')
            with mock.patch.object(runner,'selective_reader',side_effect=AssertionError('Reader called before seal check')):
                with self.assertRaises(ValueError):runner.run('evaluate',root/'evaluate_auth.json',root,'development',root/'bad_seal')
            broken=copy.deepcopy(auth);broken['inputs']['protocol']['sha256']='0'*64
            (root/'broken.json').write_text(json.dumps(broken))
            with mock.patch.object(runner,'selective_reader',side_effect=AssertionError('Reader called before hash check')):
                with self.assertRaises(ValueError):runner.run('predict',root/'broken.json',root,'development',root/'bad_hash')

    def test_no_overwrite_and_validation_decision_gate(self):
        with self.temp() as folder:
            root=Path(folder);units,auth=fixture(root)
            (root/'exists').mkdir()
            with mock.patch.object(runner,'authorize',side_effect=AssertionError('Auth read despite existing output')):
                with self.assertRaises(FileExistsError):runner.run('predict',root/'auth.json',root,'development',root/'exists')
            with mock.patch.object(runner,'selective_reader',side_effect=AssertionError('Validation reader called without decision')):
                with self.assertRaises(KeyError):runner.run('predict',root/'auth.json',root,'reserved_validation',root/'blocked_validation')
            report=root/'development_report.json';report.write_text('{"synthetic_only":true}')
            decision={'status':'development_report_frozen','unchanged_validation_procedure':True,
                'protocol_sha256':auth['inputs']['protocol']['sha256'],'report_sha256':runner.sha(report),
                'code_hashes':{r:auth['inputs'][r]['sha256'] for r in ('graph_rules','masking','runner')}}
            dp=root/'decision.json';dp.write_text(json.dumps(decision))
            for role,path in [('development_decision',dp),('development_report',report)]:auth['inputs'][role]={'path':str(path),'sha256':runner.sha(path)}
            (root/'validation_auth.json').write_text(json.dumps(auth))
            result=runner.run('predict',root/'validation_auth.json',root,'reserved_validation',root/'validation')
            self.assertEqual(result['prediction_records'],1)

    def test_selective_tokenizer_quoted_commas_and_nonfinite(self):
        line=b'"id","ignored,\xff",12.5,"NA",\xffsuffix\n'
        fields=runner.selected_fields(line,[2,3])
        self.assertEqual(runner.parse_numeric(fields[2])['value'],12.5)
        self.assertIsNone(runner.parse_numeric(fields[3])['value'])
        self.assertIsNone(runner.parse_numeric(b'Inf')['value'])
        self.assertIsNone(runner.parse_numeric(b'nan')['value'])
        with self.assertRaises(ValueError):runner.parse_numeric(b'notnumeric')

    def test_group_equal_weighting(self):
        records=[{'group':'A','error':0},{'group':'A','error':0},{'group':'B','error':6}]
        self.assertEqual(runner.grouped_mean(records,lambda r:r['error'])['value'],3)

    def test_disabled_assertions_refused_before_input_or_output(self):
        with self.temp() as folder:
            root=Path(folder)
            for flag in ('-O','-OO'):
                completed=subprocess.run([sys.executable,flag,str(Path(runner.__file__)),
                    'predict','--authorization',str(root/'absent.json'),'--root',str(root),
                    '--arm','development','--output',str(root/'must_not_exist')],capture_output=True,text=True)
                self.assertNotEqual(completed.returncode,0)
                self.assertIn('Disabled assertions',completed.stderr)
                self.assertFalse((root/'must_not_exist').exists())


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReaderTests))
    print(json.dumps({'status':'passed' if result.wasSuccessful() else 'failed','tests_run':result.testsRun,
        'failures':len(result.failures),'errors':len(result.errors),'real_source_ZIP_opened':False,
        'hashes':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('staged_analysis.py','synthetic_reader_tests.py')}}))
    raise SystemExit(0 if result.wasSuccessful() else 1)
