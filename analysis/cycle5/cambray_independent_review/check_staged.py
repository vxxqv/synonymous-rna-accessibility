"""Synthetic-only independent staged-runner challenge. No observed-input option."""
if not __debug__:
    raise SystemExit('Refuse disabled assertions before all input access')
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest import mock
import zipfile
sys.dont_write_bytecode=True
import reference as ref

HELD='5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5'
LABELS={'joint_connected':'jointly_connected','joint_only':'joint_only_obstruction',
        'individual_disconnected':'single_trait_disconnection'}
MAP={'additive':'additive','neighbor_mean':'adjacent_mean','three_mean':'three_corner_mean'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def close(a,b):
    if a is None or b is None:assert a is b,(a,b)
    else:assert math.isclose(a,b,rel_tol=3e-14,abs_tol=3e-14),(a,b)

def fixture(root,runner,poison=False):
    """Invent twelve disjoint squares, eight groups, four target orientations."""
    units=[];components=[];edges=[];raw={};expected=[]
    for n in range(12):
        target=n%4;ids=[f'ind_{n}_{i}' for i in range(4)];seed=str(10+n%8)
        # Repeated groups carry intentionally different errors.
        if n%3==0:means=[[1.,1.],[0.,2.],[2.,0.],[1.,1.]]
        elif n%3==1:means=[[2.,2.],[1.,1.],[4.,4.],[3.,3.]]
        else:means=[[2.,2.],[4.,1.],[1.,4.],[2.5,2.5]]
        reps={i:[means[i][0]+z for z in (-.125,0,.125,0)]+[means[i][1]+z for z in (-.25,0,.25)] for i in range(4)}
        for i,identifier in enumerate(ids):raw[identifier]=reps[i]+means[i]
        component=ids[:]
        if n==0:
            # Partial protein information never overrides a known growth failure.
            extra='partial_component';component.append(extra)
            raw[extra]=['NA',1.,1.,1.,-100.,-100.,-100.,'NA',-100.]
            edges.extend([(extra,ids[0]),(extra,ids[3])])
        components.append({'IDs':component})
        for a,b in ref.SQUARE_EDGES:edges.append((ids[a],ids[b]))
        units.append({'square_index':n,'seed':seed,'arm':'development',
            'coordinate_order_00_01_10_11':ids,'withheld_ID':ids[target],
            'diagonals':[[ids[0],ids[3]],[ids[1],ids[2]]],
            'component_index':n,'component_size':len(component),'isolated_observed_square':n!=0})
        known={i:means[i] for i in range(4) if i!=target};preds=ref.local_predictions(known,target)
        ep=(target^1,target^2)
        predicted={m:ref.classify({i:(v if i==target else means[i]) for i in range(4)},ref.SQUARE_EDGES,*ep)['status']=='joint_only' for m,v in preds.items()}
        expected.append({'group':seed,'truth':means[target],'predictions':preds,
            'true_class':ref.classify(dict(enumerate(means)),ref.SQUARE_EDGES,*ep)['status']=='joint_only',
            'predicted_class':predicted,'reps':reps,'means':means})
    split={'units':units,'seed_groups':[[str(10+n)] for n in range(8)]}
    masked={u['withheld_ID'] for u in units}
    mask={'arms':{'development':{'masked_IDs':sorted(masked)}},'units':[
        dict(u,prediction_eligible_from_metadata=True) for u in units]}
    dump(root/'split.json',split);dump(root/'mask.json',mask)
    (root/'components.jsonl').write_text(''.join(json.dumps(c)+'\n' for c in components))
    (root/'edges.tsv').write_text('left_ID\tright_ID\n'+''.join(a+'\t'+b+'\n' for a,b in edges))
    (root/'protocol.md').write_text('Independent synthetic fixture only.\n')
    header=['ID','ignored']+runner.REPLICATES+runner.AUTHOR_MEANS+['other']
    lines=[','.join('"'+s+'"' for s in header).encode()+b'\n']
    for identifier,values in raw.items():
        fields=[str(v).encode() for v in values]
        if poison:
            if identifier in masked or identifier=='partial_component':fields=[b'\xffSEALED']*9
            else:fields[-2:]=[b'\xffAUTHOR']*2
        lines.append(b'"'+identifier.encode()+b'",\xff,'+b','.join(fields)+b',\xff\n')
    # Other-arm/unrelated numeric bytes cannot even be decoded as text.
    lines.append(b'"unrelated",\xff,"unclosed protected numeric bytes\n')
    with zipfile.ZipFile(root/'toy.zip','w') as z:z.writestr('invented.csv',b''.join(lines))
    paths={'protocol':root/'protocol.md','split':root/'split.json','mask':root/'mask.json',
        'components':root/'components.jsonl','edges':root/'edges.tsv','archive':root/'toy.zip',
        'runner':Path(runner.__file__),'graph_rules':Path(runner.rules.__file__),'masking':Path(runner.masking.__file__)}
    auth={'schema':'Cambray_staged_v1','execution_authorized':True,'allowed_phases':['predict','evaluate'],
        'allowed_arms':['development'],'expected_records':len(raw)+1,
        'inputs':{r:{'path':str(p.resolve()),'sha256':sha(p)} for r,p in paths.items()}}
    dump(root/'auth.json',auth)
    return auth,units,expected

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--primary',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('Fresh output required')
    assert sha(args.primary)==HELD
    before={str(p):sha(p) for p in (args.primary,Path(__file__),Path(ref.__file__))}
    sys.path.insert(0,str(args.primary.resolve().parent))
    spec=importlib.util.spec_from_file_location('held_staged',args.primary);runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    checks=[]
    with tempfile.TemporaryDirectory(prefix='independent_stage_',dir=Path(__file__).parent) as tmp:
        root=Path(tmp);clean=root/'clean';clean.mkdir();auth,units,expected=fixture(clean,runner)
        runner.run('predict',clean/'auth.json',clean,'development',clean/'predict')
        pred=json.loads((clean/'predict/predictions.json').read_text())
        for p,e,u in zip(pred['records'],expected,units):
            for independent,primary in MAP.items():
                for a,b in zip(e['predictions'][independent],p['predictions'][primary]):close(a,b)
                assert e['predicted_class'][independent]==p['classifications'][primary]['joint_only_obstruction']
            target=u['coordinate_order_00_01_10_11'].index(u['withheld_ID'])
            assert p['classification_endpoint_IDs']==[u['coordinate_order_00_01_10_11'][target^i] for i in (1,2)]
        checks.append('12 independent local predictions, all 4 target orientations and classification endpoints')
        seal=clean/'predict/seal.json';auth['inputs']['prediction_seal']={'path':str(seal),'sha256':sha(seal)};dump(clean/'eval_auth.json',auth)
        runner.run('evaluate',clean/'eval_auth.json',clean,'development',clean/'evaluate')
        observed=json.loads((clean/'evaluate/evaluation.json').read_text())
        expected_metrics=ref.grouped_metrics(expected);primary_metrics=observed['prediction_evaluation']
        assert primary_metrics['scored_squares']==12 and primary_metrics['related_seed_groups']==8
        assert primary_metrics['discrimination_screen_passed']==expected_metrics['success']
        for m,p in MAP.items():
            e=expected_metrics['metrics'][m];o=primary_metrics['summaries'][p]
            for i,t in enumerate(('protein','growth')):
                close(e['MAE'][i],o['continuous'][t]['MAE']);close(e['RMSE'][i],o['continuous'][t]['RMSE'])
            close(e['misclassification'],o['group_balanced_classification_error']['value'])
            assert e['confusion']==o['pooled_confusion'];close(e['balanced_accuracy'],o['pooled_balanced_accuracy'])
        checks.append('all 3 methods, both traits: equal-group MAE/RMSE, confusion, balanced accuracy and full screen')
        for description,e in zip(observed['descriptions'],expected):
            for d,ep in zip(description['diagonals'],((0,3),(1,2))):
                r=ref.preparation_summary(e['reps'],s=ep[0],t=ep[1])
                assert d['primary']['classification']==LABELS[r['point']['status']]
                assert d['preparation_consistent_joint_obstruction']==r['all12_joint_only']
                assert [s['local']['classification'] for s in d['preparation_scenarios']]==[LABELS[x] for x in r['preparation_statuses']]
                assert (d['local_range_certificate']['status']=='certified_joint_only_obstruction')==r['range_certificate']['joint_only_certificate']
                for k in range(2):
                    values=[v[k] for v in e['means']]
                    if description['square_index']==0 and k==1:values.append(-100.)
                    close(d['point_tolerance'][k],1e-10*max(1,max(map(abs,values))))
        checks.append('24 diagonal point/range and 288 whole-preparation flags, component-derived point tolerances')
        partial=observed['descriptions'][0]['diagonals'][0]['observed_component']
        assert not partial['complete'] and partial['point'] is None
        checks.append('partial component preserves trait-specific missingness without dropping complete square')
        poison=root/'poison';poison.mkdir();fixture(poison,runner,True)
        runner.run('predict',poison/'auth.json',poison,'development',poison/'predict')
        poisonpred=json.loads((poison/'predict/predictions.json').read_text())
        assert pred['records']==poisonpred['records']
        checks.append('invalid UTF-8 targets, author means, component-only vertex and unrelated row not decoded before seal')
        # All tampering refusals must precede the phenotype reader.
        for payload in ('predictions.json','known_inputs.json'):
            p=clean/'predict'/payload;original=p.read_bytes();p.write_bytes(original+b' ')
            with mock.patch.object(runner,'selective_reader',side_effect=AssertionError('premature reader')):
                try:runner.run('evaluate',clean/'eval_auth.json',clean,'development',clean/('bad_'+payload))
                except ValueError:pass
                else:raise AssertionError('tampered payload accepted')
            p.write_bytes(original)
        checks.append('both sealed numerical payloads refuse tampering before any phenotype-reader call')
        for flag in ('-O','-OO'):
            out=root/('forbidden'+flag[1:]);p=subprocess.run([sys.executable,flag,str(args.primary.resolve()),'predict','--authorization',str(root/'absent'),'--root',str(root),'--arm','development','--output',str(out)],capture_output=True,text=True)
            assert p.returncode!=0 and 'Disabled assertions' in p.stderr and not out.exists()
        checks.append('-O and -OO fail before missing authorization read or output creation')
        try:runner.run('predict',clean/'auth.json',clean,'development',clean/'predict')
        except FileExistsError:pass
        else:raise AssertionError('overwrite accepted')
        checks.append('normal mode fresh-output-only refusal')
    after={p:sha(p) for p in before};assert before==after
    args.out.mkdir(parents=True,exist_ok=False)
    dump(args.out/'receipt.json',{'status':'passed','checks':checks,'input_hashes_before_after_equal':before,
        'real_phenotype_sources_opened':0,'independent_reference_imports_primary':False,
        'fixture_scope':'12 invented squares, 8 repeated groups, 4 target orientations, one partially measured component-only vertex',
        'independence_limit':'The staged runner is executed for comparison; all expected numerical predictions, classes and grouped metrics are from the frozen reference.'})
    print(json.dumps({'status':'passed','checks':len(checks),'receipt':str(args.out/'receipt.json')}))

if __name__=='__main__':main()
