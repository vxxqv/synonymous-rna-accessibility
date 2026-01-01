"""Independent synthetic-only review at explicitly held input hashes."""
import collections
import copy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest import mock

ROOT=Path(__file__).resolve().parent
CYCLE=ROOT.parent
PINS={
    'diagnostics.py':(CYCLE/'amendment_diagnostics/diagnostics.py','9e0307ad72b6dd0be771dee9d886bb87e63d09ed85a7ceea5409b2e3f90f897d'),
    'test_diagnostics.py':(CYCLE/'amendment_diagnostics/test_diagnostics.py','07c54d26dcae20a6f2a436e04fb4656452722e4cb19a46bc92346bd341baf5c8'),
    'amendment':(CYCLE/'cambray_pilot/amendment_01_draft.md','92c5028dff57df644ad18479f79b4be2e604b7bbe7798e30cdd3dfe5a7fa94d0'),
    'staged_schema':(CYCLE/'cambray_pilot_preflight/staged_analysis.py','5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5')}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(ok,message):
    if not ok:raise RuntimeError(message)
def write(p,value):
    with p.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def near(x,y):return math.isclose(x,float(y),rel_tol=1e-12,abs_tol=1e-12)

def reference_loss(records,groups,method,k,target_balanced,squared):
    # Construct explicit rational per-case weights, rather than invoke nested means.
    selected=[r for r in records if r['status']=='scored' and r['group'] in groups]
    active={r['group'] for r in selected};total=F(0)
    if not active:return None
    for r in selected:
        same_group=[z for z in selected if z['group']==r['group']]
        if target_balanced:
            targets={z['withheld_ID'] for z in same_group}
            denominator=len(active)*len(targets)*sum(z['withheld_ID']==r['withheld_ID'] for z in same_group)
        else:denominator=len(active)*len(same_group)
        err=F(r['methods'][method]['predicted_target'][k])-F(r['observed_target'][k])
        total+=(err*err if squared else abs(err))/denominator
    return total

def summarize_components(descriptions):
    result={key:collections.Counter() for key in ('point_component_bounds','point_component_classification',
        'point_bypass','alternative_point_bypass','range_component_bounds','range_local_bounds')}
    scenarios={}
    def bypass(local,component):
        if local['classification']!='joint_only_obstruction':return 'not_a_local_joint_only_obstruction'
        b=component['bounds']
        if b['joint_pessimistic']:return 'certified_observed_component_bypass'
        if not b['joint_optimistic']:return 'observed_component_optimistic_disconnection'
        return 'observed_component_bypass_unresolved'
    for square in descriptions:
        for diagonal in square['diagonals']:
            c=diagonal['observed_component']
            result['point_component_bounds'][c['bounds']['status']]+=1
            result['point_component_classification'][c['point']['classification'] if c['point'] else 'point_unavailable']+=1
            result['point_bypass'][bypass(diagonal['primary'],c)]+=1
            result['alternative_point_bypass'][bypass(diagonal['local_alternative_tolerance'],diagonal['observed_component_alternative_tolerance'])]+=1
            result['range_component_bounds'][diagonal['component_range_certificate']['status']]+=1
            result['range_local_bounds'][diagonal['local_range_certificate']['status']]+=1
            for s in diagonal['preparation_scenarios']:
                key=s['protein_column']+'|'+s['growth_column']
                counts=scenarios.setdefault(key,{'component_bounds':collections.Counter(),'bypass':collections.Counter(),'alternative_bypass':collections.Counter()})
                counts['component_bounds'][s['component']['bounds']['status']]+=1
                counts['bypass'][bypass(s['local'],s['component'])]+=1
                counts['alternative_bypass'][bypass(s['local_alternative_tolerance'],s['component_alternative_tolerance'])]+=1
    return {'counts':result,'per_scenario':scenarios}

def main():
    for name,(path,pin) in PINS.items():check(sha(path)==pin,'Held input drift: '+name)
    sandbox=ROOT/'synthetic_run';sandbox.mkdir(exist_ok=False)
    tmp=sandbox/'tmp';tmp.mkdir()
    for name in ('diagnostics.py','test_diagnostics.py'):shutil.copyfile(PINS[name][0],sandbox/name)
    environment=dict(os.environ,TMP=str(tmp),TEMP=str(tmp),PYTHONDONTWRITEBYTECODE='1')
    completed=subprocess.run([sys.executable,'-I','-B',str(sandbox/'test_diagnostics.py')],env=environment,capture_output=True,text=True)
    (ROOT/'author_tests_stdout.txt').write_text(completed.stdout,encoding='utf-8')
    (ROOT/'author_tests_stderr.txt').write_text(completed.stderr,encoding='utf-8')
    check(completed.returncode==0,'Held author tests failed')
    sys.dont_write_bytecode=True
    sys.path.insert(0,str(sandbox))
    tests=module('held_synthetic_fixtures',sandbox/'test_diagnostics.py');d=tests.d
    split,e=tests.fixtures();modified=copy.deepcopy(e)
    for row in modified['prediction_evaluation']['records']:
        if row['status']!='scored':continue
        base=[F(-1,2),F(5,2),F(7)][row['square_index']] if row['group']=='1' else F(2)
        for j,method in enumerate(('additive','adjacent_mean','three_corner_mean')):
            errors=[float(base*2**j),float(base*2**(j+1))]
            row['methods'][method].update(predicted_target=[y+x for y,x in zip(row['observed_target'],errors)],
                signed_error=errors,absolute_error=list(map(abs,errors)),squared_error=[x*x for x in errors])
    result=d.analyze(modified,split,'development')
    rows=[dict(r,withheld_ID=split['units'][r['square_index']]['withheld_ID']) for r in modified['prediction_evaluation']['records']]
    comparisons=0
    summaries=[('full',list(map(str,range(1,17))),result['full_arm'])]
    summaries += [(r['omitted_group'],[str(g) for g in range(1,17) if str(g)!=r['omitted_group']],r) for r in result['leave_one_group_out']]
    for label,groups,payload in summaries:
        for target,key in ((False,'case_balanced_prediction'),(True,'target_balanced_prediction')):
            for method in ('additive','adjacent_mean','three_corner_mean'):
                for k,trait in enumerate(('protein','growth')):
                    ref=reference_loss(rows,groups,method,k,target,False)
                    mse=reference_loss(rows,groups,method,k,target,True)
                    actual=payload[key]['summaries'][method]['continuous'][trait]
                    check(near(actual['MAE'],ref),label+' MAE mismatch')
                    check(near(actual['RMSE']**2,mse),label+' MSE mismatch')
                    comparisons+=2
    primary=result['full_arm']['case_balanced_prediction']['summaries']['additive']['continuous']['protein']
    target=result['full_arm']['target_balanced_prediction']['summaries']['additive']['continuous']['protein']
    check(near(primary['MAE'],F(94,45)) and near(target['MAE'],F(43,20)),'Independent hand example mismatch')
    check(near(primary['RMSE']**2,F(149,30)) and near(target['RMSE']**2,F(219,40)),'Hand squared example mismatch')
    check(near(result['full_arm']['descriptive']['equal_contributing_group_mean_fractions']['joint_only'],F(1,18)),'Equal group fraction mismatch')
    # All missing outcomes keep 16 empty group rows and 16 null deletion rows.
    empty=copy.deepcopy(e)
    for row in empty['coverage']:row['eligible']=False
    for row in empty['descriptions']:row.update(eligible=False,diagonals=[])
    for row in empty['prediction_evaluation']['records']:row.clear()
    empty['prediction_evaluation']['records']=[{'square_index':u['square_index'],'seed':u['seed'],'group':u['seed'],'status':'missing_known_replicates'} for u in split['units'] if u['arm']=='development']
    empty_result=d.analyze(empty,split,'development')
    leaf='case_balanced_prediction/summaries/additive/continuous/protein/MAE'
    check(empty_result['all_leave_one_out_minmax_ranges'][leaf]=={'minimum':None,'maximum':None,'defined_rows':0,'undefined_rows':16,'all_rows':16},'Null extrema not preserved')
    check(empty_result['full_arm']['descriptive']['contributing_groups']==0,'Empty groups count mismatch')
    missing=copy.deepcopy(e);missing['prediction_evaluation']['records'].pop()
    try:d.analyze(missing,split,'development')
    except ValueError:pass
    else:raise RuntimeError('Missing prediction record silently omitted')
    # The genuine provenance-hold schema must not call the numerical analyzer.
    sp=sandbox/'toy_split.json';sp.write_text(json.dumps(split),encoding='utf-8')
    hold={'arm':'development','status':'held_for_author_mean_provenance','interpretation_permitted':False,
        'source_hashes':{'runner':d.RUNNER_SHA256,'split':sha(sp)},'holds':['SYNTHETIC_MEAN_HOLD']}
    ep=sandbox/'toy_hold.json';ep.write_text(json.dumps(hold),encoding='utf-8')
    with mock.patch.object(d,'analyze',side_effect=RuntimeError('Analyzer accessed on provenance hold')):
        check(d.run(ep,sha(ep),sp,'development',sandbox/'held')=='provenance_hold_no_interpretation','Provenance hold failed')
    counts=summarize_components(e['descriptions'])
    check(counts['counts']['point_component_bounds']=={'inconclusive':12,'certified_jointly_connected':12,'certified_joint_only_obstruction':10},'Reference component counts mismatch')
    check(counts['counts']['point_bypass']=={'observed_component_bypass_unresolved':2,'observed_component_optimistic_disconnection':1,'certified_observed_component_bypass':2,'not_a_local_joint_only_obstruction':29},'Reference bypass counts mismatch')
    check(set(d.analyze(e,split,'development')['observed_component_diagnostics'])=={'eligible_diagonals','flip_counts','records','scope'},'Candidate unexpectedly contains an aggregate correction')
    write(ROOT/'expected_component_count_fixture.json',counts)
    write(ROOT/'independent_results.json',{'verdict':'BLOCKED_INCOMPLETE_AMENDMENT_COUNTS',
        'author_tests_passed':14,'independent_rational_numeric_comparisons':comparisons,
        'hand_values':{'case_MAE':'94/45','target_MAE':'43/20','case_MSE':'149/30','target_MSE':'219/40','equal_group_joint_fraction':'1/18'},
        'all16_null_extrema_preserved':True,'missing_prediction_record_refused':True,'provenance_hold_did_not_access_analyzer':True,
        'component_aggregate_counts_absent':True,'no_real_outcomes_read':True,'shared_project_context_not_blinded':True,
        'inputs':{name:{'sha256':pin,'unchanged_after_review':sha(path)==pin} for name,(path,pin) in PINS.items()}})
    for name,(path,pin) in PINS.items():check(sha(path)==pin,'Input changed during review: '+name)
    print(json.dumps({'verdict':'BLOCKED_INCOMPLETE_AMENDMENT_COUNTS','numeric_comparisons':comparisons,'real_outcomes_read':0}))

if __name__=='__main__':main()
