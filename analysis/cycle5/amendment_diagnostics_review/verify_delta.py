"""Narrow synthetic-only delta check. Preserve the original blocked review."""
import ast
import copy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
AUTHOR=HERE.parent/'amendment_diagnostics'
CODE='5e8c1a5638d0d2746f89f12cf26d209f11930f490df99c16a79eda3cf2628850'
TEST='cefc66dc53a45035cc4f3db1325166488055f39b922335c8b1ec485af11548bb'
OLD='9e0307ad72b6dd0be771dee9d886bb87e63d09ed85a7ceea5409b2e3f90f897d'

def check(value,message):
    if not value:raise RuntimeError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module);return module
def write(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')

def main():
    initial={p.name:sha(p) for p in (HERE/'review.md',HERE/'receipt.json',HERE/'independent_results.json',HERE/'expected_component_count_fixture.json')}
    check(sha(AUTHOR/'diagnostics.py')==CODE and sha(AUTHOR/'test_diagnostics.py')==TEST,'Corrected hash mismatch')
    check(sha(AUTHOR/'history_v1/diagnostics.py')==OLD,'Original source history lost')
    check(sha(AUTHOR/'history_v1/independent_blocked_review.md')==sha(HERE/'review.md'),'Original BLOCKED review history changed')
    scratch=HERE/'delta_scratch';scratch.mkdir(exist_ok=False);(scratch/'history_v1').mkdir()
    for name in ('diagnostics.py','test_diagnostics.py'):shutil.copyfile(AUTHOR/name,scratch/name)
    shutil.copyfile(AUTHOR/'history_v1/diagnostics.py',scratch/'history_v1/diagnostics.py')
    shutil.copyfile(HERE/'expected_component_count_fixture.json',scratch/'reviewer_expected_component_count_fixture.json')
    completed=subprocess.run([sys.executable,'-I','-B',str(scratch/'test_diagnostics.py')],capture_output=True,text=True,
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    (HERE/'delta_tests_stdout.txt').write_text(completed.stdout,encoding='utf-8')
    (HERE/'delta_tests_stderr.txt').write_text(completed.stderr,encoding='utf-8')
    check(completed.returncode==0 and 'Ran 18 tests' in completed.stderr,'18 corrected tests did not pass')
    sys.dont_write_bytecode=True
    d=load_module('diagnostics',scratch/'diagnostics.py');sys.modules['diagnostics']=d
    tests=load_module('delta_toys',scratch/'test_diagnostics.py')
    old=load_module('delta_old',scratch/'history_v1/diagnostics.py')
    ref=load_module('retained_independent_reference',HERE/'verify_independently.py')
    before_ast={n.name:n for n in ast.parse((scratch/'history_v1/diagnostics.py').read_text()).body if isinstance(n,ast.FunctionDef)}
    after_ast={n.name:n for n in ast.parse((scratch/'diagnostics.py').read_text()).body if isinstance(n,ast.FunctionDef)}
    changed=[name for name,node in before_ast.items() if ast.dump(node)!=ast.dump(after_ast[name])]
    check(changed==['diagonal_diagnostics'],'Existing numerical or guard functions changed')
    check(ast.dump(ast.Module(body=before_ast['diagonal_diagnostics'].body[:-1],type_ignores=[]))==
          ast.dump(ast.Module(body=after_ast['diagonal_diagnostics'].body[:-1],type_ignores=[])),'Old diagonal record or flip logic changed')
    split,e=tests.fixtures();original=old.analyze(e,split,'development');result=d.analyze(e,split,'development')
    for key,value in original.items():
        if key=='observed_component_diagnostics':
            for sub,subvalue in value.items():check(result[key][sub]==subvalue,'Existing diagonal output changed: '+sub)
        else:check(result[key]==value,'Existing numeric output changed: '+key)
    expected=json.loads((HERE/'expected_component_count_fixture.json').read_text())
    table=result['observed_component_diagnostics']['count_tables'];checks=0
    def agree(actual,wanted,n=34):
        nonlocal checks
        check(actual['denominator']==n and sum(actual['counts'].values())==n,'Count denominator mismatch')
        check({k:v for k,v in actual['counts'].items() if v}==wanted,'Independent expected category counts differ')
        checks+=1
    for tol,bypass_key in [('zero_tolerance','point_bypass'),('alternative_tolerance','alternative_point_bypass')]:
        point=table['primary_point_field'][tol]
        agree(point['component_bound_statuses'],expected['counts']['point_component_bounds'])
        agree(point['component_point_classifications'],expected['counts']['point_component_classification'])
        agree(point['bypass_outcomes'],expected['counts'][bypass_key])
        denominator=5 if tol=='zero_tolerance' else 4
        check(point['bypass_outcomes']['conditional_local_joint_only_denominator']==denominator,'Conditional denominator changed')
        for key,v in point['bypass_outcomes']['conditional_fractions'].items():check(v==expected['counts'][bypass_key][key]/denominator,'Conditional fraction mismatch')
        rt=table['empirical_range_field'][tol]
        agree(rt['local_certificate_statuses'],expected['counts']['range_local_bounds'])
        agree(rt['component_certificate_statuses'],expected['counts']['range_component_bounds'])
    check(set(table['named_preparation_scenarios'])==set(expected['per_scenario']) and table['pooled_across_scenarios'] is False,'Named scenario scope changed')
    for name,ex in expected['per_scenario'].items():
        for tol,bypass_key in [('zero_tolerance','bypass'),('alternative_tolerance','alternative_bypass')]:
            t=table['named_preparation_scenarios'][name][tol]
            agree(t['component_bound_statuses'],ex['component_bounds'])
            agree(t['component_point_classifications'],expected['counts']['point_component_classification'])
            agree(t['bypass_outcomes'],ex[bypass_key])
    # Independent unavailable/no-local boundary cases are not reclassified as zero effects.
    boundary=d.component_count_table([(tests.classified('joint_only_obstruction'),None),(None,None),
        (tests.classified('jointly_connected'),tests.field(connected=True))])
    check(boundary['eligible_diagonals']==3 and boundary['unavailable_component_points']==2,'Unavailable point accounting')
    check(boundary['component_bound_statuses']['counts']['certificate_unavailable']==2,'Unavailable bound accounting')
    check(boundary['bypass_outcomes']['counts']['observed_component_bypass_unresolved']==1 and
          boundary['bypass_outcomes']['counts']['local_classification_unavailable']==1 and
          boundary['bypass_outcomes']['counts']['not_a_local_joint_only_obstruction']==1,'Unavailable bypass categories')
    no_local=d.component_count_table([(tests.classified('jointly_connected'),tests.field(connected=True))])
    check(all(v is None for v in no_local['bypass_outcomes']['conditional_fractions'].values()),'No-local denominator not null')
    empty=d.aggregate_component_counts([])
    check(len(empty['named_preparation_scenarios'])==12,'Empty named scenarios dropped')
    for context in list(empty['primary_point_field'].values())+[s[t] for s in empty['named_preparation_scenarios'].values() for t in ('zero_tolerance','alternative_tolerance')]:
        check(context['eligible_diagonals']==0,'Empty context denominator')
        for key in ('local_classifications','component_point_classifications','component_bound_statuses','bypass_outcomes'):
            check(context[key]['denominator']==0 and sum(context[key]['counts'].values())==0,'Empty categories not explicit zero')
    # Preserve original408 independent exact-rational comparisons on corrected source.
    for row in e['prediction_evaluation']['records']:
        if row['status']!='scored':continue
        base=[F(-1,2),F(5,2),F(7)][row['square_index']] if row['group']=='1' else F(2)
        for j,method in enumerate(d.METHODS):
            errors=[float(base*2**j),float(base*2**(j+1))]
            row['methods'][method].update(predicted_target=[y+x for y,x in zip(row['observed_target'],errors)],
                signed_error=errors,absolute_error=list(map(abs,errors)),squared_error=[x*x for x in errors])
    result=d.analyze(e,split,'development')
    rows=[dict(r,withheld_ID=split['units'][r['square_index']]['withheld_ID']) for r in e['prediction_evaluation']['records']]
    payloads=[(list(map(str,range(1,17))),result['full_arm'])]+[
        ([str(g) for g in range(1,17) if str(g)!=r['omitted_group']],r) for r in result['leave_one_group_out']]
    numeric=0
    for groups,payload in payloads:
        for target,key in ((False,'case_balanced_prediction'),(True,'target_balanced_prediction')):
            for method in d.METHODS:
                for k,trait in enumerate(d.TRAITS):
                    actual=payload[key]['summaries'][method]['continuous'][trait]
                    check(ref.near(actual['MAE'],ref.reference_loss(rows,groups,method,k,target,False)),'Rational MAE changed')
                    check(ref.near(actual['RMSE']**2,ref.reference_loss(rows,groups,method,k,target,True)),'Rational MSE changed')
                    numeric+=2
    check(numeric==408,'Reference comparison coverage changed')
    for name,h in initial.items():check(sha(HERE/name)==h,'Original review artifact changed: '+name)
    check(sha(AUTHOR/'diagnostics.py')==CODE and sha(AUTHOR/'test_diagnostics.py')==TEST,'Corrected source drift during review')
    write(HERE/'delta_results.json',{'verdict':'PASS_EXACT_CORRECTED_REVISION','code_sha256':CODE,'tests_sha256':TEST,
        'author_tests':18,'independent_rational_comparisons':numeric,'independent_category_table_comparisons':checks,
        'old_function_AST_changes':changed,'old_numeric_record_outputs_identical':True,
        'unavailable_empty_no_local_boundaries_pass':True,'initial_blocked_history_preserved':True,
        'no_new_scientific_criteria':True,'real_outcomes_read':0,'shared_project_context_not_blinded':True})
    print(json.dumps({'verdict':'PASS_EXACT_CORRECTED_REVISION','numeric_checks':numeric,'category_table_checks':checks}))

if __name__=='__main__':main()
