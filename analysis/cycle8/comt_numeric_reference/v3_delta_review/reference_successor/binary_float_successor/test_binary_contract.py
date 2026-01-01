"""Invented non-dyadic inputs only. No observed fitting or source outcome access."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import ast
import copy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import random
import sys
sys.dont_write_bytecode = True
import reference_binary_float as ref


def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def function_sources(path):
    s=Path(path).read_text(encoding='utf8')
    return {n.name:ast.get_source_segment(s,n) for n in ast.parse(s).body if isinstance(n,ast.FunctionDef)}


def main(pilot, predecessor, invented_stage, out):
    pilot, predecessor, invented_stage, out=map(Path,(pilot,predecessor,invented_stage,out))
    if out.exists():raise FileExistsError(out)
    assert ref.sha(pilot/'proxy_models.py')=='8811af9750e813c798122c7e936b6279a0f3b0e53fea1924cb2d633162a34093'
    assert ref.sha(pilot/'implementation_v3/prediction.py')=='92d0d46fda3b7b3c49c5a4d8780c8b57bcf9c15d9c8095aeca27555f38e2dede'
    primary=load(pilot/'implementation_v3/prediction.py','frozen_primary_adapter_synthetic_only')
    kernel=load(pilot/'proxy_models.py','frozen_primary_kernel_synthetic_only')
    old=load(predecessor,'preserved_source_first_reference')
    before,after=function_sources(predecessor),function_sources(ref.__file__)
    changed=sorted(n for n in before if before[n]!=after[n])
    assert changed==['compare_downstream','rational','run','scientific_comparison','verify_provenance']
    assert ref.rational(0.1)==F(0.1) and ref.rational(0.1)!=F('0.1')
    assert ref.rational('0.1')==F(1,10) and ref.rational(F(1,10))==F(1,10)
    assert ref.rational(10)==10 and ref.rational('1/3')==F(1,3)
    for invalid in (True,False,None,float('nan'),float('inf'),float('-inf')):
        try:ref.rational(invalid)
        except ValueError:pass
        else:raise AssertionError('Invalid finite-number input accepted')
    rng=random.Random(321903)
    cases=[]
    for case in range(46):
        train=[]
        for i in range(24):
            x=rng.choice([-.7,-.3,-.1,.1,.2,.6]);z=rng.choice([-.9,-.2,.1,.3,.8])
            y=.1+.3*x-.7*z+rng.choice([-.1,0.,.2])
            if case==0:z=2*x
            if case==1:x=.1
            if case==2:z=.2
            if case==3:x,z=.1,.2
            if case==4:y=-3*x+.1*z
            if case==5:x,z=(i%4)/8,(i//4%3)/8;y=-3*x+.1*z
            train.append({'hgvs_nt':f'INVENTED_T{case}_{i}','site':str(i//2),'ROI':'ROI_1' if i%2 else 'ROI_2','R':x,'L':z,'T':y})
        assess=[{'hgvs_nt':f'INVENTED_A{j}','site':'new'+str(j),'ROI':'ROI_1' if j%2 else 'ROI_2','R':x,'L':z}
                for j,(x,z) in enumerate([(-.95,-2.),(-.1,.15),(.1,.7),(1.2,2.5)])]
        independent, fit=ref.fit_and_predict(train,assess)
        production=primary.fit_predict(train,assess,'L',kernel)
        for p in production['predictions']:
            assert independent[p['hgvs_nt']]==p['predictions']
        for model in ('RNA_linear','RNA_load_linear'):
            actual=production['fit'][model]
            assert fit[model]==(actual['intercept'],actual['rna_slope'],actual['load_slope'],actual['centered_rank'])
        cases.append({'case':case,'predictions_checked':len(assess)*5,'coefficient_rank_sets_checked':2})
    # Explicit unrelated-load extrapolation: the RNA-only model does not use that predictor.
    train=[{'hgvs_nt':str(i),'site':str(i),'ROI':'ROI_1' if i%2 else 'ROI_2','R':i*.1,'L':i*.2,'T':i*.3+.1} for i in range(8)]
    assess=[{'hgvs_nt':'INVENTED_OUTSIDE_LOAD','site':'new','ROI':'ROI_1','R':.3,'L':100.1}]
    p=primary.fit_predict(train,assess,'L',kernel)['predictions'][0]
    assert p['range_flags']['RNA_linear']=={'rna_outside':False,'load_outside':False}
    assert p['range_flags']['RNA_load_linear']=={'rna_outside':False,'load_outside':True}
    # Complete invented 125-ID/53-site CV exercise with non-dyadic values and an outlier load.
    data=json.loads((invented_stage/'normalization/contrasts.json').read_text())
    for mapping,rows in data.items():
        for i,row in enumerate(rows):
            values={'R':.1+.3*row['contrasts']['R']['center'],
                    'T':.2+.7*row['contrasts']['T']['center']+.1*(i%3)+(0 if mapping=='code' else .2*(i%2)),
                    'L':.1+.2*row['contrasts']['L']['center'],
                    'L_mid':.3+.7*row['contrasts']['L_mid']['center']}
            if i==0:values['L']=100.1
            for k,v in values.items():
                lo,hi=v-.07,v+.13
                row['contrasts'][k]={'center':v,'lower':lo,'upper':hi,'complete':True,
                    'resolved_sign':1 if lo>0 else -1 if hi<0 else 0}
    own=ref.downstream(data)
    pred=primary.cross_validate(data,kernel);evaluation=primary.evaluate(pred,data,kernel)
    all_checks=ref.compare_downstream(own,pred,evaluation)
    # Scientific exact comparisons are preserved even for invisible float-display changes.
    altered=copy.deepcopy(pred)
    record=altered['code']['heavy']['predictions'][0]['predictions']['pooled_median']
    record['exact_stored_input_arithmetic']=str(F(record['exact_stored_input_arithmetic'])+F(1,10**40))
    try:ref.compare_downstream(own,altered,evaluation)
    except AssertionError as e:assert 'Exact prediction mismatch' in str(e)
    else:raise AssertionError('Tiny exact discrepancy was silently tolerated')
    # Incorrect RNA-only and joint range flags must remain rejected.
    for model,flag in [('RNA_linear','rna_outside'),('RNA_linear','load_outside'),('RNA_load_linear','load_outside')]:
        altered=copy.deepcopy(pred)
        flags=altered['code']['heavy']['predictions'][0]['range_flags'][model]
        flags[flag]=not flags[flag]
        try:ref.compare_downstream(own,altered,evaluation)
        except AssertionError as e:assert 'flags mismatch' in str(e)
        else:raise AssertionError('Invalid range flag accepted')
    # Entirely constant non-dyadic response gives zero baselines and cannot pass the gain rule.
    constant=copy.deepcopy(data)
    for rows in constant.values():
        for row in rows:row['contrasts']['T']={'center':.1,'lower':.1,'upper':.1,'complete':True,'resolved_sign':1}
    own_constant=ref.downstream(constant)
    pp=primary.cross_validate(constant,kernel);ee=primary.evaluate(pp,constant,kernel)
    constant_checks=ref.compare_downstream(own_constant,pp,ee)
    assert own_constant['primary_operational_expansion_pass'] is False
    for mapping in ('code','geo'):
        for mode in ('heavy','mid'):
            assert not own_constant['cases'][mapping][mode]['criterion']['relative_comparison_defined']
    # Missingness remains explicit and fixed sites/IDs are not relabeled.
    missing=copy.deepcopy(data)
    for k in ('R','T'):
        missing['code'][1]['contrasts'][k]={'center':None,'lower':None,'upper':None,'complete':False,'resolved_sign':None}
    missing['geo'][2]['contrasts']['L']={'center':None,'lower':None,'upper':None,'complete':False,'resolved_sign':None}
    own_missing=ref.downstream(missing)
    pp=primary.cross_validate(missing,kernel);ee=primary.evaluate(pp,missing,kernel)
    missing_checks=ref.compare_downstream(own_missing,pp,ee)
    assert before['normalize']==after['normalize'] and before['contrast_summaries']==after['contrast_summaries']
    assert before['criterion']==after['criterion'] and before['constrained_linear']==after['constrained_linear']
    out.mkdir(parents=True)
    receipt={'status':'PASS_INVENTED_BINARY_FLOAT_AND_DIAGNOSTIC_CONTRACT_CORRECTIONS',
        'changed_common_functions':changed,'unchanged_common_functions':sorted(k for k in before if before[k]==after[k]),
        'direct_non_dyadic_fit_cases':len(cases),'direct_model_predictions_checked':sum(c['predictions_checked'] for c in cases),
        'direct_coefficient_rank_sets_checked':sum(c['coefficient_rank_sets_checked'] for c in cases),
        'non_dyadic_full_CV_checks':all_checks,'constant_zero_baseline_CV_checks':constant_checks,
        'explicit_missingness_CV_checks':missing_checks,'RNA_only_load_flag_false_joint_load_flag_true':True,
        'incorrect_range_flag_falsifiers_rejected':3,'arbitrarily_tiny_exact_mismatch_rejected':True,
        'raw_decimal_strings_and_existing_fractions_preserved':True,'actual_observed_analysis_or_replay_executed':False,
        'reference_sha256':ref.sha(ref.__file__),'test_sha256':ref.sha(__file__),
        'production_kernel_and_adapter_source_hashes':{'kernel':ref.sha(pilot/'proxy_models.py'),'adapter':ref.sha(pilot/'implementation_v3/prediction.py')}}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',required=True);p.add_argument('--predecessor',required=True);p.add_argument('--invented-stage',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.pilot,a.predecessor,a.invented_stage,a.out)
