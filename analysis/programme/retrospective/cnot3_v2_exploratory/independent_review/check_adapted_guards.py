"""Independently specified guard-only fixtures. No observed file input/model fit."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
IMPL=HERE.parent/'implementation'
PINS={'pipeline.py':'5674bb2c9f7f07b47ac4146dbf9717be2701c98849b31fc39fc7762ab4fd07f6',
      'predict.py':'be75c90e4076f29952553fa2a22529c2640db0691cc8237e5078c1d1e2cc25d7'}
CHECKS=['count_5480','distinct_stable_ids','positive_finite_paired_half_lives','positive_finite_SE_IP','STAT_IP_arithmetic',
        'no_negative_SE_among_IP_candidates','no_top1000_boundary_tie','median_log2_change_top1000_greater_rest',
        'ten_significant_digit_ratio_compatibility','all_49_original_ratio_discrepancy_genes_retained','independent_row_join_verified']

def load(name):
    path=IMPL/name
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PINS[name]:raise ValueError('Candidate changed')
    spec=importlib.util.spec_from_file_location('review_'+path.stem,path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def main():
    p,model=load('pipeline.py'),load('predict.py')
    ids=[f'ENSG{i+1:011d}' for i in range(5480)]
    before={'rows':'a'*64,'original_gate':'b'*64,'source_audit':'c'*64,'protocol':'d'*64}
    old=dict(schema='CNOT3_principal_reproduction_gate_v1',status='held_for_source_diagnosis',cohort_rows_sha256=before['rows'],
             cohort_gene_count=5480,cohort_GeneIDs=ids,principal_discrepancies=['stored_HL_ratio_inconsistent:'+g for g in ids[:49]]+['published_cohort_size_mismatch'])
    gate=dict(schema='CNOT3_operational_source_gate_v2',status='operational_cohort_verified',original_frozen_gate_status='held_for_source_diagnosis',
              historical_published_cohort_reproduced=False,reproduction_gate_passed=False,operational_source_gate_passed=True,
              unresolved_source_discrepancies=['published_cohort_size_mismatch'],cohort_gene_count=5480,cohort_rows_sha256=before['rows'],
              prediction_protocol_sha256=before['protocol'],original_frozen_gate_sha256=before['original_gate'],independent_source_audit_sha256=before['source_audit'],
              operational_checks={k:True for k in CHECKS},cohort_GeneIDs=ids,original_ratio_discrepancy_GeneIDs=ids[:49])
    p.verify_operational_gate(gate,old,before)
    tests=[{'case':'exact_operational_guard','passed':True}]
    variants=[]
    for key,value in [('historical_published_cohort_reproduced',True),('reproduction_gate_passed',True),('original_frozen_gate_status','passed'),
                      ('cohort_gene_count',5481),('cohort_rows_sha256','f'*64),('independent_source_audit_sha256','e'*64),
                      ('unresolved_source_discrepancies',[]),('original_ratio_discrepancy_GeneIDs',ids[:48])]:
        g=copy.deepcopy(gate);g[key]=value;variants.append((key,g,old))
    for key in CHECKS:
        g=copy.deepcopy(gate);g['operational_checks'][key]=False;variants.append(('failed_check_'+key,g,old))
    g=copy.deepcopy(gate);g['operational_checks']['unknown']=True;variants.append(('unknown_check',g,old))
    o=copy.deepcopy(old);o['status']='passed';variants.append(('original_reclassified_pass',gate,o))
    o=copy.deepcopy(old);o['principal_discrepancies'].append('unexpected_failure');variants.append(('unexpected_old_failure',gate,o))
    o=copy.deepcopy(old);o['cohort_rows_sha256']='e'*64;variants.append(('changed_original_rows',gate,o))
    for name,g,o in variants:
        try:p.verify_operational_gate(g,o,before)
        except (p.c.IntegrityError,ValueError):tests.append({'case':name,'passed':True})
        else:raise AssertionError('Gate accepted '+name)
    # Count/partition guard only; no feature construction, outcome values or fit.
    groups=[f'SYN_GROUP_{i:05d}' for i in range(5480)]
    records=[dict(GeneID='SYN_GENE_'+str(i),group_ID=g) for i,g in enumerate(groups)]
    train=[r for r in records if model.group_fold(r['group_ID'],'outer')!=0]
    holdout=[r for r in records if model.group_fold(r['group_ID'],'outer')==0]
    receipt=dict(status='root_authorized_operational_prediction',analysis_status='post_exposure_exploratory',
                 historical_published_cohort_reproduced=False,original_frozen_gate_status='held_for_source_diagnosis',
                 operational_source_gate_passed=True,reference_eligibility_passed=True,grouping_complete=True,
                 operational_source_gate_sha256='1'*64,original_frozen_gate_sha256='2'*64,independent_source_audit_sha256='3'*64,
                 operational_slam_overlap_count=5480,eligible_count=5480)
    model.check_design(train,holdout,receipt);tests.append({'case':'operational_prediction_design_only','passed':True})
    for key,value in [('status','root_authorized_prediction'),('analysis_status','confirmatory'),('historical_published_cohort_reproduced',True),
                      ('original_frozen_gate_status','passed'),('reproduction_gate_passed',True),('operational_source_gate_passed',False),
                      ('operational_source_gate_sha256','not-a-hash')]:
        r=copy.deepcopy(receipt);r[key]=value
        try:model.check_design(train,holdout,r)
        except ValueError:tests.append({'case':'prediction_design_refuses_'+key,'passed':True})
        else:raise AssertionError('Prediction guard accepted '+key)
    out=HERE/'guard_reference.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps({'source_pins':PINS,'tests_passed':len(tests),'tests':tests,'actual_modeling':False,'source_outcome_reads':False},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'guard_checks_passed':len(tests),'actual_fits':0}))

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Optimized Python refused')
    main()
