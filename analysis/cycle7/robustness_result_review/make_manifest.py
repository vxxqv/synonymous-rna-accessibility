"""Summarize already checked evidence and produce its computation manifest."""
from pathlib import Path
import csv
import hashlib
import json

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    audit=json.loads((BASE/'observed_v1/receipt.json').read_bytes())
    target=BASE/'preparation_radii.csv'
    if target.exists():raise RuntimeError('Refusing existing table')
    with target.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['kind','omitted_zero_based','retained_scenarios','epsilon_exact','gamma_squared_exact','gamma_lower_exact','gamma_upper_exact'])
        for row in audit['preparation_deletions']:
            a=row['radii']['additive'];m=row['radii']['multiplicative_squared'];g=m['gamma_enclosure']
            writer.writerow([row['kind'],row['omitted_zero_based'],row['retained_scenario_count'],a['parameter_exact'],m['parameter_exact'],g['decimal_lower_exact'],g['decimal_upper_exact']])
    artifacts=[]
    for name,digest in audit['input_hashes_before_and_after'].items():
        path=Path(name)
        if sha(path)!=digest:raise RuntimeError('Input drift: '+name)
        artifacts.append({'path':path.as_posix(),'sha256':digest,'sha256_after':digest})
    result_path=BASE/'observed_v1/receipt.json'
    manifest={'schema_version':2,'claim_id':'exact_observed_anchored_robustness_certificate_review',
       'repository':{'commit':'7eae1381c781c93e324f01d4e0d07692bceabea5','dirty':True},
       'command':['work/env/Scripts/python.exe','work/research_cycle7/robustness_result_review/verify_result.py',
                  '--repo','outputs/synonymous-rna-accessibility','--result','work/research_cycle7/results/robustness.json',
                  '--out','work/research_cycle7/robustness_result_review/observed_v1'],
       'environment':{'software':[{'name':'CPython','version':'3.12.14'}],'hardware':audit['platform']+'; single-process exact arithmetic'},
       'mathematics':{'assertion_tested':'All26 claimed radii have exact sharp lower and feasible upper certificates; all72 row and12 preparation deletions match the frozen domain.',
                      'coefficient_domain':'Rational scores and coefficients, integer feature design; rational enclosures of square roots',
                      'conventions':'Nonnegative closed intervals; native reference46 anchored at1; common nondecreasing response; no reoptimization.',
                      'inputs':[x['path'] for x in artifacts],
                      'bounds':{'radii':26,'row_deletions':72,'preparation_deletions':12,'search_event_certificates':220},
                      'non_claims':['No independent biological validation','No calibrated error budget','No causal or barrier attribution']},
       'randomness':{'used':False,'generator':'','seed':None},
       'run':{'started_at':audit['started_at'],'runtime_seconds':audit['runtime_seconds'],'exit_status':0,'status':'completed',
              'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':1,'max_threads':1}},
       'outputs':[{'path':result_path.relative_to(ROOT).as_posix(),'sha256':sha(result_path),'kind':'result'}],
       'checks':['Exact source and output hashes before and after','Independent interval and feature reconstruction',
                 'All primal/dual rational certificates','Every positive-weight lower comparison lifetime',
                 'All declared deletion identities and nesting checks','Reference omission separated from measured observations'],
       'result':'Implementation-independent exact verification of all declared saved finite-domain certificates completed successfully.',
       'residual_risks':['Unmeasured recovery bias is not calibrated','Scenario dependence not modeled by Cartesian interval widening',
                         'Review is separately authored with inherited context'],
       'input_artifacts':artifacts,'declared_results':[result_path.relative_to(ROOT).as_posix()]}
    dest=BASE/'computation_manifest.json'
    if dest.exists():raise RuntimeError('Refusing existing manifest')
    dest.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'manifest_sha256':sha(dest),'table_sha256':sha(target),'review_sha256':sha(BASE/'review.md')}))

if __name__=='__main__':main()
