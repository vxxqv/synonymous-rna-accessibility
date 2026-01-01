"""Record the completed synthetic preflight; no observed source is opened."""
from pathlib import Path
import hashlib
import json
import platform

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def relative(path):return path.relative_to(ROOT).as_posix()

def main():
    result_path=BASE/'synthetic_v1.json'
    result=json.loads(result_path.read_bytes())
    inputs=[]
    for name,digest in result['source_hashes'].items():
        path=Path(name)
        if sha(path)!=digest:raise RuntimeError('Source drift')
        inputs.append({'path':relative(path),'sha256':digest,'sha256_after':sha(path)})
    manifest={
      'schema_version':2,'claim_id':'finite_scenario_accessibility_quantifiers_v1',
      'repository':{'commit':'9ed05ad8ab75c596e87a5ce7b63eb1bb38e71ba5','dirty':True},
      'command':['work/env/Scripts/python.exe','work/research_cycle7/scenario_geometry/verify_synthetic.py',
                 '--out','work/research_cycle7/scenario_geometry/synthetic_v1.json'],
      'environment':{'software':[{'name':'CPython','version':'3.12.14'}],
                     'hardware':platform.machine()+'; '+platform.platform()+'; single-process execution'},
      'mathematics':{'assertion_tested':result['assertion_tested'],
                     'coefficient_domain':'Exact rational scores with symbolic infinities for disconnected/unknown bounds',
                     'conventions':'Paths include endpoints; acceptance >= threshold; nonnegative common absolute thresholds; all distinct unordered pairs.',
                     'inputs':[item['path'] for item in inputs],
                     'bounds':{'complete_vertices':3,'complete_graphs':8,'complete_scenarios':2,'score_alphabet':[0,1,2],
                               'unknown_vertices':3,'unknown_scenarios':2,'unknown_alphabet':[0,1,None],
                               'binary_completions':[0,1],'full_shape_constant_fixture':result['full_shape']},
                     'non_claims':['No new observed summaries','No independently authored review','No novelty or causal claim','No biological uncertainty or independence estimate']},
      'randomness':{'used':False,'generator':'','seed':None},
      'run':{'started_at':result['started_at'],'runtime_seconds':result['runtime_seconds'],
             'exit_status':0,'status':'completed',
             'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':1,'max_threads':1}},
      'outputs':[{'path':relative(result_path),'sha256':sha(result_path),'kind':'result'}],
      'checks':['Independent simple-path exact capacity oracle','Direct path/scenario acceptance oracle',
                'Unknown-completion containment and certificate checks','All threshold ties and empty intervals retained',
                'All72/27/2556 coverage fixture','Source hashes unchanged during run'],
      'result':'Implementation and finite assertion verified in the explicitly stated synthetic range.',
      'residual_risks':['Same author for implementation and oracle','Unknown-value finite tests use only binary completions; general guarantee relies on proof',
                        'Observed runner and independent review remain pending'],
      'input_artifacts':inputs,'declared_results':[relative(result_path)]}
    candidate={
      'status':'candidate_awaiting_independent_review_and_root_freeze',
      'observed_execution_authorized':False,'new_observed_summaries_computed':False,
      'files':{p.name:sha(p) for p in [BASE/'protocol.md',BASE/'derivation.md',BASE/'scenario_geometry.py',
                                      BASE/'verify_synthetic.py',result_path,BASE/'preflight_report.md']},
      'expected_future_source':{'compressed_sha256':'e3c3d061052643046de58dc789866899b0fb14ade83f06170e2cf32e9e54e145',
                                'decompressed_sha256':'3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466',
                                'source_not_opened_during_preflight':True},
      'coverage':{'vertices':72,'edges':204,'scenarios':27,'unordered_pairs':2556,'scenario_pair_entries':69012},
      'failure_rule':'Any false certificate, lost case, unexpected schema/hash drift or algebraic counterexample blocks observed evaluation.'}
    for name,obj in [('computation_manifest.json',manifest),('candidate_manifest.json',candidate)]:
        dest=BASE/name
        if dest.exists():raise RuntimeError('Refusing existing receipt')
        dest.write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({name:sha(BASE/name) for name in ('candidate_manifest.json','computation_manifest.json')}))

if __name__=='__main__':main()
