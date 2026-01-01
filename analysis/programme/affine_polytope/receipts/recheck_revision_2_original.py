"""Replay the corrected corpus using the reviewed exact independent oracle."""
from pathlib import Path
from fractions import Fraction
import hashlib
import importlib.util
import json
from time import perf_counter

BASE = Path(__file__).resolve().parent
OUT = BASE/'revision_2'
sha = lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def encode(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value

receipt=OUT/'exact_recheck.json'
if receipt.exists():raise FileExistsError('Refusing to replace receipt')
ledger=json.loads((OUT/'revision_ledger.json').read_text())
inputs={name:sha(OUT/name) for name in ledger['input_sha256']}
assert inputs==ledger['input_sha256']
oracle_path=BASE/'independent_review/check_exact.py'
assert sha(oracle_path)==ledger['reference_code_sha256']
production=module('production_v2',OUT/'affine_polytope.py')
reference=module('independent_oracle',oracle_path)
start=perf_counter()
fields=['fixed','guaranteed','possible','barrier_minimum','barrier_maximum']
comparisons=0
results=[]
old={r['id']:r['exact'] for r in json.loads((BASE/'recorded_check/results.json').read_text())['cases']}
for case in json.loads((OUT/'cases.json').read_text()):
    answer=production.analyze(case)
    check=reference.check_case(case)
    assert answer['status']==check['status'],case['id']
    if answer['status']=='connected':
        for name in fields:
            assert answer[name]==check[name],(case['id'],name)
            comparisons+=1
    old_name='two_dimensional_corner_failure' if case['id']=='rank_one_sum_two_parameters' else case['id']
    if old_name in old:
        assert answer['status']==old[old_name]['status']
        if answer['status']=='connected':
            assert all(answer[name]==Fraction(old[old_name][name]) for name in fields)
    results.append({'id':case['id'],'exact':answer,'reference':check})
assert inputs=={name:sha(OUT/name) for name in inputs}
result={'status':'passed','cases':len(results),'exact_value_comparisons':comparisons,
        'all_original_values_preserved':True,'production_algorithm_unchanged':True,
        'seconds':perf_counter()-start,'inputs_sha256':inputs,
        'reference_code_sha256':sha(oracle_path),'recheck_code_sha256':sha(Path(__file__)),
        'arithmetic':'exact fractions.Fraction, no tolerance',
        'independence':'Root replay of separately authored oracle; not a new blinded review',
        'cases_detail':results}
receipt.write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='cases_detail'}))
