"""Portable exact replay of affine uncertainty and route-quantifier checks.

The independent oracle uses vertex subsets and Fourier-Motzkin elimination.
It does not use the production score-order arrangement or widest-path routine.
"""
from pathlib import Path
from fractions import Fraction
from time import perf_counter
import argparse
import hashlib
import importlib.util
import json
import platform

HERE=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def encode(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value

def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True,help='New output directory')
    args=parser.parse_args()
    out=args.out.resolve()
    if out.exists():raise FileExistsError('Output directory must be new')
    out.mkdir(parents=True)
    paths=[HERE/name for name in ['claims.md','protocol.json','cases.json','affine_polytope.py',
            'reference/check_exact.py','historical_v1/recorded_check/results.json']]+[Path(__file__)]
    before={p.relative_to(HERE).as_posix():digest(p) for p in paths}
    ledger=json.loads((HERE/'revision_ledger.json').read_text(encoding='utf-8'))
    assert all(before[name]==value for name,value in ledger['input_sha256'].items())
    assert before['reference/check_exact.py']==ledger['reference_code_sha256']
    (out/'inputs.json').write_text(json.dumps(before,indent=2)+'\n',encoding='utf-8')
    production=load('affine_production',HERE/'affine_polytope.py')
    reference=load('affine_independent',HERE/'reference/check_exact.py')
    start=perf_counter()
    cases=json.loads((HERE/'cases.json').read_text(encoding='utf-8'))
    old={r['id']:r['exact'] for r in json.loads((HERE/'historical_v1/recorded_check/results.json').read_text())['cases']}
    fields=['fixed','guaranteed','possible','barrier_minimum','barrier_maximum']
    comparisons=0;results=[];preserved=0
    for case in cases:
        answer=production.analyze(case);check=reference.check_case(case)
        assert answer['status']==check['status'],case['id']
        if answer['status']=='connected':
            for field in fields:
                assert answer[field]==check[field],(case['id'],field)
                comparisons+=1
        old_name='two_dimensional_corner_failure' if case['id']=='rank_one_sum_two_parameters' else case['id']
        if old_name in old:
            assert answer['status']==old[old_name]['status']
            if answer['status']=='connected':
                assert all(answer[k]==Fraction(old[old_name][k]) for k in fields)
            preserved+=1
        results.append({'id':case['id'],'exact':answer,'reference':check})
    assert len(cases)==82 and comparisons==400 and preserved==81
    graphs=reference.graph_suite();formulas=reference.formula_suite()
    after={p.relative_to(HERE).as_posix():digest(p) for p in paths}
    assert before==after,'Sources changed during run'
    result={'status':'passed','cases':len(cases),'exact_value_comparisons':comparisons,
            'original_cases_preserved':preserved,'graphs':graphs,'formulas':formulas,
            'arithmetic':'exact rational; no tolerance','inputs_sha256':before,
            'python':platform.python_version(),'seconds':perf_counter()-start,
            'scope':'Finite corroboration plus retained proofs. No statistical calibration or priority claim.',
            'cases_detail':results}
    (out/'results.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:encode(v) for k,v in result.items() if k not in ['inputs_sha256','cases_detail']}))

if __name__=='__main__':main()
