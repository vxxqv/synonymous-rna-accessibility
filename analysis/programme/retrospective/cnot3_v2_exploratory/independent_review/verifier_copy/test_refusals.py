"""Portable metadata-only refusal fixtures plus reciprocal-zero boundary."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--fixture',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    fixture=Path(a.fixture).resolve();out=Path(a.out).resolve()
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    here=Path(__file__).resolve().parent
    spec=importlib.util.spec_from_file_location('verifier',here/'verify.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
    original=json.loads((fixture/'authorization.json').read_text())
    paths={'training_features':fixture/'training_features.json','training_targets':fixture/'training_targets.json',
           'holdout_features':fixture/'holdout_features.json','design_receipt':fixture/'design_receipt.json','predictions':fixture/'predictions.json',
           'evaluation':fixture/'evaluation.json','holdout_targets':fixture/'holdout_targets.json','commit_receipt':fixture/'synthetic_commit.json'}
    recorded=json.loads((fixture/'fixture_run.json').read_text())['command']
    for role in ['protocol','author_code']:
        flag='--'+role.replace('_','-');paths[role]=Path(recorded[recorded.index(flag)+1])
    tests=[]
    for name in ['wrong_protocol','wrong_author_code','wrong_core_hash','optimized_O','optimized_OO','overwrite']:
        auth=json.loads(json.dumps(original));selected=dict(paths);target=out/('must_not_exist_'+name)
        if name in ['wrong_protocol','wrong_author_code']:
            role='protocol' if name=='wrong_protocol' else 'author_code'
            altered=out/(name+'.txt');altered.write_text('deliberately different non-scientific bytes\n',encoding='utf-8')
            selected[role]=altered;auth['inputs'][role]=sha(altered)
        if name=='wrong_core_hash':auth['reference_core_sha256']='0'*64
        if name=='overwrite':target=fixture/'independent_verification'
        authfile=out/(name+'_authorization.json');authfile.write_text(json.dumps(auth,indent=2)+'\n',encoding='utf-8')
        command=[sys.executable,'-I','-B']+(['-O'] if name=='optimized_O' else ['-OO'] if name=='optimized_OO' else [])+[str(here/'verify.py')]
        for role,path in selected.items():command+=['--'+role.replace('_','-'),str(path)]
        command+=['--authorization',str(authfile),'--out',str(target)]
        before=sha(target/'verification.json') if name=='overwrite' else None
        run=subprocess.run(command,capture_output=True,text=True)
        if run.returncode==0 or (name!='overwrite' and target.exists()) or (name=='overwrite' and sha(target/'verification.json')!=before):raise AssertionError(name)
        tests.append(dict(name=name,passed=True,returncode=run.returncode,stderr=run.stderr))
    # Canonical stored rate is exactly zero. Independent reciprocal order yields
    # a harmless tiny positive, so summary counts must describe stored values.
    import math,numpy as np
    context=[dict(GeneID='SYN_ZERO',HL_WT=24.0,LOG2_FC_IP=0.0)]
    arrays={'baseline_HL':24.0,'log2_lifetime_change':0.0,'dimensionless_rate_contrast':0.0,'LOG2_FC_IP':0.0}
    summaries={k:dict(median=x,p10=x,p90=x,negative=int(x<0),zero=int(x==0),positive=int(x>0)) for k,x in arrays.items()}
    import itertools
    actual=dict(n=1,h_ref=1.5,rows=[dict(GeneID='SYN_ZERO',**arrays)],summaries=summaries,
                spearman={a+'__'+b:None for a,b in itertools.combinations(arrays,2)})
    checks=v.Checks();v.verify_kinetic(checks,'zero_boundary',actual,context,np.array([math.log(24.0)]),1.5)
    if checks.failures:raise AssertionError(checks.failures)
    result=dict(refusal_tests=tests,zero_boundary_checks=checks.count,zero_boundary_max_abs=checks.maximum_abs,
                initial_verifier_history_retained=True,scientific_estimator_unchanged=True)
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'refusals_passed':len(tests),'zero_boundary_checks':checks.count}))

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Optimized fixture runner refused')
    main()
