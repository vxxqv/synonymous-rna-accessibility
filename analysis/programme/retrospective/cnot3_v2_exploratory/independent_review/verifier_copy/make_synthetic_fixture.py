"""Build the already reviewed artificial320-row fixture and candidate phase files.

The author is imported only to CREATE synthetic comparison artifacts. verify.py
never imports it and independently recomputes all numerical references.
"""
import sys
if not __debug__:
    raise RuntimeError("Optimized Python refused")
sys.dont_write_bytecode=True
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def import_path(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def write(p,obj):
    with Path(p).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(obj,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--author-code',required=True);p.add_argument('--protocol',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out).resolve()
    if out.exists():raise FileExistsError('Fixture directory must be fresh')
    verifier=import_path('verifier_configuration',HERE/'verify.py')
    author_path,protocol=Path(args.author_code).resolve(),Path(args.protocol).resolve()
    if sha(author_path)!=verifier.AUTHOR_SHA or sha(protocol)!=verifier.PROTOCOL_SHA or sha(HERE/'reference_core.py')!=verifier.CORE_SHA:
        raise ValueError('Frozen candidate/protocol/core mismatch')
    core=import_path('independent_fixture_core',HERE/'reference_core.py')
    author=import_path('synthetic_fixture_author',author_path)
    train,test,truth=core.create_data()
    training={r['GeneID']:truth[r['GeneID']] for r in train}; heldout={r['GeneID']:truth[r['GeneID']] for r in test}
    design={'status':'synthetic_only'}
    out.mkdir(parents=True)
    paths={'training_features':out/'training_features.json','training_targets':out/'training_targets.json',
           'holdout_features':out/'holdout_features.json','design_receipt':out/'design_receipt.json','predictions':out/'predictions.json',
           'evaluation':out/'evaluation.json','holdout_targets':out/'holdout_targets.json','commit_receipt':out/'synthetic_commit.json',
           'protocol':protocol,'author_code':author_path}
    for role,value in [('training_features',train),('training_targets',training),('holdout_features',test),('design_receipt',design),('holdout_targets',heldout)]:write(paths[role],value)
    artifact=author.predict(train,training,test,design)
    artifact['file_provenance']={'implementation_sha256':verifier.AUTHOR_SHA,'freeze_sha256':'synthetic fixture, no research commit',
                                 'inputs':{k:sha(paths[k]) for k in ['training_features','training_targets','holdout_features','design_receipt','protocol']}}
    write(paths['predictions'],artifact)
    commit={'status':'synthetic_only','predictions_sha256':sha(paths['predictions']),'holdout_targets_sha256':sha(paths['holdout_targets'])}
    write(paths['commit_receipt'],commit)
    evaluation=author.evaluate(artifact,heldout)
    evaluation['provenance']={'predictions_sha256':sha(paths['predictions']),'holdout_targets_sha256':sha(paths['holdout_targets']),
                              'commit_receipt_sha256':sha(paths['commit_receipt']),'implementation_sha256':verifier.AUTHOR_SHA}
    write(paths['evaluation'],evaluation)
    auth={'status':'synthetic_only','verifier_sha256':sha(HERE/'verify.py'),'reference_core_sha256':verifier.CORE_SHA,
          'inputs':{k:sha(v) for k,v in paths.items()}}
    write(out/'authorization.json',auth)
    command=[sys.executable,'-I','-B',str(HERE/'verify.py')]
    for k,v in paths.items():command+=['--'+k.replace('_','-'),str(v)]
    command+=['--authorization',str(out/'authorization.json'),'--out',str(out/'independent_verification')]
    run=subprocess.run(command,capture_output=True,text=True)
    write(out/'fixture_run.json',{'status':'passed' if run.returncode==0 else 'failed','synthetic_only':True,'command':command,
                                 'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                                 'seed':912731,'total_genes':len(train)+len(test),'author_used_for_fixture_generation_only':True})
    print(run.stdout,end='')
    if run.returncode:
        print(run.stderr,file=sys.stderr)
        raise SystemExit(run.returncode)

if __name__=='__main__':main()
