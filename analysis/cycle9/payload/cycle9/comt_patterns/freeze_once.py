"""Seal approved protocol and tested implementation before observed comparisons."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

out=HERE/'execution_freeze.json'
if out.exists():raise FileExistsError(out)
if (HERE/'observed_v1').exists():raise RuntimeError('Observed output already exists')
receipt=json.loads((HERE/'synthetic_receipt.json').read_text())
if receipt['status']!='PASS' or receipt['tests']!=18:raise ValueError('Synthetic acceptance failed')
for n,h in receipt['source_hashes'].items():
    if sha(HERE/n)!=h:raise ValueError('Source changed after synthetic checks')
files=['analyze.py','verify.py','test_synthetic.py','synthetic_receipt.json','proposed_protocol.md','closest_prior.md',
       'metadata_inventory.py','metadata_inventory.json','freeze_once.py']
inputs={'input':ROOT/'work/research_cycle8/comt_results/actual_v3/normalize/contrasts.json',
        'inventory':HERE/'metadata_inventory.json',
        'contract':ROOT/'work/research_cycle8/comt_pilot/implementation_v3/selective_contract.json'}
pins={k:sha(p) for k,p in inputs.items()}
if pins['input']!='5654df32df45840415e483a53127a12744887ac9c6b49a18f0a9387bfa04f083':raise ValueError('Outcome source drift')
if sha(HERE/'proposed_protocol.md')!='6ebf332b56642766f3b78345e0c3a9cc85d995628be98272f0a3aee8efce594d':raise ValueError('Approved protocol drift')
freeze={'status':'FROZEN_SYNTHETIC_PASS_ROOT_AUTHORIZED_EXPLORATORY_RUN',
        'approval_basis':'Parent explicitly approved bounded implementation and once-only complete134pair run, preserving all128jointcases and frozen prediction failures.',
        'timestamp_utc':datetime.now(timezone.utc).isoformat(),
        'repository_revision_observed':'057546edc0a8263d348b20440ba51705cfffe361',
        'owned_files':{n:sha(HERE/n) for n in files},'input_hashes':pins,
        'input_paths':{k:p.relative_to(ROOT).as_posix() for k,p in inputs.items()},
        'expected_pairs':134,'expected_joint_cases':17152,'expected_SNV_pairs':97,
        'python':platform.python_version(),'platform':platform.platform(),'machine':platform.machine(),
        'observed_comparisons_computed_before_freeze':False,'prior_exposure':'Cycle8dataalreadyexposed;followupexploratory,notprospectivevalidation.'}
out.write_text(json.dumps(freeze,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'freeze_sha256':sha(out),'source_sha256':freeze['owned_files']['analyze.py'],'expected_joint_cases':17152}))
