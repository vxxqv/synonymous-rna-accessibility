"""Test the portable adapter in an owned, relocated minimal archive subset."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import subprocess
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
REPO=ROOT/'outputs/synonymous-rna-accessibility'
SCRIPT=BASE/'verify_portable.py'

def need(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(raw):return hashlib.sha256(raw).hexdigest()

def main():
    fixture=BASE/'portable_fixture';need(not fixture.exists(),'Fixture already exists')
    fixture.mkdir()
    freeze_path=REPO/'analysis/cycle7/freeze.json';freeze=json.loads(freeze_path.read_bytes())
    relative=['analysis/cycle7/freeze.json']+['analysis/cycle7/'+name for name in freeze['files']]
    relative+=['analysis/cycle6/fas_pilot/ordinal_geometry.py','analysis/cycle6/fas_pilot/source_contract.json',
               'results/cycle6/full_result_v1.json.gz']
    copied={}
    for name in relative:
        raw=(REPO/name).read_bytes();dest=fixture/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
        copied[name]=sha(raw)
    raw_result=(ROOT/'work/research_cycle7/results/robustness.json').read_bytes()
    need(sha(raw_result)=='861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa','Wrong successor source')
    compressed=gzip.compress(raw_result,mtime=0)
    successor=fixture/'results/cycle7/robustness.json.gz';successor.parent.mkdir(parents=True,exist_ok=True);successor.write_bytes(compressed)
    copied['results/cycle7/robustness.json.gz']=sha(compressed)
    local_script=fixture/'analysis/cycle7/robustness_result_review/verify_portable.py'
    local_script.parent.mkdir(parents=True);local_script.write_bytes(SCRIPT.read_bytes())
    env=os.environ.copy();env['PATH']='';env['PYTHONDONTWRITEBYTECODE']='1'
    def invoke(args,flags=()):
        return subprocess.run([sys.executable,*flags,str(local_script),*map(str,args)],cwd=BASE,env=env,capture_output=True,text=True)
    first=BASE/'portable_default_output'
    p=invoke(['--root',fixture,'--out',first]);need(p.returncode==0,'Default compressed replay failed: '+p.stderr)
    reference=json.loads((BASE/'observed_v1/receipt.json').read_bytes())
    actual=json.loads((first/'receipt.json').read_bytes())
    ignored=('input_hashes_before_and_after','verifier_sha256','started_at','runtime_seconds')
    def payload(x):return {k:v for k,v in x.items() if k not in ignored}
    need(payload(reference)==payload(actual),'Portable scientific payload differs')
    second=BASE/'portable_raw_override_output'
    p=invoke(['--repo',fixture,'--result',ROOT/'work/research_cycle7/results/robustness.json','--out',second])
    need(p.returncode==0,'Raw override failed: '+p.stderr)
    need(payload(json.loads((second/'receipt.json').read_bytes()))==payload(reference),'Raw override payload differs')
    tests=['default_gzip_root_alias','explicit_raw_repo_alias']
    for flag in ('-O','-OO'):
        dest=BASE/('rejected_'+flag.lstrip('-'))
        p=invoke(['--root',fixture,'--out',dest],flags=(flag,))
        need(p.returncode!=0 and not dest.exists() and 'Optimized Python' in p.stderr,'Optimization guard')
        tests.append('reject_'+flag)
    p=invoke(['--root',BASE/'missing_archive','--out',first])
    need(p.returncode!=0 and 'Refusing existing output' in p.stderr,'Overwrite-before-input guard')
    tests.append('reject_existing_output_before_input')
    try:
        successor.write_bytes(gzip.compress(b'{}',mtime=0))
        dest=BASE/'rejected_decoded_hash'
        p=invoke(['--root',fixture,'--out',dest])
        need(p.returncode!=0 and not dest.exists() and 'Decoded successor result hash mismatch' in p.stderr,'Wrong decoded result accepted')
        tests.append('reject_wrong_decoded_successor_before_output')
    finally:successor.write_bytes(compressed)
    need(not (fixture/'.git').exists(),'Unexpected Git metadata')
    for name,digest in copied.items():need(sha((fixture/name).read_bytes())==digest,'Fixture byte drift')
    receipt={'status':'PASS','portable_sha256':sha(SCRIPT.read_bytes()),'original_sha256':sha((BASE/'verify_result.py').read_bytes()),
             'tests':tests,'all_scientific_output_fields_equal':True,
             'ignored_runtime_and_location_metadata_fields':list(ignored),
             'radii':actual['radii_checked'],'row_deletions':actual['row_deletion_certificates_checked'],
             'tested_event_certificates':actual['searched_event_certificates_checked'],
             'git_metadata_absent':True,'empty_PATH':True,'network_required':False,'scientific_packages_imported':False,
             'fixture_input_hashes':copied,'portable_output_sha256':sha((first/'receipt.json').read_bytes()),
             'test_source_sha256':sha(Path(__file__).read_bytes()),
             'scope':'Replays the same frozen saved certificates. No reoptimization, new outcome, or altered scientific procedure.'}
    dest=BASE/'portable_receipt.json';need(not dest.exists(),'Existing test receipt')
    dest.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','tests':len(tests),'receipt_sha256':sha(dest.read_bytes())}))

if __name__=='__main__':main()
