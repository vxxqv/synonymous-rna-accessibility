"""Exercise offline replay and fail-closed archive guards in the owned fixture."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

BASE=Path(__file__).resolve().parent
REPO=BASE/'r'
SCRIPT=REPO/'analysis/cycle6/result_reviews/replay.py'

def sha(raw):return hashlib.sha256(raw).hexdigest()
def need(ok,message):
    if not ok:raise RuntimeError(message)

def invoke(output, flags=(), repo=REPO):
    env=os.environ.copy()
    env['PATH']=''
    env['PYTHONDONTWRITEBYTECODE']='1'
    return subprocess.run([sys.executable,*flags,str(SCRIPT),'--repo',str(repo),'--out',str(output)],
                          cwd=BASE,env=env,capture_output=True,text=True)

def main():
    receipt_path=BASE/'portability_receipt.json'
    need(not receipt_path.exists(),'Refusing existing portability receipt')
    need(not (REPO/'.git').exists(),'Relocated fixture has Git metadata')
    checks=[]
    for flag in ('-O','-OO'):
        out=BASE/('reject_'+flag.lstrip('-'))
        p=invoke(out,(flag,),BASE/'nonexistent_repository')
        need(p.returncode!=0 and 'Optimized Python is not allowed' in p.stderr and not out.exists(),'Optimization guard')
        checks.append({'test':'optimization_'+flag,'status':'PASS','no_output_created':True})
    p=invoke(BASE/'replayed',repo=BASE/'nonexistent_repository')
    need(p.returncode!=0 and 'Refusing existing output directory' in p.stderr,'Overwrite before input guard')
    checks.append({'test':'overwrite_before_inputs','status':'PASS'})
    manifest=REPO/'analysis/cycle6/result_reviews/replay_manifest.json'
    raw=manifest.read_bytes()
    try:
        manifest.write_bytes(raw+b' ')
        out=BASE/'reject_manifest'
        p=invoke(out)
        need(p.returncode!=0 and 'Pinned replay manifest mismatch' in p.stderr and not out.exists(),'Manifest tamper guard')
        checks.append({'test':'manifest_tamper_before_output','status':'PASS'})
    finally:
        manifest.write_bytes(raw)
    cfg=json.loads(raw)
    selected=[
        next(x for x in cfg['files'] if x['path'].startswith('data/')),
        next(x for x in cfg['files'] if x['path'].endswith('verify_counts_offline.py')),
        next(x for x in cfg['files'] if x['path'].endswith('.gz')),
    ]
    for index,item in enumerate(selected):
        path=REPO/item['path'];original=path.read_bytes()
        try:
            path.write_bytes(original+b' ')
            out=BASE/f'reject_input_{index}'
            p=invoke(out)
            need(p.returncode!=0 and 'Pinned input mismatch' in p.stderr and not out.exists(),'Input tamper guard')
            checks.append({'test':'input_tamper_before_output','input':item['path'],'status':'PASS'})
        finally:
            path.write_bytes(original)
    path=REPO/selected[0]['path'];moved=path.with_suffix(path.suffix+'.temporarily_absent')
    need(not moved.exists(),'Missing-input scratch collision')
    try:
        path.rename(moved)
        out=BASE/'reject_missing'
        p=invoke(out)
        need(p.returncode!=0 and 'FileNotFoundError' in p.stderr and not out.exists(),'Missing input guard')
        checks.append({'test':'missing_input_before_output','status':'PASS'})
    finally:
        moved.rename(path)
    for item in cfg['files']:
        need(sha((REPO/item['path']).read_bytes())==item['sha256'],'Fixture restoration mismatch')
    out=BASE/'replayed_empty_path'
    p=invoke(out)
    need(p.returncode==0,'Offline replay failed: '+p.stdout+p.stderr)
    audit=json.loads((out/'replay_receipt.json').read_bytes())
    need(audit['status']=='PASS' and audit['certificates_checked']==30 and audit['pair_records_checked']==74124,'Incomplete replay')
    checks.append({'test':'relocated_no_git_metadata_empty_PATH_full_replay','status':'PASS',
                   'certificates':30,'pair_records':74124,'selected_field_comparisons':3816})
    result={'status':'PASS','scope':'Portability and guard verification only; original scientific algorithms unchanged.',
            'checks':checks,'git_metadata_absent':True,'process_PATH_empty':True,
            'no_network_code_or_requests':'Source inspection: neither replay path imports network clients or makes requests.',
            'wrapper_sha256':sha(SCRIPT.read_bytes()),'manifest_sha256':sha(raw),
            'replay_receipt_sha256':sha((out/'replay_receipt.json').read_bytes()),
            'preparation_receipt_sha256':sha((BASE/'preparation_receipt.json').read_bytes()),
            'verifier_sha256':sha(Path(__file__).read_bytes())}
    receipt_path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','guard_and_replay_checks':len(checks),'receipt_sha256':sha(receipt_path.read_bytes())}))

if __name__=='__main__':main()
