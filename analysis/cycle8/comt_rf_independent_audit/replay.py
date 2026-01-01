"""CLI-only saved-RF trace verifier. Does not train or read archive bodies."""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Optimized Python refused before inputs')
P=Path(__file__).resolve().parent
RUNTIME='a97c2057080245e691640ca955af6248aa93102b58341b76f37f473461231c95'
COMMON=('protocol','amendment','integration_addendum','contract','code_manifest')
FILES={'extract':('rows_1.tsv','rows_2.tsv'), 'rf':('model_1.rds','model_2.rds','predictions_1.tsv','predictions_2.tsv','row_status_1.tsv','row_status_2.tsv','rejections.tsv','status.txt')}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')
def extended(p):
    s=str(p.resolve())
    if os.name=='nt' and not s.startswith('\\\\?\\'):s='\\\\?\\UNC\\'+s[2:] if s.startswith('\\\\') else '\\\\?\\'+s
    return Path(s)
def execute(authfile,out):
    authfile=Path(authfile);out=Path(out)
    if out.exists():raise FileExistsError('Fresh output required before input reads')
    auth=json.loads(authfile.read_text())
    if auth.get('schema')!='COMT_RF_trace_authorization_v1' or auth.get('status')!='approved_fixed_trace_only':raise ValueError('Missing explicit trace authorization')
    if auth.get('mode') not in ('synthetic','actual'):raise ValueError('Declared mode required')
    if auth.get('mode')=='actual' and auth.get('actual_RF_result_access_authorized') is not True:raise ValueError('Actual access closed')
    roles={}
    for n,x in auth['inputs'].items():
        path=Path(x['path']).resolve()
        if sha(path)!=x['sha256']:raise ValueError('Input hash mismatch '+n)
        roles[n]=path
    if not set(COMMON+('extract_manifest','rf_manifest','runtime_manifest','verifier_r'))<=set(roles):raise ValueError('Missing trace pins')
    if sha(roles['runtime_manifest'])!=RUNTIME:raise ValueError('Changed runtime specification')
    if roles['verifier_r']!=P/'replay_forest.R':raise ValueError('Unexpected numerical verifier')
    if auth['verifier_python_sha256']!=sha(Path(__file__)):raise ValueError('Verifier Python changed')
    stages={}
    for stage in ('extract','rf'):
        path=roles[stage+'_manifest'];m=json.loads(path.read_text())
        if m.get('schema')!='COMT_A1_stage_seal_v1' or m.get('phase')!=stage:raise ValueError('Wrong stage')
        for n in COMMON:
            if m['input_hashes'].get(n)!=sha(roles[n]):raise ValueError('Changed common stage pin '+n)
        if m['input_hashes'].get('runtime_manifest')!=RUNTIME:raise ValueError('Stage runtime differs')
        if not set(FILES[stage])<=set(m['files']):raise ValueError('Unsealed required payload')
        for n,h in m['files'].items():
            if Path(n).name!=n or sha(path.parent/n)!=h:raise ValueError('Changed sealed stage payload '+n)
        stages[stage]=m
    if stages['rf']['input_hashes'].get('extract_manifest')!=sha(roles['extract_manifest']):raise ValueError('RF/extraction association differs')
    pilot=Path(auth['pilot_root']).resolve();runtime=json.loads(roles['runtime_manifest'].read_text())
    if sha(Path(runtime['Rscript']))!=runtime['Rscript_sha256']:raise ValueError('R binary changed')
    for n,h in runtime['package_files'].items():
        if sha(extended(pilot/'runtime/library'/n))!=h:raise ValueError('R package changed '+n)
    out.mkdir(parents=True)
    run=subprocess.run([runtime['Rscript'],'--vanilla',str(roles['verifier_r']),str(pilot),str(roles['extract_manifest'].parent),str(roles['rf_manifest'].parent),str(out/'trace')],capture_output=True,text=True,timeout=600)
    (out/'execution.log').write_text(run.stdout+'\n'+run.stderr,encoding='utf8')
    receipt={'status':'PASS_SAVED_RF_TRACE' if run.returncode==0 else 'FAILED_TRACE_PRESERVED','authorization_sha256':sha(authfile),'mode':auth['mode'],'input_hashes':{n:sha(p) for n,p in roles.items()},'returncode':run.returncode,'python':sys.version,'probability_absolute_tolerance':1e-14,'claim':'Replay of saved states and source masks, not independent training, historical seed authentication or biological calibration.','outputs':{str(f.relative_to(out)).replace('\\','/'):sha(f) for f in out.rglob('*') if f.is_file()}}
    dump(out/'verification.json',receipt)
    if run.returncode:raise RuntimeError('RF trace failed; retained outputs')
    return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--authorization',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();print(json.dumps(execute(a.authorization,a.output),indent=2))
