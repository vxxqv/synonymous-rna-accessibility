"""One-shot root authorizations after the successor freeze is publicly committed."""
from pathlib import Path
import argparse,subprocess,json,hashlib,sys
H=Path(__file__).resolve().parent;P=H/'comt_pilot';I=P/'implementation_v3';O=H/'comt_results/actual_v3';OLD=H/'comt_results/actual_v1';R=H.parents[1]/'outputs/synonymous-rna-accessibility'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
a=argparse.ArgumentParser();a.add_argument('phase',choices=['normalize','predict','evaluate']);a.add_argument('--public-revision',required=True);args=a.parse_args()
freeze=P/'normalization_successor_freeze_v3.json';f=json.loads(freeze.read_text());code=P/'code_freeze_v3.json'
assert sha(code)==f['code_freeze_sha256']
subprocess.run(['git','merge-base','--is-ancestor',args.public_revision,'HEAD'],cwd=R,check=True)
assert hashlib.sha256(subprocess.check_output(['git','show',args.public_revision+':analysis/cycle8/comt_pilot/normalization_successor_freeze_v3.json'],cwd=R)).hexdigest()==sha(freeze)
roles={'protocol':P/'protocol_v2.md','amendment':P/'prospective_reconstruction_amendment.md','integration_addendum':I/'integration_addendum.md','contract':I/'selective_contract.json','code_manifest':code,'normalization_amendment':P/'normalization_contract_amendment_v3.md','predecessor_code_manifest':P/'code_freeze_v1.json'}
roles.update({'normalize':{'extract_manifest':OLD/'extract/sealed_manifest.json','rf_manifest':OLD/'rf/sealed_manifest.json'},'predict':{'normalization_manifest':O/'normalize/sealed_manifest.json'},'evaluate':{'normalization_manifest':O/'normalize/sealed_manifest.json','prediction_manifest':O/'predict/sealed_manifest.json'}}[args.phase])
for k,h in f['input_pins'].items():assert sha(roles[k])==h
if args.phase=='normalize':
 for k in ['extract_manifest','rf_manifest']:assert sha(roles[k])==f[k+'_sha256']
O.mkdir(exist_ok=True);auth=O/f'authorization_{args.phase}.json';out=O/args.phase
if auth.exists() or out.exists():raise FileExistsError('Preserve original successor execution; no overwrite')
d={'schema':'COMT_A1_authorization_v1','phase':args.phase,'status':'approved_frozen_normalization_successor','historical_reproduction_claim':False,'outcome_access_authorized':True,'public_freeze_revision':args.public_revision,'successor_freeze_sha256':sha(freeze),'scope':f['scientific_delta'],'inputs':{k:{'path':str(v.resolve()),'sha256':sha(v)} for k,v in roles.items()}}
if args.phase=='evaluate':d['sealed_prediction_sha256']=sha(roles['prediction_manifest'])
write(auth,d);log=O/f'{args.phase}_execution.log'
with log.open('w',encoding='utf8') as stream:r=subprocess.run([sys.executable,'-B',str(I/'driver.py'),'--phase',args.phase,'--authorization',str(auth),'--out',str(out)],stdout=stream,stderr=subprocess.STDOUT)
rec={'phase':args.phase,'returncode':r.returncode,'authorization_sha256':sha(auth),'log_sha256':sha(log),'sealed_manifest_exists':(out/'sealed_manifest.json').exists(),'public_freeze_revision':args.public_revision}
if rec['sealed_manifest_exists']:rec['sealed_manifest_sha256']=sha(out/'sealed_manifest.json')
write(O/f'{args.phase}_execution_receipt.json',rec);print(json.dumps(rec))
if r.returncode:print(log.read_text()[-5000:])
raise SystemExit(r.returncode)
