"""Builds trace-only authorization for existing wholly invented RF fixtures."""
import copy,json,shutil,subprocess,sys
from pathlib import Path
from replay import sha,dump,execute
P=Path(__file__).resolve().parent
PILOT=P.parent/'comt_pilot'
SOURCE=PILOT/'implementation_v2'/'synthetic_RF_v1'
OUT=P/'invented_v1'
if OUT.exists():raise FileExistsError('Fresh invented test root required')
OUT.mkdir()
for folder in ('prepared','RF'):shutil.copytree(SOURCE/folder,OUT/folder)
common={'protocol':PILOT/'protocol_v2.md','amendment':PILOT/'prospective_reconstruction_amendment.md','integration_addendum':PILOT/'implementation_v2/integration_addendum.md','contract':PILOT/'implementation_v2/selective_contract.json','code_manifest':OUT/'synthetic_code.json'}
dump(common['code_manifest'],{'scope':'Explicitly invented trace fixture, not an actual phase freeze','verifier_python_sha256':sha(P/'replay.py'),'verifier_r_sha256':sha(P/'replay_forest.R')})
pins={n:sha(p) for n,p in common.items()};pins['runtime_manifest']=sha(PILOT/'runtime/runtime_manifest.json')
def seal(folder,phase,more=None):
    d={'schema':'COMT_A1_stage_seal_v1','phase':phase,'input_hashes':{**pins,**(more or {})},'files':{p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.name!='sealed_manifest.json'}}
    dump(folder/'sealed_manifest.json',d)
seal(OUT/'prepared','extract')
seal(OUT/'RF','rf',{'extract_manifest':sha(OUT/'prepared/sealed_manifest.json')})
roles={**common,'runtime_manifest':PILOT/'runtime/runtime_manifest.json','extract_manifest':OUT/'prepared/sealed_manifest.json','rf_manifest':OUT/'RF/sealed_manifest.json','verifier_r':P/'replay_forest.R'}
a={'schema':'COMT_RF_trace_authorization_v1','status':'approved_fixed_trace_only','mode':'synthetic','actual_RF_result_access_authorized':False,'pilot_root':str(PILOT.resolve()),'verifier_python_sha256':sha(P/'replay.py'),'inputs':{n:{'path':str(p.resolve()),'sha256':sha(p)} for n,p in roles.items()}}
dump(OUT/'authorization.json',a)
try:r=execute(OUT/'authorization.json',OUT/'trace_baseline')
except Exception:
    raise # Failed fixture remains without fallback or retraining.
if r['status']!='PASS_SAVED_RF_TRACE':raise AssertionError('Invented trace failed')
run=subprocess.run([sys.executable,'-O',str(P/'replay.py'),'--help'],capture_output=True,text=True)
if run.returncode==0 or 'Optimized Python' not in run.stderr:raise AssertionError('Optimization guard failed')
bad=copy.deepcopy(a);bad['mode']='actual'
dump(OUT/'closed_actual_authorization.json',bad)
try:execute(OUT/'closed_actual_authorization.json',OUT/'forbidden_actual')
except ValueError:pass
else:raise AssertionError('Closed actual access accepted')
if (OUT/'forbidden_actual').exists():raise AssertionError('Closed actual made output')
dump(OUT/'receipt.json',{'status':'PASS_INVENTED_TRACE_AND_GUARDS','actual_source_rows':0,'actual_simulation_rows':0,'saved_forests_retrained':False,'baseline_verification_sha256':sha(OUT/'trace_baseline/verification.json'),'summary':(OUT/'trace_baseline/trace/summary.tsv').read_text(),'source_sha256':sha(Path(__file__))})
print((OUT/'receipt.json').read_text())
