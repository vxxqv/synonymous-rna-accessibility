"""Replay only frozen sequence features. No response access, fitting or Git."""
if not __debug__:
 raise RuntimeError('Scientific assertions are required; Python -O and -OO are refused')

from pathlib import Path
import argparse, hashlib, importlib.util, json, math, os, platform, subprocess, sys, time, traceback, warnings
sys.dont_write_bytecode=True
warnings.simplefilter('error')

def absolute(p):
 p=Path(p).resolve()
 return Path('\\\\?\\'+str(p)) if os.name=='nt' and not str(p).startswith('\\\\?\\') else p
HERE=absolute(__file__).parent
TOLERANCE=1e-10
FIELDS=('SequenceID','RNA.sequence','CDS.sequence.longest','Barcode')
PROTOCOL='analysis/programme/retrospective/persist_v2/protocol.json'
AUTHOR='analysis/programme/retrospective/persist_v2/persist_prediction.py'
PARAMETERS='analysis/programme/retrospective/persist_v2/Turner2004.par'
EXPECTED='results/programme/persist_v2/train/features.json'
METADATA='results/programme/persist_v2/feature_replay/selected_metadata.json'
def require(ok,message):
 if not ok:raise RuntimeError(message)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def text_sha(s):return hashlib.sha256(s.encode('ascii')).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write_new(p,value):
 with Path(p).open('x',encoding='utf-8',newline='\n') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def save_warning_checked_row(out,result):
 write_new(out/'row.json',result)
 require(not result['python_warnings'],'Worker Python warning; row retained and execution stopped')
def accept_process(proc,key):
 require(not proc.stdout and not proc.stderr,'Worker stdout/stderr is nonempty; stop without retry: '+key)
 require(proc.returncode==0,'Worker failed with no retry/fallback: '+key)
def module_at(path):
 spec=importlib.util.spec_from_file_location('preserved_persist_prediction',path)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def versions():
 require(platform.python_implementation()=='CPython' and sys.version_info[:2]==(3,12),'CPython3.12 required; no installation is performed')
 import numpy as np
 import RNA
 require(np.__version__=='2.5.4','NumPy2.5.4 required')
 require(RNA.__version__=='2.7.2' and RNA.pf_float_precision()==0,'ViennaRNA2.7.2 double-precision PF required')
 return np,RNA
def identity_check(rows,protocol):
 ids=protocol['all_ids'];require(len(ids)==37 and len(set(ids))==37,'Expected37 unique frozen IDs')
 require(isinstance(rows,list) and [r['SequenceID'] for r in rows]==ids,'Selected metadata order/membership mismatch')
 cohort={r['SequenceID']:r for r in protocol['cohort']};require(set(cohort)==set(ids),'Protocol cohort mismatch')
 for row in rows:
  key=row['SequenceID'];entry=cohort[key];require(set(row)==set(FIELDS),'Unexpected selected metadata field')
  for field,digest,length in [('RNA.sequence','deposited_full_sha256','full_length'),('CDS.sequence.longest','CDS_sha256','CDS_length'),('Barcode','barcode_sha256','barcode_length')]:
   require(type(row[field]) is str and len(row[field])==entry[length] and text_sha(row[field])==entry[digest],'Original protocol sequence hash/length mismatch: '+key+'/'+field)
  require(text_sha(row['RNA.sequence'].replace('T','U'))==entry['normalized_full_sha256'],'Normalized sequence hash mismatch: '+key)
def authenticated(repo):
 manifest=load(HERE/'manifest.json')
 for rel,h in manifest['code_sha256'].items():require(sha(HERE/rel)==h,'Replay code hash mismatch: '+rel)
 for rel,h in manifest['repository_inputs'].items():require(sha(repo/rel)==h,'Archived input hash mismatch: '+rel)
 protocol=load(repo/PROTOCOL);rows=load(repo/METADATA);identity_check(rows,protocol)
 return manifest,protocol,rows
def worker(repo,out,key):
 require(out.is_dir() and not any(out.iterdir()),'Internal worker output must be an empty new directory')
 require(repo!=out and repo not in out.parents,'Worker output must be outside archived input repository')
 manifest,protocol,rows=authenticated(repo);np,RNA=versions();row=next((r for r in rows if r['SequenceID']==key),None);require(row is not None,'Unselected worker identifier')
 module=module_at(repo/AUTHOR);require(len(module.COMPOSITION_NAMES+module.CODON_NAMES+(module.MUP_NAME,))==83,'Unexpected feature definition')
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  # This configuration-only call records actual MD and exported global parameters.
  # fold_one below still executes its own unchanged configured_md call and single rescaling/PF calculation.
  loaded,md=module.configured_md(protocol['vienna'])
  snapshot=module.snapshot_md(md,loaded)
  require(snapshot==protocol['vienna']['effective_md'],'Effective MD differs')
  require(bool(loaded.params_save(str(out/'Turner2004.exported.par'))),'Parameter export failed')
  require((out/'Turner2004.exported.par').read_bytes()==(repo/PARAMETERS).read_bytes(),'Exported parameter contents differ')
  diagnostic=module.fold_one(row['RNA.sequence'],protocol['vienna'])
  values=module.sequence_features(row,diagnostic['value'])
 result={'SequenceID':key,'values':values,'diagnostic':diagnostic,'effective_md':snapshot,'parameter_export_sha256':sha(out/'Turner2004.exported.par'),'python_warnings':[{'category':w.category.__name__,'message':str(w.message),'filename':w.filename,'line':w.lineno} for w in caught]}
 save_warning_checked_row(out,result)
def compare_value(observed,expected,path,fail,errors):
 if isinstance(expected,(bool,str,int)) or expected is None:
  if type(observed) is not type(expected) or observed!=expected:fail.append({'path':path,'expected':expected,'actual':observed,'kind':'exact_identity_or_discrete'})
 elif isinstance(expected,float):
  if not isinstance(observed,(float,int)) or isinstance(observed,bool) or not math.isfinite(observed):fail.append({'path':path,'expected':expected,'actual':observed,'kind':'nonfinite_or_type'});return
  error=abs(observed-expected);errors.append((error,path))
  if error>TOLERANCE:fail.append({'path':path,'expected':expected,'actual':observed,'absolute_error':error,'tolerance':TOLERANCE})
 elif isinstance(expected,list):
  if not isinstance(observed,list) or len(observed)!=len(expected):fail.append({'path':path,'kind':'list_shape'});return
  for j,(a,b) in enumerate(zip(observed,expected)):compare_value(a,b,path+'/'+str(j),fail,errors)
 elif isinstance(expected,dict):
  if not isinstance(observed,dict) or set(observed)!=set(expected):fail.append({'path':path,'kind':'object_keys'});return
  for k in expected:compare_value(observed[k],expected[k],path+'/'+k,fail,errors)
 else:raise RuntimeError('Unsupported expected comparison type')
def run(repo,out):
 require(not out.exists(),'Refusing overwrite: '+str(out))
 require(repo!=out and repo not in out.parents,'Output must be outside archived input repository')
 manifest,protocol,rows=authenticated(repo);np,RNA=versions()
 before={rel:sha(repo/rel) for rel in manifest['repository_inputs']};out.mkdir(parents=True,exist_ok=False)
 started=time.monotonic();write_new(out/'STARTED.json',{'status':'started','tolerance_fixed_before_execution':TOLERANCE,'row_count':37,'feature_count':83,'input_sha256':before})
 environment={'platform':platform.platform(),'machine':platform.machine(),'python':platform.python_version(),'python_executable_sha256':sha(sys.executable),'numpy':np.__version__,'ViennaRNA':RNA.__version__,'pf_float_precision_code':RNA.pf_float_precision(),'numpy_binary_sha256':sha(np._core._multiarray_umath.__file__),'ViennaRNA_binary_sha256':sha(RNA._RNA.__file__),'original_Windows_binary_identity_required':False}
 write_new(out/'environment.json',environment)
 completed=[]
 try:
  for row in rows:
   key=row['SequenceID'];destination=out/key;destination.mkdir(exist_ok=False)
   cmd=[sys.executable,'-I','-B',str(HERE/'replay.py'),'--repo-root',str(repo),'--worker',key,'--output',str(destination)]
   proc=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,encoding='utf-8',errors='replace',timeout=180)
   (destination/'stdout.txt').write_text(proc.stdout,encoding='utf-8');(destination/'stderr.txt').write_text(proc.stderr,encoding='utf-8')
   accept_process(proc,key)
   require((destination/'row.json').is_file(),'Missing row output: '+key)
   completed.append(key)
  expected=load(repo/EXPECTED);module=module_at(repo/AUTHOR);names=list(module.COMPOSITION_NAMES+module.CODON_NAMES+(module.MUP_NAME,))
  values={};diagnostics={};warning_records=[];fail=[];errors=[]
  for key in protocol['all_ids']:
   row=load(out/key/'row.json');values[key]=row['values'];diagnostics[key]=row['diagnostic']
   require(len(values[key])==83,'Feature count mismatch: '+key)
   require(row['effective_md']==protocol['vienna']['effective_md'],'Worker MD mismatch: '+key)
   require(row['parameter_export_sha256']==manifest['repository_inputs'][PARAMETERS],'Worker parameter mismatch: '+key)
   warning_records.append({'SequenceID':key,'python_warnings':row['python_warnings'],'stderr':(out/key/'stderr.txt').read_text(encoding='utf-8')})
  observed={'ids':protocol['all_ids'],'names':names,'values':values,'diagnostics':diagnostics,'protocol_sha256':manifest['repository_inputs'][PROTOCOL]}
  compare_value(observed,expected,'features',fail,errors);write_new(out/'features.json',observed)
  write_new(out/'comparison.json',{'absolute_tolerance':TOLERANCE,'float_comparisons':len(errors),'maximum_absolute_error':max((e for e,_ in errors),default=0),'maximum_error_locations':[p for e,p in errors if e==max((x for x,_ in errors),default=0)],'discrepancies':fail,'ids_exact':observed['ids']==expected['ids'],'feature_names_exact':names==expected['names'],'warnings':warning_records})
  after={rel:sha(repo/rel) for rel in before};require(after==before,'Input changed during replay')
  require(not any(x['python_warnings'] or x['stderr'] for x in warning_records),'Warning present; no successful receipt permitted')
  require(not fail,'Feature/diagnostic mismatch; failed evidence retained, no retry')
  write_new(out/'receipt.json',{'status':'passed_all37_rows_all83_features','rows':len(completed),'features_per_row':83,'feature_values':37*83,'float_comparisons':len(errors),'absolute_tolerance':TOLERANCE,'maximum_absolute_error':max((e for e,_ in errors),default=0),'all_original_sequence_hashes_verified':True,'discrete_diagnostics_exact':True,'effective_md_exact':True,'parameter_contents_exact_all37':True,'source_sha256':sha(HERE/'replay.py'),'input_sha256_before':before,'input_sha256_after':after,'inputs_immutable':True,'historical_git_reauthenticated':False,'historical_access_order_reauthenticated':False,'row_processes':37,'worker_retries':0,'fallbacks':0,'exclusions':0,'runtime_seconds':time.monotonic()-started,'environment':environment,'result_sha256':sha(out/'features.json'),'warning_rows':sum(bool(x['python_warnings'] or x['stderr']) for x in warning_records),'scope':'Replay of specified model features only. No target comparisons, fitting, measured-structure validation or biological validity claim.'})
  print(json.dumps({'status':'passed','rows':37,'features_per_row':83,'maximum_absolute_error':max((e for e,_ in errors),default=0),'output':str(out)},indent=2))
 except BaseException as exc:
  write_new(out/'FAILURE.json',{'status':'failed_without_retry','completed_ids':completed,'exception':repr(exc),'traceback':traceback.format_exc(),'inputs_after':{rel:sha(repo/rel) for rel in before}})
  raise
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--worker')
 args=ap.parse_args();repo,out=absolute(args.repo_root),absolute(args.output)
 if args.worker:worker(repo,out,args.worker)
 else:run(repo,out)
if __name__=='__main__':main()
