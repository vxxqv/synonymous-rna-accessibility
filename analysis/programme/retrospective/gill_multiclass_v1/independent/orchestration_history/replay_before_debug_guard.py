"""Offline orchestration of byte-preserved independent scientific implementations."""
from pathlib import Path
import argparse, gzip, hashlib, importlib.util, json, os, platform, shutil, sys, time

sys.dont_write_bytecode=True
def absolute(path):
 p=Path(path).resolve()
 return Path('\\\\?\\'+str(p)) if os.name=='nt' and not str(p).startswith('\\\\?\\') else p
HERE=absolute(__file__).parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
 b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def replace_once(source,old,new):
 require(source.count(old)==1,'Frozen adapter binding not unique: '+old)
 return source.replace(old,new,1)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo-root',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True)
 args=ap.parse_args();repo=absolute(args.repo_root);output=absolute(args.output)
 require(not output.exists(),'Refusing to overwrite existing output: '+str(output))
 require(platform.python_implementation()=='CPython' and platform.python_version()=='3.12.14','Requires CPython3.12.14; no installation is performed')
 import numpy as np
 require(np.__version__=='2.5.4','Requires NumPy2.5.4; no installation is performed')
 package=load(HERE/'package_manifest.json')
 # All preserved source/history bytes and external scientific inputs are authenticated before executing either scientific module.
 for rel,digest in package['package_files'].items():
  require(sha(HERE/rel)==digest,'Package file hash mismatch: '+rel)
 before={}
 for rel,digest in package['external_inputs'].items():
  p=repo/rel;require(p.is_file(),'Missing archived repository input: '+rel)
  before[rel]=sha(p);require(before[rel]==digest,'External input hash mismatch: '+rel)
 require(output!=repo and repo not in output.parents,'Output must be outside the archived input repository')
 output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 write(output/'STARTED.json',{'status':'started','historical_git_reauthenticated':False,'external_input_hashes':before})
 scratch=output/'adapter';scratch.mkdir()
 # Preserve the frozen module __file__ identity used by its own source-hash assertion.
 for name in ['verifier.py','freeze.json','compare_outputs.py']:
  shutil.copyfile(HERE/'frozen'/name,scratch/name)
 spec=importlib.util.spec_from_file_location('verifier',scratch/'verifier.py')
 module=importlib.util.module_from_spec(spec);sys.modules['verifier']=module;spec.loader.exec_module(module)
 module.HERE=scratch;module.ROOT=repo
 module.BASE=repo/'analysis/programme/retrospective/gill_multiclass_v1'
 module.TRAIN=repo/'results/programme/gill_multiclass_v1/train'
 vdir=repo/'results/programme/gill_multiclass_v1/validate'
 evout=scratch/'evaluation_run_v1/calculation'
 opts=argparse.Namespace(validation_receipt=str(vdir/'receipt.json'),validation_receipt_sha=package['external_inputs']['results/programme/gill_multiclass_v1/validate/receipt.json'],validation_counts=str(vdir/'heldout_counts.json.gz'),validation_counts_sha=package['external_inputs']['results/programme/gill_multiclass_v1/validate/heldout_counts.json.gz'],output=str(evout))
 module.evaluate(opts)
 independent_sha=sha(evout/'independent_results.json.gz')
 # Only environment bindings and the current run's output-integrity digest are adapted in memory.
 # The frozen comparator's formulas, coverage, budgets and comparison logic remain byte-preserved on disk.
 src=(scratch/'compare_outputs.py').read_text(encoding='utf-8')
 src=replace_once(src,"H=Path(__file__).resolve().parent; R=H.parents[2]; P=R/'work/gill_multiclass_cycle3'",f"H=Path(__file__).resolve().parent; R=Path({str(repo)!r}); P=R/'results/programme/gill_multiclass_v1'")
 source_pin=package['frozen_independent_result_sha256']
 src=replace_once(src,"assert sha(H/'evaluation_run_v1/calculation/independent_results.json.gz')=='"+source_pin+"'","assert sha(H/'evaluation_run_v1/calculation/independent_results.json.gz')=='"+independent_sha+"'")
 # Exact, enumerated historical path-token remapping. The input hash values are unchanged.
 require(src.count('validation_v1/')==7 and src.count('training_v1/')==3,'Unexpected frozen input path-token count')
 src=src.replace('validation_v1/','validate/').replace('training_v1/','train/')
 binding_record={'adapter':'explicit directory globals, historical input path tokens, current replay output integrity digest only','frozen_verifier_sha256':sha(scratch/'verifier.py'),'frozen_comparison_sha256':sha(scratch/'compare_outputs.py'),'executed_comparison_binding_sha256':hashlib.sha256(src.encode()).hexdigest(),'historical_independent_result_sha256':source_pin,'current_independent_result_sha256':independent_sha,'historical_independent_bytes_match':independent_sha==source_pin,'mathematical_source_changed':False,'git_reauthenticated':False}
 write(output/'binding_record.json',binding_record)
 exec(compile(src,str(scratch/'compare_outputs.py'),'exec'),{'__file__':str(scratch/'compare_outputs.py'),'__name__':'__portable_comparison__'})
 compared=load(scratch/'comparison_v3.json')
 require(compared['status']=='pass' and not compared['discrepancies'],'Scientific comparison failed; all artifacts retained')
 expected=load(HERE/'expected_comparison.json')
 require(compared['checks']==expected['checks'],'Comparison coverage differs')
 require(len(compared['headlines'])==len(expected['headlines']),'Headline coverage differs')
 # Same scientific tolerance as the frozen comparison; byte identity is reported separately.
 for x,y in zip(compared['headlines'],expected['headlines']):
  for key in ['pseudocount','partition','classes','vertices','strata','target_pass']:require(x[key]==y[key],'Headline contract differs: '+key)
  for key in ['residual_cosine','residual_prediction_gain_vs_zero']:
   require(abs(x[key]-y[key])<=1e-9,'Scientific headline differs: '+key)
  require(max(abs(a-b) for a,b in zip(x['leave_one_stratum_range'],y['leave_one_stratum_range']))<=1e-9,'Influence headline differs')
 after={rel:sha(repo/rel) for rel in before};require(after==before,'Input mutated during replay')
 require(sha(evout/'independent_results.json.gz')==independent_sha,'Independent output mutated during comparison')
 # Final replay products live outside the source tree; no source/result overwrite occurs.
 shutil.copyfile(scratch/'comparison_v3.json',output/'comparison.json')
 receipt={'status':'complete_computational_agreement_scientific_target_failed','checks':sum(compared['checks'].values()),'unresolved_discrepancies':0,'historical_git_reauthenticated':False,'historical_git_status':'Historical commit identities are provenance claims copied from the original audit, not independently authenticated in this offline replay. No Git or network call is made.','input_hashes_before':before,'input_hashes_after':after,'inputs_unchanged':True,'python':platform.python_version(),'numpy':np.__version__,'runtime_seconds':time.monotonic()-start,'independent_result_sha256':independent_sha,'comparison_sha256':sha(output/'comparison.json'),'headlines':compared['headlines'],'binding_record':binding_record,'limitations':'Numerical replay of archived exported counts; no workbook/FASTQ extraction, biological calibration or historical access-order reauthentication.'}
 write(output/'replay_receipt.json',receipt)
 print(json.dumps({'status':receipt['status'],'checks':receipt['checks'],'unchanged_inputs':True,'historical_independent_bytes_match':binding_record['historical_independent_bytes_match'],'output':str(output)},indent=2))

if __name__=='__main__':main()
