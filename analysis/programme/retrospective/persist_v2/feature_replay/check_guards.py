"""Fast synthetic wrapper guard checks. Never imports RNA or opens reporter data."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, subprocess, sys
from types import SimpleNamespace
sys.dont_write_bytecode=True
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 here=Path(__file__).resolve().parent
 spec=importlib.util.spec_from_file_location('feature_replay_guard_test',here/'replay.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 out=m.absolute(a.output)
 if out.exists():raise RuntimeError('Refusing test overwrite')
 out.mkdir(parents=True,exist_ok=False);checks=[]
 def rejected(name,call):
  try:call()
  except RuntimeError as exc:checks.append({'name':name,'rejected':True,'message':str(exc)})
  else:raise RuntimeError('Expected rejection failed: '+name)
 for field in ['stdout','stderr']:
  p=SimpleNamespace(returncode=0,stdout='',stderr='');setattr(p,field,'synthetic warning')
  rejected('nonempty_'+field,lambda p=p:m.accept_process(p,'synthetic'))
 rejected('nonzero_return',lambda:m.accept_process(SimpleNamespace(returncode=1,stdout='',stderr=''),'synthetic'))
 m.accept_process(SimpleNamespace(returncode=0,stdout='',stderr=''),'synthetic');checks.append({'name':'quiet_zero_return','accepted':True})
 warning=out/'synthetic_warning';warning.mkdir()
 rejected('recorded_python_warning',lambda:m.save_warning_checked_row(warning,{'python_warnings':[{'category':'UserWarning','message':'synthetic only'}]}))
 if not (warning/'row.json').is_file():raise RuntimeError('Warning evidence not retained')
 existing=out/'existing';existing.mkdir();(existing/'sentinel').write_text('unchanged')
 rejected('worker_existing_output',lambda:m.worker(out/'nonexistent_repo',existing,'synthetic'))
 rejected('parent_existing_output',lambda:m.run(out/'nonexistent_repo',existing))
 for flag in ['-O','-OO']:
  destination=out/('not_created_'+flag[1:]);p=subprocess.run([sys.executable,flag,str(here/'replay.py'),'--repo-root',str(out/'absent_repo'),'--output',str(destination)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
  if p.returncode==0 or destination.exists() or 'Python -O and -OO are refused' not in p.stderr:raise RuntimeError('Optimization guard failed')
  checks.append({'name':flag,'rejected':True,'output_created':False})
 if 'RNA' in sys.modules or 'numpy' in sys.modules:raise RuntimeError('Guard tests imported scientific dependencies')
 record={'status':'passed','checks':checks,'folds':0,'reporter_inputs_opened':False,'source_sha256':hashlib.sha256((here/'replay.py').read_bytes()).hexdigest(),'synthetic_warning_row_retained':True,'initial_wrapper_folded_rows':0,'initial_warning_rows':'not applicable: initial wrapper never executed'}
 (out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'status':'passed','checks':len(checks),'folds':0},indent=2))
if __name__=='__main__':main()
