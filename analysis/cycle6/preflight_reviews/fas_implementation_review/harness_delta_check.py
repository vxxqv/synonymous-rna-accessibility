"""Check only the optimization guard and iterator regression amendment."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'fas_pilot'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected={'ordinal_geometry.py':'8872d6452afe4ccbf84e61450d4b1bfbdadf26469f20c0e24452f3abc205f9ac',
          'verify_synthetic.py':'a098549d48dcce62b4292cd433266e09265ba3f5c46bd08cf845c1a34e998656'}
for name,value in expected.items():
 if sha(SOURCE/name)!=value:raise ValueError('Source drift')
folder=HERE/'harness_delta';folder.mkdir(exist_ok=False)
runs=[]
for label,flag in [('normal',[]),('O',['-O']),('OO',['-OO'])]:
 d=folder/label;d.mkdir()
 for name in expected:shutil.copyfile(SOURCE/name,d/name)
 command=[sys.executable]+flag+[str((d/'verify_synthetic.py').resolve())]
 result=subprocess.run(command,capture_output=True,text=True,timeout=30)
 output=d/'synthetic_verification.json'
 if flag:
  if result.returncode==0 or output.exists() or 'requires Python assertions enabled' not in result.stderr:
   raise RuntimeError('Optimized invocation did not fail before evidence output')
 else:
  if result.returncode or not output.exists():raise RuntimeError('Normal suite failed')
  report=json.loads(output.read_text());prior=json.loads((HERE/'candidate_v1/synthetic_verification.json').read_text())
  retained=[c for c in report['checks'] if c['name']!='edge_iterator_is_materialized_before_reuse']
  if retained!=prior['checks'] or report['positive_fixture_certificate']!=prior['positive_fixture_certificate']:
   raise RuntimeError('Previous test results changed')
  if len(report['checks'])!=17:raise RuntimeError('Unexpected check count')
 runs.append({'label':label,'command':command,'exit_code':result.returncode,
              'output_exists':output.exists(),'stderr':result.stderr,
              'output_sha256':sha(output) if output.exists() else None})
for name,value in expected.items():
 if sha(SOURCE/name)!=value:raise ValueError('Source changed during review')
payload={'status':'PASS_NO_OPEN_BLOCKERS_IN_REVIEWED_DOMAIN','source_hashes':expected,'runs':runs,
         'normal_previous16results_unchanged':True,'normal_total_checks':17,
         'optimized_O_and_OO_refused_before_evidence_output':True,
         'production_numeric_source_unchanged_since_independent_checks':True,
         'scope':'Narrow synthetic harness delta only, no biological inputs',
         'script_sha256':sha(Path(__file__).resolve())}
(HERE/'delta_receipt.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':payload['status'],'receipt_sha256':sha(HERE/'delta_receipt.json')}))
