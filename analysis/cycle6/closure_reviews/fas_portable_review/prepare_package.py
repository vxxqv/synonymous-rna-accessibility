"""Create owned offline adapters and a minimal relocated archive fixture."""
from pathlib import Path
import ast
import difflib
import hashlib
import json
import os

BASE = Path(__file__).resolve().parent
WORKSPACE = BASE.parents[2]
REPO = WORKSPACE / 'outputs/synonymous-rna-accessibility'

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def native(path):
    value = str(path.resolve())
    if os.name == 'nt' and not value.startswith('\\\\?\\'):
        value = '\\\\?\\UNC\\' + value[2:] if value.startswith('\\\\') else '\\\\?\\' + value
    return Path(value)

def fresh_write(path, raw):
    path = native(path)
    if path.exists():
        raise RuntimeError('Refusing existing file: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)

def main():
    numeric_path = WORKSPACE / 'work/research_cycle6/fas_result_numerical_audit/verify_counts.py'
    original = numeric_path.read_bytes()
    if digest(original) != 'f5a932166d3cd8ddefc47020e190610b0b845aea9a72a1d6a3531f9693ffd840':
        raise RuntimeError('Original numerical source changed')
    text = original.decode('utf-8')
    old = ''' commit=subprocess.check_output(['git','-C',str(repo),'rev-parse',COMMIT+'^{commit}'],text=True).strip()
 need(commit==COMMIT,'Amendment commit unavailable')
 # Verify both immutable freeze manifests and each scientific source, before count decoding.
 core_freeze_raw=(src/'freeze.json').read_bytes();format_freeze_raw=(src/'format_freeze.json').read_bytes()
 for name in ['freeze.json','format_freeze.json','run_headerless.py','format_amendment_v1.md']:
  committed=subprocess.check_output(['git','-C',str(repo),'show',f'{COMMIT}:analysis/cycle6/fas_pilot/{name}'])
  need(committed==(src/name).read_bytes(),'Committed amendment byte drift: '+name)
'''
    new = ''' # Offline replay verifies the archived amendment bytes, not local Git history.
 core_freeze_raw=(src/'freeze.json').read_bytes();format_freeze_raw=(src/'format_freeze.json').read_bytes()
 offline_pins = {
  'freeze.json':'7bddb75dcdb4b401a120dcb518bd34fcfb5f6a801a492ee1a59f1fced90be3eb',
  'format_freeze.json':'d3b69ed44015a712d4dbe436cdb7399c009e521f3444a6a5736c547a04094041',
  'run_headerless.py':'11addfe6ee9c02f7fc17b9b6a414167cc0d98baccbab71bf7fa899fd7b713bb6',
  'format_amendment_v1.md':'73a89257046695a1882c01b98971018d78aba95a144afb5e6752a9e10d2b9a62'}
 for name,expected_hash in offline_pins.items():
  need(sha((src/name).read_bytes())==expected_hash,'Archived amendment byte drift: '+name)
'''
    if text.count(old) != 1:
        raise RuntimeError('Expected Git-only block not unique')
    adapted = text.replace('import subprocess\n', '').replace(old, new)
    # Preserve every scientific expression and recorded result. Add only provenance metadata.
    old_return = "return {'status':'EXACT_REFERENCE_AGREEMENT','amendment_commit':COMMIT,'result_sha256':RESULT_SHA,"
    new_return = "return {'status':'EXACT_REFERENCE_AGREEMENT','git_history_checked_in_this_replay':False,'amendment_commit':COMMIT,'result_sha256':RESULT_SHA,"
    if adapted.count(old_return) != 1:
        raise RuntimeError('Expected result provenance insertion not unique')
    adapted = adapted.replace(old_return, new_return)
    old_nodes = {node.name:ast.dump(node,include_attributes=False) for node in ast.parse(text).body if isinstance(node,ast.FunctionDef)}
    new_nodes = {node.name:ast.dump(node,include_attributes=False) for node in ast.parse(adapted).body if isinstance(node,ast.FunctionDef)}
    unchanged = [name for name in old_nodes if old_nodes[name] == new_nodes[name]]
    changed = [name for name in old_nodes if old_nodes[name] != new_nodes[name]]
    if changed != ['observed']:
        raise RuntimeError('Unexpected function change: ' + repr(changed))
    # Compare all statements after the gate up to the return object exactly.
    def numerical_body(source):
        body = next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='observed').body
        start = next(i for i,n in enumerate(body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='freeze' for t in n.targets))
        return [ast.dump(n,include_attributes=False) for n in body[start:-1]]
    if numerical_body(text) != numerical_body(adapted):
        raise RuntimeError('Observed numerical body changed')
    items = {
        'numerical/verify_counts.py':original,
        'numerical/verify_counts_offline.py':adapted.encode('utf-8'),
        'numerical/observed_receipt.json':(numeric_path.parent/'observed_receipt.json').read_bytes(),
        'numerical/synthetic_receipt.json':(numeric_path.parent/'synthetic_receipt.json').read_bytes(),
        'mathematical/verify_result.py':(WORKSPACE/'work/research_cycle6/fas_result_math_audit/verify_result.py').read_bytes(),
        'mathematical/verification.json':(WORKSPACE/'work/research_cycle6/fas_result_math_audit/observed_v1/verification.json').read_bytes(),
    }
    if digest(items['mathematical/verify_result.py']) != '047846795d9383eeef60ad8f74887783f706406fde0f4a75f916bd33af667db2':
        raise RuntimeError('Original mathematical source changed')
    for name,raw in items.items():
        fresh_write(BASE/'package'/name,raw)
    fresh_write(BASE/'adapter.diff',''.join(difflib.unified_diff(text.splitlines(True),adapted.splitlines(True),fromfile='numerical/verify_counts.py',tofile='numerical/verify_counts_offline.py')).encode('utf-8'))
    contract = json.loads((REPO/'analysis/cycle6/fas_pilot/source_contract.json').read_bytes())
    freeze = json.loads((REPO/'analysis/cycle6/fas_pilot/freeze.json').read_bytes())
    paths = ['analysis/cycle6/fas_pilot/'+name for name in freeze['frozen_files']]
    paths += ['analysis/cycle6/fas_pilot/'+name for name in ['freeze.json','format_freeze.json','run_headerless.py','format_amendment_v1.md']]
    paths += ['data/cycle6/fas_counts/'+item['local_name'] for item in contract['files']]
    paths += ['results/cycle6/full_result_v1.json.gz']
    manifest_files = []
    for name in paths:
        raw = native(REPO/name).read_bytes()
        manifest_files.append({'path':name,'sha256':digest(raw),'bytes':len(raw)})
        fresh_write(BASE/'r'/name,raw)
    for name,raw in items.items():
        dest='analysis/cycle6/result_reviews/'+name
        manifest_files.append({'path':dest,'sha256':digest(raw),'bytes':len(raw)})
        fresh_write(BASE/'r'/dest,raw)
    manifest = {'schema':'fas_independent_offline_replay_v1','files':manifest_files,
                'decompressed_result_sha256':'3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466',
                'decompressed_result_bytes':29145849,
                'original_amendment_commit':'9ed05ad8ab75c596e87a5ce7b63eb1bb38e71ba5',
                'git_history_required':False,'network_required':False}
    raw=(json.dumps(manifest,indent=2)+'\n').encode('utf-8')
    fresh_write(BASE/'package/replay_manifest.json',raw)
    fresh_write(BASE/'r/analysis/cycle6/result_reviews/replay_manifest.json',raw)
    receipt={'unchanged_functions':unchanged,'changed_functions':changed,'observed_numerical_body_ast_identical':True,
             'changes':'Git provenance gate replaced by exact archived-byte hashes; one explicit Git-history-not-checked metadata field.',
             'original_sources_preserved':True,'manifest_sha256':digest(raw),'input_files':len(manifest_files),
             'files':{name:digest(raw) for name,raw in items.items()}}
    fresh_write(BASE/'preparation_receipt.json',(json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(receipt))

if __name__=='__main__':
    main()
