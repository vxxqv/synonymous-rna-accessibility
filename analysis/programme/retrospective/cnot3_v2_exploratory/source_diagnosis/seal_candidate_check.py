import sys
if not __debug__: raise RuntimeError('Optimization refused')
sys.dont_write_bytecode=True
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
previous=json.loads((HERE/'manifest.json').read_text())
for name,value in previous['files_sha256'].items():
    if sha(HERE/name)!=value: raise RuntimeError('Prior sealed evidence changed: '+name)
receipt=json.loads((HERE/'candidate_universe_check.json').read_text())
if sha(HERE/'check_candidate_universe.py')!=receipt['code_sha256']:raise RuntimeError('Candidate source changed')
manifest={'status':'full_operational_candidate_membership_verified',
          'previous_manifest_sha256':sha(HERE/'manifest.json'),'previous_sealed_evidence_preserved':True,
          'added_files_sha256':{name:sha(HERE/name) for name in ['check_candidate_universe.py','candidate_universe_check.json','candidate_universe_addendum.md','seal_candidate_check.py']},
          'scientific_result':{'finite_IP_candidates':7359,'noncandidate_rows':15489,'candidate_sets_equal':True,'operational_paired_cohort':5480,'publication_count5481_still_unresolved':True},
          'failed_gate_modified':False,'scientific_protocol_modified':False,'modeling_or_grouping_run':False}
with (HERE/'candidate_universe_manifest.json').open('x',encoding='utf-8') as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
print(json.dumps({'manifest_sha256':sha(HERE/'candidate_universe_manifest.json'),'script_sha256':sha(HERE/'check_candidate_universe.py'),'receipt_sha256':sha(HERE/'candidate_universe_check.json')}))
