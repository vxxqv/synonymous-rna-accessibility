from pathlib import Path
import hashlib,json,shutil
H=Path(__file__).resolve().parent;P=H/'comt_pilot';I=P/'implementation_v3';R=H.parents[1]/'outputs/synonymous-rna-accessibility';A=R/'analysis/cycle8'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
review=H/'comt_numeric_reference/v3_delta_review/review_receipt.json'
assert sha(review)=='6a9bd81313dddb0d134aaf7774f2ea4c03b7351658a4b1c6cfb49c9210b8377e'
rev=json.loads(review.read_text());assert rev['status']=='PASS_BOUNDED_V3_DELTA_REVIEW' and not rev['remaining_blocking_findings']
candidate=I/'candidate_manifest_v3.json';assert sha(candidate)==rev['accepted_candidate_manifest_sha256']
for n,d in json.loads(candidate.read_text())['files'].items():assert sha(I/n)==d,n
roles={'protocol':P/'protocol_v2.md','amendment':P/'prospective_reconstruction_amendment.md','integration_addendum':I/'integration_addendum.md','contract':I/'selective_contract.json','normalization_amendment':P/'normalization_contract_amendment_v3.md','predecessor_code_manifest':P/'code_freeze_v1.json'}
files={p.name:sha(p) for p in sorted(I.glob('*.py'))}
for n in ['apply_source_rf.R','../proxy_models.py','../source_helpers.R','../modeling_utils.r','../summarization_utils.r']:files[n]=sha(I/n)
code={'schema':'COMT_A1_code_freeze_v1','version':'3_source_contract_successor','files':files,'named_input_hashes':{k:sha(v) for k,v in roles.items()}}
codepath=P/'code_freeze_v3.json';freezepath=P/'normalization_successor_freeze_v3.json'
assert not codepath.exists() and not freezepath.exists()
write(codepath,code)
freeze={'schema':'COMT_source_contract_successor_freeze_v3','status':'FROZEN_AWAITING_PUBLIC_COMMIT','candidate_sha256':sha(candidate),'code_freeze_sha256':sha(codepath),'review_receipt_sha256':sha(review),'input_pins':code['named_input_hashes'],
 'original_freeze_sha256':sha(P/'final_preoutcome_freeze.json'),'original_failure_log_sha256':sha(H/'comt_results/actual_v1/normalize_execution.log'),
 'extract_manifest_sha256':sha(H/'comt_results/actual_v1/extract/sealed_manifest.json'),'rf_manifest_sha256':sha(H/'comt_results/actual_v1/rf/sealed_manifest.json'),
 'exposure':'Extraction/RF plus identities, DP and fixed masks inspected; no normalized contrasts or performance evaluated. This is a post-diagnostic successor, not original prospective execution.',
 'scientific_delta':'Restore source sorted-VAR_ID first-match depth selection. Preserve all target/fold/model/threshold rules and the negative-WT guard.',
 'historical_reproduction_claim':False,'source_defined_WT_is_physical_count_claim':False,
 'preserved_primary_rule':'At least100variants/40sites and10percent improvement over everyfourbaselines underbothmaps; midsensitivitycannotrescue.'}
write(freezepath,freeze)
allowed={'.py','.r','.md','.json','.csv','.tsv','.txt','.log','.diff','.rds','.gz'}
def cp(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
def tree(p,q):
 for x in sorted(p.rglob('*')):
  if x.is_file() and '__pycache__' not in x.parts and x.suffix.lower() in allowed:cp(x,q/x.relative_to(p))
tree(I,A/'comt_pilot/implementation_v3')
for p in [roles['normalization_amendment'],codepath,freezepath]:cp(p,A/'comt_pilot'/p.name)
tree(H/'comt_depth_failure_diagnosis',A/'comt_depth_failure_diagnosis')
tree(H/'comt_numeric_reference/failure_audit',A/'comt_numeric_reference/failure_audit')
tree(H/'comt_numeric_reference/v3_delta_review',A/'comt_numeric_reference/v3_delta_review')
tree(H/'comt_rf_independent_audit',A/'comt_rf_independent_audit')
for p in (H/'comt_results/actual_v1').iterdir():
 if p.is_file():cp(p,A/'comt_results/actual_v1'/p.name)
for stage in ['extract','rf']:
 cp(H/f'comt_results/actual_v1/{stage}/sealed_manifest.json',A/f'comt_results/actual_v1/{stage}/sealed_manifest.json')
for n in ['model_1.rds','model_2.rds','status.txt','execution_log.txt','rejections.tsv']:
 cp(H/'comt_results/actual_v1/rf'/n,A/'comt_results/actual_v1/rf'/n)
for p in [H/'freeze_comt_v3.py',H/'run_comt_successor.py']:
 if p.exists():cp(p,A/p.name)
write(A/'comt_results/actual_v1/PAYLOAD_RECONSTRUCTION.json',{'note':'Large extracted rows and full RF row tables remain reproducible from the externally pinned source archive and saved forests. Their complete stage hashes are retained. Original failure and independent trace review are retained. No publishing assets or source archive bundled.','source_archive_sha256':'e3fc127f268c10efbb1d1c736423c51bc5aa900e9e363d66cc8683490516cbe2'})
print(json.dumps({'freeze_sha256':sha(freezepath),'code_sha256':sha(codepath),'review_sha256':sha(review)}))
