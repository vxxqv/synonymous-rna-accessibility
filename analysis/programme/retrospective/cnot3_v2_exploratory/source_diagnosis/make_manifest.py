"""Finalize bounded source-diagnosis provenance, no new endpoint evaluation."""
import sys
if not __debug__: raise RuntimeError('Optimization refused')
sys.dont_write_bytecode=True
import hashlib
import json
from pathlib import Path
import platform

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
latest=json.loads((HERE/'observed_v2/summary.json').read_text())
old=json.loads((HERE/'observed_v1/summary.json').read_text())
if sha(HERE/'audit_source.py')!=latest['code_sha256'] or sha(HERE/'audit_source.initial.py')!=old['code_sha256']: raise RuntimeError('Source/receipt mismatch')
identical={name:sha(HERE/'observed_v1'/name)==sha(HERE/'observed_v2'/name) for name in ('ratio_mismatches.csv','candidate_exclusions.json')}
if not all(identical.values()): raise RuntimeError('Unexpected arithmetic/selection change across interface cleanup')
manifest={'status':'bounded_source_diagnosis_complete','latest_evidence':'observed_v2/summary.json',
          'python':platform.python_version(),'synthetic_test_command':'python -I -B test_audit.py','synthetic_tests':5,'synthetic_failures':0,
          'independent_cohort_n':5480,'frozen_ratio_failure_n':49,'frozen_failure_retained':True,
          'diagnostic_10sig_all5480_compatible':True,'undocumented_generating_algorithm_authenticated':False,
          'published5481membership_recovered':False,'initial_to_v2_numerical_artifacts_identical':identical,
          'external_sha256':{'failed_principal_gate':sha(PARENT/'execution/principal_gate.json'),
                             'published_definition_documentation':sha(PARENT/'data_eligibility/Carneiro_gene_level_TableS2.xlsx.documentation.json'),
                             'published_article_text':sha(PARENT/'prior_review/carneiro_article.txt')},
          'input_sha256_before':latest['input_sha256_before'],'input_sha256_after':latest['input_sha256_after'],
          'files_sha256':{p.relative_to(HERE).as_posix():sha(p) for p in sorted(HERE.rglob('*')) if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts},
          'scope':{'authorized_candidate_rows':7359,'data_columns_decoded':['A','B','C','D','I','J','K'],'F_G_decoded':False,'other_endpoint_columns_decoded':False,'modeling_or_grouping':False,'scientific_protocol_changed':False,'repository_edits':False}}
with (HERE/'manifest.json').open('x',encoding='utf-8') as stream: json.dump(manifest,stream,indent=2); stream.write('\n')
print(json.dumps({'manifest_sha256':sha(HERE/'manifest.json'),'report_sha256':sha(HERE/'report.md'),'csv_sha256':sha(HERE/'observed_v2/ratio_mismatches.csv'),'script_sha256':sha(HERE/'audit_source.py')}))
