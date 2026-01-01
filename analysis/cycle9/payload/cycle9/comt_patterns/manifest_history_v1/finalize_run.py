"""Seal existing results and evidence only; no outcome analysis."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

ROOT=Path(__file__).resolve().parents[3]
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def load(p):return json.loads(p.read_text(encoding='utf-8'))

freeze=load(P/'execution_freeze.json');run=load(P/'observed_v1/receipt.json');review=load(P/'reconstruction_v1/receipt.json')
for n,h in freeze['owned_files'].items():
    if sha(P/n)!=h:raise ValueError('Frozen source changed: '+n)
for role,h in run['input_hashes_after'].items():
    if sha(ROOT/freeze['input_paths'][role])!=h:raise ValueError('Input changed: '+role)
if review['maximum_exact_discrepancy']!='0':raise ValueError('Reference failure')
files=['observed_v1/result.json','observed_v1/receipt.json','reconstruction_v1/receipt.json',
       'tables_v1/all_134_pairs.csv','tables_v1/all_sites.csv','tables_v1/SNV_sites.csv',
       'results_report.md','REPLAY.md','execution_freeze.json','synthetic_receipt.json']
inputs=[{'path':freeze['input_paths'][k],'sha256':v,'sha256_after':run['input_hashes_after'][k]}
        for k,v in run['input_hashes_before'].items()]
manifest={'schema_version':2,'claim_id':'COMT_local_pair_order_exploratory_cycle9',
 'repository':{'commit':freeze['repository_revision_observed'],'dirty':False},
 'command':['work/env/Scripts/python.exe','-B','work/research_cycle9/comt_patterns/analyze.py']+
           [v for k,val in run['argv'].items() for v in ('--'+k,val)],
 'environment':{'software':[{'name':'CPython','version':platform.python_version()}],
                'hardware':platform.machine()+'; '+platform.platform()+'; standard-library sequential execution'},
 'mathematics':{'assertion_tested':'For each of134prespecified same-codon pairs, one fixed strict opposite orientation across all128 declared map/preparation combinations.',
                'coefficient_domain':'Rational representations of stored binary floats',
                'conventions':'Lexical allele orientation; shared D within each R/T case; no cross-assay label pairing; strict zero threshold; missing and mask mismatch unresolved.',
                'inputs':[x['path'] for x in inputs],
                'bounds':{'alleles':125,'sites':53,'pairs':134,'joint_cases':17152,'SNV_pairs':97},
                'non_claims':['calibrated biological uncertainty','physical protein amount','causal compensation','evolutionary trade-off','independent biological validation','predictor rescue']},
 'randomness':{'used':False,'generator':'','seed':None},
 'run':{'started_at':datetime.fromtimestamp((P/'observed_v1/run_started.json').stat().st_mtime,timezone.utc).isoformat(),
        'runtime_seconds':run['elapsed_seconds'],'exit_status':0,'status':'completed',
        'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':None,'max_threads':None}},
 'outputs':[{'path':rel(P/f),'sha256':sha(P/f),'kind':'result' if f.endswith('result.json') else 'evidence'} for f in files],
 'checks':['18syntheticchecks','200inventedreferencefields','20665directreconstructionchecks','17152jointcasesexact','allinputhashesunchanged'],
 'result':'1of134localpairs uniformly reverses; site149 c.330C>G versus c.330C>T.133unresolved.1of97SNVpairs. Conditional empirical radius0.003713370313814135nats.',
 'residual_risks':['Same-source exposed observations; shared normalization and biological preparations; representation uncertainty is not calibrated noise.',
                   'Direct verifier is separately formulated by the same agent, not a fresh reviewer.',
                   'Started timestamp uses the filesystem write time of the run_started marker written before input decoding.'],
 'input_artifacts':inputs,'declared_results':[]}
target=P/'computation_manifest.json'
if target.exists():raise FileExistsError(target)
target.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
all_files=sorted(p for p in P.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='final_manifest.json')
for p in all_files:
    if p.suffix in ('.md','.py') and '\u2014' in p.read_text(encoding='utf-8'):raise ValueError('Em dash in owned writing')
final={'status':'BOUNDED_COMPLETE_EXPLORATORY_RESULT_DISTINCT_RECONSTRUCTION_PASS',
       'files':{p.relative_to(P).as_posix():sha(p) for p in all_files},
       'observed_scientific_runs':1,'source_outcomes_newly_acquired':False,'frozen_sources_unchanged':True,
       'next_action':'Fresh separately authored review before integrating the narrow conditional result. No branch expansion.'}
(P/'final_manifest.json').write_text(json.dumps(final,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'final_manifest_sha256':sha(P/'final_manifest.json'),'report_sha256':sha(P/'results_report.md'),
                  'verifier_receipt_sha256':sha(P/'reconstruction_v1/receipt.json')}))
