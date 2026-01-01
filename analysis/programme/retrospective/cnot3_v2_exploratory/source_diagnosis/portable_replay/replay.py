"""Portable replay of two immutable source-diagnosis scripts. Standard library only."""
import sys
if not __debug__: raise RuntimeError('Python -O/-OO is refused before input access')
sys.dont_write_bytecode=True
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time
import traceback

SOURCE='analysis/programme/retrospective/cnot3_v2_exploratory/source_diagnosis/'
REPRO='results/programme/cnot3_v1/reproduction/'
DATA='data/programme/cnot3_v1/'
PINS={
    SOURCE+'audit_source.py':'fa3dfe63521b26fc23a3c8850ab2b5e2d5f99be47a23a6d030a3ce3858256be0',
    SOURCE+'check_candidate_universe.py':'21f02e6e2fce5e4a6ca39e2ff1dafb58a3a20e2f70de1d449f4bdab2af0f576a',
    DATA+'Carneiro_gene_level_TableS2.xlsx':'a338e7211fcb0ad6bed0b1bfdb57433a5a43266f3a8d93e473b1bf50c583d712',
    DATA+'reanalysis_gene_ids_only.json':'85fb60552f3c12477d28dbeeafbd171d52828ecb66e30591cefa60f77f482c52',
    REPRO+'result.json':'67fa85f42354c8f6413b167d38a55ca3630f7ea12d3b2dad6658b302bb3f2376',
    REPRO+'exposure.json':'b8bae6cc9e71f7740f54a304b6182c108fa2e69d6fc577660e2522ec86490159',
    REPRO+'cohort_rows.json':'cfd4cededb5be6719795fa88262793ca580eb267649127cdf9cf300f75d9df1e',
    SOURCE+'observed_v2/summary.json':'b4279c33b1f425229a3267ac92b21d9589cc65904b93be69db0846125dbcc635',
    SOURCE+'observed_v2/ratio_mismatches.csv':'ca8a1e6c0e2f056316912e29a0b1b391bf03cb9d5a5f2c56fae969944ebe7ae4',
    SOURCE+'observed_v2/candidate_exclusions.json':'498995d9dc0f77708547181d0d592780deab1680b10c518cfb18ce8d1bef9278',
    SOURCE+'candidate_universe_check.json':'8e2788e581f5729af27b1ca0bcc473b3233792dee82a9608f75a94bfc49667b5',
}
STAGE_MAP={
    SOURCE+'audit_source.py':'source_discrepancy_review/audit_source.py',
    SOURCE+'check_candidate_universe.py':'source_discrepancy_review/check_candidate_universe.py',
    DATA+'Carneiro_gene_level_TableS2.xlsx':'data_eligibility/Carneiro_gene_level_TableS2.xlsx',
    DATA+'reanalysis_gene_ids_only.json':'data_eligibility/reanalysis_gene_ids_only.json',
    REPRO+'result.json':'execution/reproduction_v1/result.json',
    REPRO+'exposure.json':'execution/reproduction_v1/exposure.json',
    REPRO+'cohort_rows.json':'execution/reproduction_v1/cohort_rows.json',
}

def absolute(path):
    resolved=Path(path).resolve()
    if os.name=='nt' and not str(resolved).startswith('\\\\?\\'):
        return Path('\\\\?\\'+str(resolved))
    return resolved

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk:=stream.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))

def save(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,allow_nan=False); stream.write('\n')

def compare_csv(actual,expected):
    def rows(path):
        with Path(path).open(encoding='utf-8',newline='') as stream: return list(csv.reader(stream))
    a,b=rows(actual),rows(expected)
    if a!=b: raise RuntimeError('Complete mismatch CSV differs from archived evidence')
    return {'all_headers_fields_rows_equal':True,'data_rows':len(a)-1}

def run(repo,output):
    # Resolve aliases before containment checks. Refuse before reading archived inputs.
    repo_normal=Path(repo).resolve(); output_normal=Path(output).resolve()
    if output_normal.exists(): raise ValueError('Output already exists; no overwrite')
    if output_normal==repo_normal or output_normal.is_relative_to(repo_normal): raise ValueError('Output must be outside the input repository')
    if not repo_normal.is_dir(): raise ValueError('Input repository/archive directory is absent')
    if platform.python_implementation()!='CPython' or platform.python_version()!='3.12.14': raise RuntimeError('Preinstalled CPython3.12.14 is required')
    repo,output=absolute(repo_normal),absolute(output_normal)
    before={name:sha(repo/name) for name in PINS}
    if before!=PINS: raise RuntimeError('Immutable archived input or original script hash mismatch')
    output.mkdir(parents=True)
    start=time.monotonic(); stage=output/'stage/decay_cycle4'
    save(output/'STARTED.json',{'status':'started','archived_inputs_sha256':before,'wrapper_sha256':sha(absolute(__file__))})
    try:
        staged={}
        for source,destination in STAGE_MAP.items():
            target=stage/destination; target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(repo/source,target)
            staged[destination]=sha(target)
            if staged[destination]!=PINS[source]: raise RuntimeError('Stage copy changed frozen bytes')
        review=stage/'source_discrepancy_review'
        commands=[('ratio_audit',[sys.executable,'-I','-B',str(review/'audit_source.py'),'--out',str(review/'observed_v2')]),
                  ('candidate_universe',[sys.executable,'-I','-B',str(review/'check_candidate_universe.py')])]
        processes=[]
        for name,argv in commands:
            proc=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,encoding='utf-8',errors='replace')
            (output/(name+'.stdout.txt')).write_text(proc.stdout,encoding='utf-8')
            (output/(name+'.stderr.txt')).write_text(proc.stderr,encoding='utf-8')
            processes.append({'name':name,'returncode':proc.returncode,'stdout_sha256':sha(output/(name+'.stdout.txt')),'stderr_sha256':sha(output/(name+'.stderr.txt'))})
            if proc.returncode!=0 or proc.stderr: raise RuntimeError('Original audited script failed or emitted stderr: '+name)
        comparisons={}
        for actual,expected,label in [
            (review/'observed_v2/summary.json',repo/(SOURCE+'observed_v2/summary.json'),'complete_ratio_summary'),
            (review/'candidate_universe_check.json',repo/(SOURCE+'candidate_universe_check.json'),'complete_candidate_summary'),
            (review/'observed_v2/candidate_exclusions.json',repo/(SOURCE+'observed_v2/candidate_exclusions.json'),'complete_candidate_exclusion_ledger'),
        ]:
            if load(actual)!=load(expected): raise RuntimeError('Complete JSON content mismatch: '+label)
            comparisons[label]={'all_fields_equal':True,'provenance_fields_ignored':[],'byte_identical':sha(actual)==sha(expected)}
        comparisons['complete_ratio_mismatch_CSV']=compare_csv(review/'observed_v2/ratio_mismatches.csv',repo/(SOURCE+'observed_v2/ratio_mismatches.csv'))
        after={name:sha(repo/name) for name in PINS}
        staged_after={destination:sha(stage/destination) for destination in STAGE_MAP.values()}
        if before!=after or staged_after!=staged: raise RuntimeError('Repository or staged immutable input mutation')
        summary=load(review/'observed_v2/summary.json'); universe=load(review/'candidate_universe_check.json')
        receipt={'status':'portable_source_audit_replay_passed','comparisons':comparisons,
                 'paired_cohort':summary['paired_cohort_independent'],'frozen_ratio_failures':summary['cohort_ratio_failures'],
                 'finite_IP_candidates':universe['candidate_rows'],'excluded_IP_rows':universe['excluded_rows'],
                 'original_failed_gate_not_revised':True,'original_script_hashes_preserved':True,
                 'input_sha256_before':before,'input_sha256_after':after,'staged_frozen_files_sha256':staged_after,
                 'original_processes':processes,'runtime_seconds':time.monotonic()-start,
                 'environment':{'implementation':platform.python_implementation(),'python':platform.python_version(),'platform':platform.platform()},
                 'wrapper_sha256':sha(absolute(__file__)),'runtime_Git':False,'runtime_network':False,'runtime_installation':False,
                 'copied_source_files':len(STAGE_MAP),'archived_evidence_files_compared':4,
                 'no_new_columns_models_grouping_or_thresholds':True}
        save(output/'receipt.json',receipt)
        return receipt
    except BaseException as exc:
        save(output/'FAILURE.json',{'status':'portable_replay_failed','exception':repr(exc),'traceback':traceback.format_exc(),'scientific_outputs_must_not_be_accepted':True})
        raise

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--repo',required=True); parser.add_argument('--out',required=True)
    args=parser.parse_args(); receipt=run(args.repo,args.out)
    print(json.dumps({'status':receipt['status'],'paired_cohort':receipt['paired_cohort'],'ratio_failures':receipt['frozen_ratio_failures']}))

if __name__=='__main__': main()
