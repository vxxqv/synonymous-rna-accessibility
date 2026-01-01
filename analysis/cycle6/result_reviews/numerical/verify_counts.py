"""Independent exact count/normalization audit, using only the standard library.

The observed mode requires explicit authorization and hash-pinned saved inputs.
No production analysis module is imported. Synthetic mode never opens raw data.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import re
import subprocess
import sys
import time

WT='GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
ANCESTOR='GATCCAGATCTAACTTGCTGTGGTTGTGTCTCCTGCTTCTCCCGATTCTAGTAATTGTTTGGG'
DESIGN=re.compile(r'GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG')
COMMIT='9ed05ad8ab75c596e87a5ce7b63eb1bb38e71ba5'
CONTRACT_SHA='f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38'
RESULT_SHA='3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'


def need(condition,message):
 if not condition:raise RuntimeError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
def enc(v):return None if v is None else str(v)


def parse(raw,expected_blob):
 need(blob(raw)==expected_blob,'Git blob mismatch before decode')
 lines=raw.decode('ascii').splitlines();need(bool(lines),'Empty source')
 need(len(lines[0].split())==3 and lines[0].split()[0]==WT,'First source row is not native WT')
 result={}
 for line in lines:
  if not line.strip():continue
  row=line.split();need(len(row)==3,'Malformed row')
  seq,count,distance=row
  need(DESIGN.fullmatch(seq) is not None and seq not in result,'Unexpected/duplicate genotype')
  need(re.fullmatch('[0-9]+',count) is not None and re.fullmatch('[0-9]+',distance) is not None,'Bad integer')
  need(int(distance)==sum(a!=b for a,b in zip(seq,ANCESTOR)),'Ancestor annotation mismatch')
  result[seq]=int(count)
 need(bool(result) and sum(result.values())>0,'Empty/nonpositive library')
 return result


def midpoint_order_stat(values):
 if not values:return None
 s=sorted(values);n=len(s)
 return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2


def reference(libraries,sequences,wt=WT):
 inputs=[libraries[f'TR{i}'] for i in range(1,4)]
 outputs=[libraries[f'BR{i}'] for i in range(1,10)]
 universe=set(inputs[0]);all_libs=inputs+outputs
 # Directly condition each full library on the TR1 genotype universe.
 aligned=[sum(lib.get(s,0) for s in universe) for lib in all_libs]
 whole=[sum(lib.values()) for lib in all_libs]
 medians={}
 for s in set(sequences)|{wt}:
  available=[]
  for k,inp in enumerate(inputs):
   if s in inp and aligned[k]>0:available.append(Q(inp[s],aligned[k]))
  medians[s]=midpoint_order_stat(available)
 records=[]
 for index,s in enumerate(sequences):
  ratios=[];undefined=[]
  for b,out in enumerate(outputs,1):
   for t,inp in enumerate(inputs,1):
    if s not in out or wt not in out or s not in inp or wt not in inp:
     value=None;reason='missing_required_count'
    elif out[wt]==0 or inp[s]==0 or inp[wt]==0:
     value=None;reason='zero_required_denominator_or_input_reference'
    else:
     # Calculate two count enrichments and divide, rather than using cross-products.
     value=Q(out[s],inp[s])/Q(out[wt],inp[wt]);reason=None
    ratios.append(value)
    if reason:undefined.append({'BR':b,'TR':t,'reason':reason})
  source=[]
  for b,out in enumerate(outputs):
   total=aligned[b+3]
   if s not in universe or wt not in universe or s not in out or wt not in out or total<=0 or out[wt]<=0 or medians[s] is None or medians[wt] is None or medians[s]<=0 or medians[wt]<=0:
    source.append(None)
   else:
    # Explicit normalized output frequencies and median input enrichments.
    focal_enrichment=Q(out[s],total)/medians[s]
    control_enrichment=Q(out[wt],total)/medians[wt]
    source.append(focal_enrichment/control_enrichment)
  complete=all(v is not None for v in ratios)
  lower=min(ratios) if complete else Q(0)
  upper=max(ratios) if complete else None
  point=sum(source,Q(0))/9 if all(v is not None for v in source) else None
  if upper is not None and point is not None:need(lower<=point<=upper,'Corrected point outside complete envelope')
  records.append({'vertex':index,'sequence':s,'output_counts':[out.get(s) for out in outputs],
                  'input_counts':[inp.get(s) for inp in inputs],
                  'scenario_ratios':[enc(v) for v in ratios],'lower':enc(lower),'upper':enc(upper),
                  'author_relative_by_output':[enc(v) for v in source],'point':enc(point),
                  'undefined_scenarios':undefined})
 return {'records':records,'normalization':{'TR1_aligned_totals':aligned,'whole_file_totals':whole,
          'total_order':'TR1..TR3, BR1..BR9','aligned_genotype_count':len(universe),
          'undefined_aligned_frequency_indices':[i for i,t in enumerate(aligned) if t<=0]}}


def synthetic():
 libraries={'TR1':{'W':2,'A':1,'X':1},'TR2':{'W':4,'A':4,'X':0,'Z':999},
            'TR3':{'W':8,'A':2,'X':6,'Z':333}}
 libraries.update({f'BR{i}':{'W':2,'A':3,'X':5,'Z':1000} for i in range(1,10)})
 r=reference(libraries,['W','A'],'W');a=r['records'][1]
 need(a['scenario_ratios']==['3','3/2','6']*9 and a['point']=='3','Known-answer ratios/point')
 need(a['lower']=='3/2' and a['upper']=='6','Known envelope')
 need(r['normalization']['TR1_aligned_totals']==[4,8,16]+[10]*9,'TR1 universe alignment')
 need(r['normalization']['whole_file_totals']!=r['normalization']['TR1_aligned_totals'],'Outside universe did not affect whole totals')
 checks=['known27ratios','median_source_point','TR1_alignment','whole_vs_aligned_totals']
 for case in ['missing_TR2','zero_TR2','zero_output','missing_output','zero_WT_output','missing_TR1']:
  d={k:v.copy() for k,v in libraries.items()}
  if case=='missing_TR2':del d['TR2']['A']
  if case=='zero_TR2':d['TR2']['A']=0
  if case=='zero_output':d['BR1']['A']=0
  if case=='missing_output':del d['BR1']['A']
  if case=='zero_WT_output':d['BR1']['W']=0
  if case=='missing_TR1':del d['TR1']['A']
  a=reference(d,['W','A'],'W')['records'][1]
  if case=='zero_output':need(a['scenario_ratios'][:3]==['0']*3 and a['upper'] is not None and a['point']=='8/3','Zero output is observed zero')
  else:need(a['upper'] is None and a['lower']=='0','Missing/denominator-zero not conservative unknown')
  if case in ['missing_output','zero_WT_output','missing_TR1']:need(a['point'] is None,'Undefined point not preserved')
  if case=='missing_TR2':need(a['point']=='4','Two-input median not exact midpoint')
  checks.append(case)
 return {'status':'PASS','checks':checks,'scope':'Pure synthetic dictionaries; no raw source or result opened'}


def observed(args):
 need(args.authorized_observed,'Observed mode requires root authorization flag')
 need(args.repo and args.input_directory and args.result,'Explicit observed paths required')
 repo=args.repo.resolve();src=repo/'analysis/cycle6/fas_pilot'
 commit=subprocess.check_output(['git','-C',str(repo),'rev-parse',COMMIT+'^{commit}'],text=True).strip()
 need(commit==COMMIT,'Amendment commit unavailable')
 # Verify both immutable freeze manifests and each scientific source, before count decoding.
 core_freeze_raw=(src/'freeze.json').read_bytes();format_freeze_raw=(src/'format_freeze.json').read_bytes()
 for name in ['freeze.json','format_freeze.json','run_headerless.py','format_amendment_v1.md']:
  committed=subprocess.check_output(['git','-C',str(repo),'show',f'{COMMIT}:analysis/cycle6/fas_pilot/{name}'])
  need(committed==(src/name).read_bytes(),'Committed amendment byte drift: '+name)
 freeze=json.loads(core_freeze_raw);amendment=json.loads(format_freeze_raw)
 need(amendment['scientific_freeze_sha256']==sha(core_freeze_raw),'Science freeze mismatch')
 pins=dict(freeze['frozen_files']);pins.update(amendment['files'])
 for name,value in pins.items():need(sha((src/name).read_bytes())==value,'Pinned implementation changed: '+name)
 contract_raw=(src/'source_contract.json').read_bytes();need(sha(contract_raw)==CONTRACT_SHA,'Contract drift')
 cfg=json.loads(contract_raw);sequences=cfg['sequences']
 need(len(sequences)==72 and len(set(sequences))==72 and WT in sequences,'Wrong full domain')
 result_raw=args.result.read_bytes();need(sha(result_raw)==RESULT_SHA,'Saved result drift')
 result=json.loads(result_raw)
 libraries={};receipts=[]
 for item in cfg['files']:
  path=args.input_directory/item['local_name'];raw=path.read_bytes()
  need(len(raw)==item['bytes'],'Original source length drift')
  libraries[item['label']]=parse(raw,item['git_blob'])
  receipts.append({'label':item['label'],'git_blob':item['git_blob'],'sha256':sha(raw),
                   'bytes':len(raw),'observed_rows':len(libraries[item['label']])})
 need(set(libraries)=={f'BR{i}' for i in range(1,10)}|{f'TR{i}' for i in range(1,4)},'Source labels')
 expected=reference(libraries,sequences)
 need(result['source_receipts']==receipts,'Source receipts disagree')
 need(result['freeze_sha256']==sha(core_freeze_raw),'Result scientific freeze hash mismatch')
 need(result['normalization']==expected['normalization'],'Normalization metadata discrepancy')
 need(len(result['vertices_data'])==72,'Saved row count mismatch')
 checks=0;undefined=[]
 fields=['vertex','sequence','output_counts','input_counts','scenario_ratios','lower','upper','author_relative_by_output','point']
 for calculated,actual in zip(expected['records'],result['vertices_data']):
  for field in fields:
   need(actual[field]==calculated[field],f'Exact disagreement at vertex{calculated["vertex"]}, field{field}')
   checks+=len(calculated[field]) if isinstance(calculated[field],list) else 1
  if calculated['undefined_scenarios']:undefined.append({'vertex':calculated['vertex'],'reasons':calculated['undefined_scenarios']})
 finite=sum(r['upper'] is not None for r in expected['records'])
 need(result['finite_envelopes']==finite and result['unknown_envelopes']==72-finite,'Eligibility counts differ')
 scenario_labels=[f'BR{b}_TR{t}' for b in range(1,10) for t in range(1,4)]+['mean_author_per_output_relative_ES']
 need([r['scenario'] for r in result['scenario_results']]==scenario_labels,'Scenario indexing changed')
 for index,row in enumerate(result['scenario_results']):
  unknown=sum((r['scenario_ratios'][index] if index<27 else r['point']) is None for r in expected['records'])
  need(row['unknown_vertices']==unknown,'Scenario dropout differs')
 # No result changes or source changes during verification.
 need(args.result.read_bytes()==result_raw,'Result changed during verification')
 for item,receipt in zip(cfg['files'],receipts):need(sha((args.input_directory/item['local_name']).read_bytes())==receipt['sha256'],'Source changed during audit')
 return {'status':'EXACT_REFERENCE_AGREEMENT','amendment_commit':COMMIT,'result_sha256':RESULT_SHA,
         'contract_sha256':CONTRACT_SHA,'scientific_freeze_sha256':sha(core_freeze_raw),
         'format_freeze_sha256':sha(format_freeze_raw),'source_receipts':receipts,
         'selected_genotypes':72,'scenario_ratios_checked':72*27,'per_output_corrected_points_checked':72*9,
         'corrected_mean_points_checked':72,'selected_count_entries_checked':72*12,
         'envelope_endpoints_checked':72*2,'selected_field_scalar_comparisons':checks,
         'normalization_metadata_exact':True,'scenario_unknown_counts_checked':28,
         'finite_envelopes':finite,'unknown_envelopes':72-finite,'undefined_cases':undefined,
         'maximum_exact_numerical_discrepancy':'0',
         'scope':'Independent exact rational count/normalization/eligibility audit. Ordinal decisions and graph certificates are not recalculated here.',
         'point_interpretation':'Corrected headerless source-formula reconstruction, not authenticated historical source result.',
         'exposure':'Verifier decoded all twelve pinned source count files for TR1-aligned totals; all72 selected records verified. Only compact audit metadata reported, no full data dump.',
         'independence':'27 scenarios share nine outputs and three technical inputs; they are not27 independent biological replicates.',
         'no_new_models_or_selection':True}


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['synthetic','observed'],required=True)
 parser.add_argument('--out',type=Path,required=True);parser.add_argument('--repo',type=Path)
 parser.add_argument('--input-directory',type=Path);parser.add_argument('--result',type=Path)
 parser.add_argument('--authorized-observed',action='store_true');args=parser.parse_args()
 need(not args.out.exists(),'Output already exists');start=time.perf_counter()
 result=synthetic() if args.mode=='synthetic' else observed(args)
 result['verifier_sha256']=sha(Path(__file__).read_bytes());result['runtime_seconds']=time.perf_counter()-start
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:result[k] for k in ['status','runtime_seconds']}
                  |{'receipt_sha256':sha(args.out.read_bytes())}))


if __name__=='__main__':main()
