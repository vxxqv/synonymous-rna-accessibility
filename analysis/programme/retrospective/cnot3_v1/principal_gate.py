"""Fixed published-reproduction gate. No alternate selection or model fitting."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

SECONDARY={'PTMD_overlap_not_available_on_full_IP_cohort','published_PTMD_cutoff_tie_rule_unreported'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def positive(x):return type(x) in (int,float) and math.isfinite(x) and x>0

def assess(result,rows):
    ids=[r.get('GeneID') for r in rows]
    issues=result['discrepancies']
    primary=[x for x in issues if x not in SECONDARY]
    ratio=result['half_life_ratio_diagnostics']
    stat=result['IP_statistic_diagnostics']
    checks={
      'count_5481':len(rows)==result['n_cohort']==5481 and result['n_high']==1000 and result['n_remainder']==4481,
      'distinct_stable_ids':all(isinstance(x,str) and re.fullmatch(r'ENSG[0-9]{11}',x) for x in ids) and len(set(ids))==len(ids) and sorted(ids)==sorted(result['cohort_ids']),
      'positive_finite_paired_half_lives':all(positive(r.get(k)) for r in rows for k in ['HL_WT','HL_CNOT3KO']),
      'stored_ratio_arithmetic':len(ratio)==len(rows) and all(x.get('consistent') is True for x in ratio) and sorted(x['GENE_ID'] for x in ratio)==sorted(ids),
      'positive_finite_SE_IP':all(positive(r.get('SE_IP')) for r in rows),
      'STAT_IP_arithmetic':len(stat)==len(rows) and all(x.get('consistent') is True for x in stat) and sorted(x['GENE_ID'] for x in stat)==sorted(ids),
      'no_top1000_boundary_tie':result['selection_ties'].get('status')=='selected' and result['selection_ties'].get('boundary_tie_crosses') is False,
      'median_log2_change_top1000_greater_rest':positive(result['group_effects']['HL_FC_LOG2'].get('median_difference')),
      'no_negative_SE_among_IP_candidates':not any(x.startswith('negative_SE_IP:') for x in issues),
    }
    fields={'GeneID','HL_WT','HL_CNOT3KO','LOG2_FC_IP','SE_IP'}
    malformed=not all(set(r)==fields and type(r.get('LOG2_FC_IP')) in (int,float) and math.isfinite(r['LOG2_FC_IP']) for r in rows)
    if malformed:primary.append('malformed_normalized_cohort')
    return {'schema':'CNOT3_principal_reproduction_gate_v1','status':'passed' if all(checks.values()) and not primary else 'held_for_source_diagnosis',
      'cohort_gene_count':len(rows),'cohort_GeneIDs':sorted(ids),'principal_checks':checks,
      'principal_discrepancies':primary,'secondary_warnings':[x for x in issues if x in SECONDARY],
      'unresolved_source_rules':result['unresolved_source_rules'],
      'interpretation':'Operational Figure2 reproduction only. Original read-coverage and missingness rules remain unverified. A passed gate permits only the separately frozen conditional predictive comparison.'}

def execute(directory,protocol,out):
    directory=Path(directory);out=Path(out)
    if out.exists():raise FileExistsError('Do not overwrite a gate decision')
    receipt=load(directory/'receipt.json')
    for key,name in [('result_sha256','result.json'),('cohort_sha256','cohort_rows.json'),('exposure_sha256','exposure.json')]:
        if receipt[key]!=sha(directory/name):raise ValueError('Reproduction receipt mismatch: '+name)
    result=assess(load(directory/'result.json'),load(directory/'cohort_rows.json'))
    result.update(cohort_rows_sha256=sha(directory/'cohort_rows.json'),reproduction_result_sha256=sha(directory/'result.json'),reproduction_receipt_sha256=sha(directory/'receipt.json'),prediction_protocol_sha256=sha(protocol),implementation_sha256=sha(__file__))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reproduction',required=True);p.add_argument('--protocol',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();r=execute(a.reproduction,a.protocol,a.out)
    print(json.dumps({'status':r['status'],'principal_checks':r['principal_checks'],'principal_issue_count':len(r['principal_discrepancies']),'secondary_warning_count':len(r['secondary_warnings'])},indent=2))
    if r['status']!='passed':raise SystemExit(2)
