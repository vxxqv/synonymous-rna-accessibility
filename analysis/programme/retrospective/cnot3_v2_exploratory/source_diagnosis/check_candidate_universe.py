"""Independently authenticate every operational IP candidate. No HL-cell access."""
import sys
if not __debug__: raise RuntimeError('Optimization refused before input access')
sys.dont_write_bytecode=True
from collections import Counter
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
HELPER=HERE/'audit_source.py'
HELPER_SHA='fa3dfe63521b26fc23a3c8850ab2b5e2d5f99be47a23a6d030a3ce3858256be0'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding='utf-8'))
def payload_sha(value): return hashlib.sha256(json.dumps(value,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()

def main():
    target=HERE/'candidate_universe_check.json'
    if target.exists(): raise FileExistsError(target)
    if sha(HELPER)!=HELPER_SHA: raise RuntimeError('Previously audited independent SAX helper changed')
    spec=importlib.util.spec_from_file_location('independent_source_helper',HELPER)
    a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
    inventory_path=PARENT/'data_eligibility/reanalysis_gene_ids_only.json'
    exposure_path=PARENT/'execution/reproduction_v1/exposure.json'
    before={'workbook':sha(a.WORKBOOK),'inventory':sha(inventory_path),'exposure':sha(exposure_path),'helper':sha(HELPER)}
    if before['workbook']!=a.PINNED_WORKBOOK or before['exposure']!=a.PINNED_EXPOSURE or before['inventory']!='85fb60552f3c12477d28dbeeafbd171d52828ecb66e30591cefa60f77f482c52': raise RuntimeError('Authorized source hashes changed')
    inventory=load(inventory_path)
    source_rows={entry['row']:entry['id'] for entry in inventory['entries']}
    if len(source_rows)!=22848 or set(source_rows)!=set(range(2,22850)): raise RuntimeError('Not the complete authorized22848-row inventory')
    # Configure ONLY this newly imported, hash-pinned independent helper instance.
    # Original file bytes, original seven-column output and frozen failure stay unchanged.
    a.FIELDS={'A':'GENE_ID','I':'LOG2_FC_IP','J':'SE_IP','K':'STAT_IP'}
    cells=a.independent_cells(set(source_rows)|{1})
    for col,header in a.FIELDS.items():
        if cells[(1,col)]['decoded']!=header: raise RuntimeError('Header mismatch')
    if any(col not in a.FIELDS for _,col in cells): raise AssertionError('Protected column decoded')
    identifiers={row:cells.get((row,'A'),{}).get('decoded') for row in source_rows}
    if identifiers!=source_rows: raise AssertionError('Independent source GeneIDs differ from authenticated gene-only inventory')
    multiplicity=Counter(identifiers.values()); candidates=[]; exclusions=[]; reasons=Counter(); missing_patterns=Counter(); decimal_float_disagreements=[]
    for row in sorted(source_rows):
        gene=identifiers[row]; why=[]; nonfinite=[]
        if not isinstance(gene,str) or not gene or gene.strip()!=gene: why.append('invalid_gene_id')
        elif multiplicity[gene]!=1: why.append('duplicate_gene_id')
        for col in ('I','J','K'):
            value=a.number(cells.get((row,col)))
            decimal_finite=value is not None
            try: float_finite=value is not None and math.isfinite(float(value))
            except (OverflowError,ValueError): float_finite=False
            if decimal_finite!=float_finite: decimal_float_disagreements.append({'row':row,'column':col})
            if not float_finite: why.append('missing_or_nonfinite_'+a.FIELDS[col]); nonfinite.append(col)
        if why:
            exclusions.append(row); reasons.update(why); missing_patterns[','.join(nonfinite)]+=1
        else: candidates.append(row)
    # Compare only after independently deciding membership for all authorized rows.
    exposure=load(exposure_path)
    stages=[s for s in exposure['stages'] if s['columns']==['B','C','D']]
    if len(stages)!=1: raise RuntimeError('Missing frozen candidate selector')
    frozen=set(stages[0]['row_selector'])-{1}
    independently_selected=set(candidates)
    if independently_selected!=frozen: raise AssertionError({'newly_found_rows':sorted(independently_selected-frozen),'frozen_only_rows':sorted(frozen-independently_selected)})
    if len(candidates)!=7359 or len(exclusions)!=15489 or decimal_float_disagreements: raise AssertionError('Unexpected candidate completeness discrepancy')
    previous=load(HERE/'observed_v2/summary.json')
    after={'workbook':sha(a.WORKBOOK),'inventory':sha(inventory_path),'exposure':sha(exposure_path),'helper':sha(HELPER)}
    if after!=before: raise RuntimeError('Input mutation')
    receipt={'status':'independent_full_operational_candidate_universe_verified','source_data_rows':22848,
             'candidate_rows':len(candidates),'excluded_rows':len(exclusions),'candidate_row_set_exactly_matches_frozen':True,
             'candidate_gene_ID_count':len({identifiers[row] for row in candidates}),
             'candidate_row_set_sha256':payload_sha(sorted(candidates)),
             'candidate_gene_ID_set_sha256':payload_sha(sorted(identifiers[row] for row in candidates)),
             'missing_or_nonfinite_patterns':dict(missing_patterns),'exclusion_reason_counts':dict(reasons),
             'all_IDs_match_original_gene_only_extraction':True,'decimal_vs_float_finite_classification_differences':decimal_float_disagreements,
             'full_operational_cohort_n_authenticated_by_prior_HL_audit':previous['paired_cohort_independent'],
             'newly_decoded_columns':['A','I','J','K'],'B_C_D_outside_frozen_candidates_decoded':False,'F_G_or_other_columns_decoded':False,
             'new_biological_outcomes_exposed_beyond_original_reproduction_stage1':False,
             'method':'Hash-pinned independently authored SAX parser; membership decided over all22848 rows before frozen-selector comparison',
             'input_sha256_before':before,'input_sha256_after':after,'helper_source_sha256':HELPER_SHA,'code_sha256':sha(Path(__file__)),
             'prior_observed_v2_preserved':True,'failed_principal_gate_changed':False,'models_grouping_or_protocol_changes':False}
    with target.open('x',encoding='utf-8') as stream: json.dump(receipt,stream,indent=2); stream.write('\n')
    print(json.dumps({'status':receipt['status'],'candidate_rows':len(candidates),'excluded_rows':len(exclusions),'patterns':dict(missing_patterns),'receipt_sha256':sha(target)}))

if __name__=='__main__': main()
