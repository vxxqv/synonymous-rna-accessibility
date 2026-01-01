"""Freeze sequence-only eligibility, no phenotype workbook reads."""
from pathlib import Path
from collections import defaultdict,Counter
import json,hashlib
HERE=Path(__file__).resolve().parent;GATE=HERE.parent/'cprime_gate'
raw=json.loads((GATE/'selected_assay_metadata.json').read_text())['multi_metadata']
motifs=json.loads((GATE/'search_sequence_metadata.json').read_text())
byid={};bygroup=defaultdict(set);counts=Counter()
for row in raw:
    key=row['combined_var_id'];counts[key,row['cell_line']]+=1
    if 'synonymous_variant' not in (row['consequence'] or ''):continue
    info={k:row[k] for k in ['combined_var_id','var_id','hgvsc','hgvsp','exon_number','syn_var','consequence']}
    if key in byid and byid[key]!=info:raise ValueError('Conflicting metadata')
    byid[key]=info;bygroup[row['var_id']].add(key)
eligible={k:sorted(v) for k,v in sorted(bygroup.items()) if len(v)>=2}
ids=sorted({v for group in eligible.values() for v in group})
if len(eligible)!=30 or len(ids)!=66 or len(bygroup)!=63:raise ValueError('B0 metadata mismatch')
states={}
for key in ids:
    if any(counts[key,c]!=3 for c in ['K562','KBM7']):raise ValueError('Unexpected metadata multiplicity')
    words={m:sorted({r['motif'] for r in motifs if r['id']==key and r['modality']==m}) for m in ['DNA','RNA']}
    if any(len(v)!=1 for v in words.values()):raise ValueError('Unauthenticated local motif')
    states[key]={**byid[key],'motifs':words}
inputs={n:hashlib.sha256((GATE/n).read_bytes()).hexdigest() for n in ['selected_assay_metadata.json','search_sequence_metadata.json','measurement_contract.json','supp_Supplementary Data 8.xlsx']}
out={'status':'metadata_only','groups':eligible,'states':states,'cell_contexts':['K562','KBM7'],
 'expected_state_cell_rows':3,'groups_count':30,'states_count':66,'all_synonymous_groups':63,
 'single_marker_ineligible_groups':sorted(set(bygroup)-set(eligible)),'input_sha256':inputs,
 'exons':sorted({int(r['exon_number']) for r in states.values()}),'phenotype_values_inspected':False}
(HERE/'eligibility.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'groups':len(eligible),'states':len(ids),'exons':out['exons'],'outcomes_read':False}))
