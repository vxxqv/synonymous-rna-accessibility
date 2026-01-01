"""Synthetic-only count-reader and normalization checks before data access."""
if not __debug__:
    raise RuntimeError('Verification requires Python assertions enabled')
from pathlib import Path
from fractions import Fraction as F
import json
from count_analysis import parse_counts, ratio, source_relative, analyze, WT, ANCESTOR

here=Path(__file__).resolve().parent
seq=json.loads((here/'source_contract.json').read_text())['sequences']
base={s:100+i for i,s in enumerate(seq)}
inputs=[{s:v*k for s,v in base.items()} for k in (1,3,7)]
outputs=[{s:v*k for s,v in base.items()} for k in range(2,11)]
for a in outputs:
    for b in inputs:
        assert all(ratio(a,b,s)==1 for s in seq)
source,norm=source_relative(outputs,inputs,seq)
assert all(v==1 for row in source.values() for v in row)
assert len(norm['TR1_aligned_totals'])==12

def raw(lib):
    return ('Genotype\tOccurrences\tMutationCount\n'+'\n'.join(f'{s}\t{v}\t{sum(a!=b for a,b in zip(s,ANCESTOR))}' for s,v in lib.items())+'\n').encode()
assert parse_counts(raw(base))==base
guards=0
for wrong in (raw(base).replace(b'Occurrences',b'Reads',1),raw(base)+raw(base).splitlines()[1]+b'\n',raw(base).replace(b'\t100\t',b'\t-1\t',1),raw(base).replace(b'\t100\t',b'\t1.2\t',1)):
    try:parse_counts(wrong)
    except ValueError:guards+=1
assert guards==4
try:parse_counts(raw(base),'0'*40)
except ValueError:guards+=1
assert guards==5
s=next(s for s in seq if s!=WT)
zero=dict(outputs[0]);zero[s]=0
assert ratio(zero,inputs[0],s)==0
zero_wt=dict(outputs[0]);zero_wt[WT]=0
assert ratio(zero_wt,inputs[0],s) is None
zero_in=dict(inputs[0]);zero_in[s]=0
assert ratio(outputs[0],zero_in,s) is None
missing=dict(inputs[0]);del missing[s]
assert ratio(outputs[0],missing,s) is None

# Source reproduction uses TR1's rows, not every row in each other file.
extra=dict(outputs[0]);extra[ANCESTOR]=100000
rep,norm2=source_relative([extra]+outputs[1:],inputs,seq)
assert norm2['TR1_aligned_totals'][3] != norm2['whole_file_totals'][3]
assert rep==source

# Median of normalized inputs is not generally median of their WT ratios.
special=[{WT:1,s:1,ANCESTOR:98},{WT:3,s:1,ANCESTOR:6},{WT:2,s:10,ANCESTOR:88}]
src,_=source_relative([{WT:5,s:10}]*9,special,[s,WT])
assert src[s][0]==F(2,5)
assert src[s][0] != sorted(ratio({WT:5,s:10},z,s) for z in special)[1]

result=analyze(outputs,inputs,seq)
assert result['finite_envelopes']==72 and result['primary_codon_null']['status']=='feasible'
assert result['conservative_barrier_pairs']==0 and all(x['bounded_pair_count']==0 for x in result['scenario_results'])
summary={'status':'PASS','scope':'Synthetic only; source outcomes unopened','guard_cases':guards,
         'ratio_checks':72*27,'full_pipeline_constant_fixture_vertices':72,
         'full_pipeline_constant_fixture_scenarios':27,'source_TR1_alignment_verified':True,
         'mean_and_median_normalization_distinctions_verified':True}
(here/'count_reader_verification.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary))
