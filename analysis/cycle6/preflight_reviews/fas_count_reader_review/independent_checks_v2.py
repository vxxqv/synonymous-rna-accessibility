"""Synthetic-only independent count contract references, no real source data."""
from pathlib import Path
from fractions import Fraction as F
import collections
import hashlib
import importlib.util
import json
import random
import sys
from statistics import median
from sympy import Matrix

HERE = Path(__file__).resolve().parent
# Import the held candidate snapshot, never observed data or mutable work products.
spec = importlib.util.spec_from_file_location('ordinal_geometry', HERE/'corrected_ordinal_geometry.py')
geometry = importlib.util.module_from_spec(spec);sys.modules['ordinal_geometry']=geometry;spec.loader.exec_module(geometry)
spec = importlib.util.spec_from_file_location('held_count_analysis', HERE/'corrected_count_analysis.py')
candidate = importlib.util.module_from_spec(spec);spec.loader.exec_module(candidate)
config=json.loads((HERE/'corrected_source_contract.json').read_text(encoding='utf8'))
seq=config['sequences'];wt=candidate.WT;v=next(s for s in seq if s!=wt)
rng=random.Random(612730)

def reference_source(outputs, inputs, selected):
    universe=list(inputs[0])
    totals=[sum(lib.get(s,0) for s in universe) for lib in inputs]
    result={s:[] for s in selected}
    def input_med(s):
        vals=[F(lib[s],total) for lib,total in zip(inputs,totals) if s in lib and total>0]
        return median(vals) if vals else None
    for out in outputs:
        for s in selected:
            left=input_med(s);right=input_med(wt)
            if s not in universe or wt not in universe or s not in out or wt not in out or out[wt]<=0 or left is None or right is None or left<=0 or right<=0:
                result[s].append(None)
            else:
                # Form both normalized enrichments explicitly before dividing.
                total=sum(out.get(x,0) for x in universe)
                if total<=0:result[s].append(None)
                else:result[s].append((F(out[s],total)/left)/(F(out[wt],total)/right))
    return result

checks=collections.Counter();failures=[]
for trial in range(20):
    # Include ancestral DNA beyond the fixed human-protein class; totals are source-level.
    inp=[{s:rng.randrange(1,300) for s in seq+[candidate.ANCESTOR]} for _ in range(3)]
    out=[{s:rng.randrange(0,300) for s in seq+[candidate.ANCESTOR]} for _ in range(9)]
    for z in out:z[wt]=rng.randrange(1,300)
    if trial%3==0:del inp[1][v]
    if trial%4==0:del out[4][v]
    if trial%5==0:inp[0].pop(candidate.ANCESTOR)
    expected=reference_source(out,inp,seq)
    got,_=candidate.source_relative(out,inp,seq)
    if got!=expected:failures.append({'scope':'random normalized-median source reproduction','trial':trial})
    checks['source_output_vertex_comparisons']+=72*9
    for b in out:
        for t in inp:
            for s in seq:
                expected_ratio=(F(b[s]*t[wt],b[wt]*t[s]) if s in b and s in t and t[s]>0 and t[wt]>0 else None)
                if candidate.ratio(b,t,s)!=expected_ratio:failures.append({'scope':'relative ratio','trial':trial})
                checks['relative_ratio_comparisons']+=1

# Independently build treatment-coded features rather than call author's feature helper.
def independent_features(rows):
    levels=[sorted({s[i] for s in rows}) for i in range(len(rows[0]))]
    return Matrix([[1]+[int(row[k]==a) for k,vals in enumerate(levels) for a in vals[1:]] for row in rows])
codons=[tuple(s[k:k+3] for k in (23,38,29,32,41)) for s in seq]
positions=[i for i in range(63) if len({s[i] for s in seq})>1]
nts=[tuple(s[i] for i in positions) for s in seq]
cx=independent_features(codons);nx=independent_features(nts)
span={'codon_rank':cx.rank(),'nucleotide_rank':nx.rank(),'joint_rank':cx.row_join(nx).rank(),
      'same_column_space':cx.rank()==nx.rank()==cx.row_join(nx).rank(),
      'variable_nt_positions1':[k+1 for k in positions]}
if not span['same_column_space']:failures.append({'scope':'expected design equivalence failed'})

# Valid nonempty partial libraries: missing TR1-aligned input universe must not
# erase primary unknown-state results. No source numeric data are involved.
inputs=[{wt:1,v:1},{candidate.ANCESTOR:1},{wt:1,v:1}]
outputs=[{wt:2,v:3} for _ in range(9)]
expected=reference_source(outputs,inputs,[wt,v])
try:
    result,_=candidate.source_relative(outputs,inputs,[wt,v])
    degenerate={'returned':True,'matches_NA_omission_reference':result==expected}
    if result!=expected:failures.append({'scope':'corrected normalization mismatches independent degenerate reference'})
except Exception as ex:
    degenerate={'returned':False,'exception_type':type(ex).__name__,'message':str(ex),
                'expected_source_relative_by_output':{k:list(map(str,vals)) for k,vals in expected.items()},
                'primary_TR2_ratios_all_unknown':all(candidate.ratio(o,inputs[1],s) is None for o in outputs for s in seq)}
    failures.append({'scope':'valid missing TR1-aligned input aborts analysis','details':degenerate})

# Zero reference behavior for ratio() is correct; all affected vertices unknown.
out=dict.fromkeys(seq,5);inp=dict.fromkeys(seq,7)
for target,which in [(out,'zero_reference_output'),(inp,'zero_reference_input')]:
    original=target[wt];target[wt]=0
    if not all(candidate.ratio(out,inp,s) is None for s in seq):failures.append({'scope':which})
    checks[which]=72;target[wt]=original
out[v]=0
if candidate.ratio(out,inp,v)!=0:failures.append({'scope':'valid output zero'})
checks['zero_output']=1

receipt={'adapter_history':'Copy of held independent verifier, input snapshot names and output filename changed; added explicit mismatch failure on newly reachable degenerate success branch.', 'scope':'independent synthetic/reference checks only; no observed count values',
         'candidate_sha256':hashlib.sha256((HERE/'corrected_count_analysis.py').read_bytes()).hexdigest(),
         'checks':dict(checks),'design_column_space':span,
         'degenerate_TR1_alignment':degenerate,'failures':failures,
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'independent_checks_v2.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(receipt,indent=2))
