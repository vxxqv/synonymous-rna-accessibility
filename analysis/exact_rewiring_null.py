"""Exhaustive conditional pair-rewiring null for targets with <=6 pairs.

Preserves the exact paired endpoint set, pair count, observed compatibility,
variable-endpoint phase strata, fixed-context status, and broad span bins.
No MCMC or sampled disconnection is used. Larger targets are explicitly gated.
"""
from pathlib import Path
import sys,json,itertools,math
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from accessibility import *

def span_bin(d):return 0 if d<10 else 1 if d<30 else 2 if d<75 else 3

def targets(block,condition):
    pairs=block['touching_pairs'][condition];context=block['context'];v=block['variable_start0']
    def phase(i):return str((i-v)%3) if v<=i<v+15 else 'F'
    def category(i,j):return tuple(sorted((phase(i),phase(j))))+(span_bin(j-i),)
    counts=Counter(category(i,j) for i,j in pairs);endpoints=sorted(i for edge in pairs for i in edge)
    output=[]
    def recurse(available,remaining,current):
        if not available:
            output.append(tuple(current));return
        i=available[0]
        for j in available[1:]:
            cat=category(i,j)
            if j-i<4 or context[i]+context[j] not in PAIRS or not remaining[cat] or phase(i)==phase(j)=='F':continue
            remaining[cat]-=1
            recurse([k for k in available if k!=i and k!=j],remaining,current+[(i,j)])
            remaining[cat]+=1
    recurse(endpoints,counts,[])
    assert tuple(map(tuple,pairs)) in output
    return output

def evaluate(block,condition):
    pairings=targets(block,condition);seqs=enumerate_space(block['peptide'],cap=1024);v=block['variable_start0'];context=block['context']
    result=[]
    for pairs in pairings:
        keep=[s for s in seqs if compatible(context[:v]+s+context[v+15:],pairs)]
        internal=[(a-v,b-v) for a,b in pairs if v<=a<b<v+15]
        counts=Counter(signature(s,internal) for s in keep)
        result.append({'n':len(keep),'a':sum(n*n for n in counts.values())/len(keep)**2,'connection_numerator':str(sum(n*n for n in counts.values())),'connection_denominator':str(len(keep)**2),'pairings':pairs})
    native=json.loads((ROOT/'results/native'/f"{block['id']}_{condition}.json").read_text())
    return {'id':block['id'],'transcript_id':block['transcript_id'],'gene':block['gene'],'condition':condition,
        'null_population':len(result),'native_a':native['a_compatible'],'mean_null_a':sum(r['a'] for r in result)/len(result),
        'native_n':native['n_compatible'],'mean_null_n':sum(r['n'] for r in result)/len(result),'complete_null':result}

if __name__=='__main__':
    blocks=json.loads((ROOT/'data/native_windows.json').read_text());result=[];gated=0
    for block in blocks:
        for condition in ('incell','cellfree'):
            if len(block['touching_pairs'][condition])>6:gated+=1;continue
            result.append(evaluate(block,condition))
    out={'scope':'complete conditional null on <=6-pair targets; larger targets not inferred','max_pairs':6,'eligible':len(result),'gated':gated,
        'more_than_one_matching':sum(r['null_population']>1 for r in result),
        'metric_invariant_populations':sum(all(q['a']==r['native_a'] for q in r['complete_null']) for r in result),
        'interpretation':'The conditional null is nondiscriminating for connection probability if every matching in each population has the native value. Structural variation does not imply statistical information for this metric.', 'results':result}
    (ROOT/'results/exact_rewiring_null.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='results'}))
