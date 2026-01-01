"""Exact codon-block transfer for two-arm shifted-register Hamming enumerators.

No main-project code imported. Supports shifts 1,2,3. Polynomial coefficients,
minimum costs, and degeneracies use integers. Codons in Q are reversed RNA codons.
"""
from collections import Counter,defaultdict
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import random
import time
from independent_oracle import CODONS,PAIRS,vertices,register,distance

def ham(a,b): return sum(x!=y for x,y in zip(a,b))

@lru_cache(None)
def block(a,b):
    left=CODONS[a]; right=tuple(c[::-1] for c in CODONS[b])
    target=[(u,v) for u in left for v in right if all(u[j]+v[j] in PAIRS for j in range(3))]
    data=[]
    for l in left:
        for q in right:
            hist=Counter(ham(u,l)+ham(v,q) for u,v in target)
            if hist:
                z=min(hist)
                data.append((l,q,tuple(sorted(hist.items())),z,hist[z]))
    return len(target),tuple(data)

@lru_cache(None)
def within(l,q,d):
    return all(l[j]+q[j-d] in PAIRS for j in range(d,3))

@lru_cache(None)
def across(qprev,l,d):
    return all(l[j]+qprev[3+j-d] in PAIRS for j in range(d))

def states(left,right,d):
    assert d in (1,2,3) and len(left)==len(right) and len(left)>0
    return [tuple(z for z in block(a,b)[1] if within(z[0],z[1],d)) for a,b in zip(left,right[::-1])]

def shifted_count(left,right,d):
    """Count Cd even when C0 is empty; no target compatibility filtering."""
    assert d in (1,2,3) and len(left)==len(right) and len(left)>0
    dp=None
    for a,b in zip(left,right[::-1]):
        nxt={}
        for l in CODONS[a]:
            for codon in CODONS[b]:
                q=codon[::-1]
                if not within(l,q,d):continue
                if dp is None:nxt[q]=nxt.get(q,0)+1
                else:
                    for prev,c in dp.items():
                        if across(prev,l,d):nxt[q]=nxt.get(q,0)+c
        dp=nxt
    return sum(dp.values())

def minimum(left,right,d):
    """Return (delta, minimizer-pair count, |C0|, |Cd|); delta=None if empty."""
    bs=states(left,right,d)
    c0=math.prod(block(a,b)[0] for a,b in zip(left,right[::-1]))
    if not c0: return None,0,0,shifted_count(left,right,d)
    dp={}; counts={}
    for l,q,hist,z,num in bs[0]:
        counts[q]=counts.get(q,0)+1
        old=dp.get(q)
        if old is None or z<old[0]: dp[q]=(z,num)
        elif z==old[0]: dp[q]=(z,old[1]+num)
    for choices in bs[1:]:
        nxt={}; nc={}
        for l,q,hist,z,num in choices:
            for prev,(cost,ways) in dp.items():
                if not across(prev,l,d): continue
                new=cost+z; waysnew=ways*num; old=nxt.get(q)
                if old is None or new<old[0]: nxt[q]=(new,waysnew)
                elif new==old[0]: nxt[q]=(new,old[1]+waysnew)
                nc[q]=nc.get(q,0)+counts[prev]
        dp,counts=nxt,nc
    if not dp: return None,0,c0,0
    val=min(c for c,k in dp.values())
    return val,sum(k for c,k in dp.values() if c==val),c0,sum(counts.values())

def enumerator(left,right,d):
    bs=states(left,right,d)
    dp={}
    for l,q,hist,z,num in bs[0]:
        for k,v in hist: dp.setdefault(q,Counter())[k]+=v
    for choices in bs[1:]:
        nxt={}
        for l,q,hist,z,num in choices:
            for prev,poly in dp.items():
                if not across(prev,l,d): continue
                target=nxt.setdefault(q,Counter())
                for i,v in poly.items():
                    for j,w in hist: target[i+j]+=v*w
        dp=nxt
    result=Counter()
    for poly in dp.values(): result.update(poly)
    return dict(sorted(result.items()))

def brute(left,right,d):
    n=3*len(left); peptide=tuple(left)+tuple(right)
    a=vertices(peptide,register(n,0)); b=vertices(peptide,register(n,d))
    hist=Counter(distance(x,y) for x in a for y in b)
    return dict(sorted(hist.items())),len(a),len(b)

def validate():
    counts=Counter(); rng=random.Random(20261024)
    pairs=list(itertools.product(CODONS,repeat=2))
    # Exhaustive 400 peptide pairs and all three shifts for n=3.
    cases=[(a,b) for a,b in pairs]
    # All constant two-codon arms plus 300 unrestricted random four-aa peptides.
    cases += [(a*2,b*2) for a,b in pairs]
    cases += [(''.join(rng.choices(tuple(CODONS),k=2)),''.join(rng.choices(tuple(CODONS),k=2))) for _ in range(300)]
    # All all-genotype fibre fixtures provide rich intersections, including Gamma4.
    examples=json.loads((Path(__file__).parent/'register_fibonacci_examples.json').read_text())
    for d in ('1','2','3'):
        for fixture in examples[d].values():
            p=fixture['peptide']; cases.append((p[:2],p[2:]))
    for l,r in cases:
        for d in (1,2,3):
            expected,c0,cd=brute(l,r,d)
            actual=enumerator(l,r,d)
            assert actual==expected,(l,r,d,expected,actual)
            delta,ties,c0b,cdb=minimum(l,r,d)
            assert c0b==c0
            assert cdb==cd
            assert delta==min(expected,default=None)
            assert ties==expected.get(delta,0)
            counts['cases']+=1;counts['ordered_endpoint_pairs']+=sum(expected.values())
    return dict(counts)

def exact_short_ensemble():
    """Uniform peptides; enumerate all C0-feasible pairs, plus exact excluded counts."""
    compatible=[(a,b) for a,b in itertools.product(CODONS,repeat=2) if block(a,b)[0]]
    assert len(compatible)==84
    output={}
    for m in (1,2):
        for d in (1,2,3):
            hist=Counter(); zeros=0; totalpairs=0
            for chunks in itertools.product(compatible,repeat=m):
                l=''.join(a for a,b in chunks);r=''.join(b for a,b in chunks[::-1])
                delta,ties,c0,cd=minimum(l,r,d)
                if delta is not None: hist[delta]+=1
                else: zeros+=1
                totalpairs+=1
            output[f'm{m}_d{d}']={'peptide_pairs_total':20**(2*m),'C0_feasible':totalpairs,'both_feasible':sum(hist.values()),'Cd_empty_given_C0':zeros,'delta_counts':dict(sorted(hist.items()))}
    return output

def sampled_ensemble(samples=4000,maxm=12):
    """IID uniform peptide pairs conditional on C0 nonempty, no outcome selection.

    C0 feasibility factorizes into 84 aligned peptide-pair blocks. Uniform iid
    selection of those blocks is exact conditional sampling. Retain every draw;
    Cd failure is recorded, and delta distribution conditions on both nonempty.
    """
    seed=2026102401;rng=random.Random(seed)
    compatible=[(a,b) for a,b in itertools.product(CODONS,repeat=2) if block(a,b)[0]]
    output=[]
    for m in range(3,maxm+1):
        examples=[tuple(rng.choice(compatible) for _ in range(m)) for _ in range(samples)]
        for d in (1,2,3):
            hist=Counter(); minima=Counter()
            for chunks in examples:
                l=''.join(a for a,b in chunks);r=''.join(b for a,b in chunks[::-1])
                delta,ties,c0,cd=minimum(l,r,d)
                if delta is not None:
                    hist[delta]+=1
                    minima['pairs_at_minimum']+=ties
            accepted=sum(hist.values())
            # Wilson score interval for Cd feasible conditional on C0, z=1.96.
            p=accepted/samples;den=1+1.96**2/samples
            ctr=(p+1.96**2/(2*samples))/den
            radius=1.96*math.sqrt(p*(1-p)/samples+1.96**2/(4*samples**2))/den
            output.append({'m':m,'shift':d,'draws_given_C0':samples,'both_feasible':accepted,'delta_counts':dict(sorted(hist.items())),'both_given_C0_fraction':p,'both_given_C0_wilson95':[ctr-radius,ctr+radius],'C0_unconditional_probability':(84/400)**m})
    return {'seed':seed,'distribution':'uniform over all ordered peptide pairs conditional C0 nonempty; same draws reused for shifts 1,2,3 at a given m','results':output}

def main():
    root=Path(__file__).parent;t=time.perf_counter()
    valid=validate();print('validation',valid,flush=True)
    (root/'transfer_validation.json').write_text(json.dumps(valid,indent=2)+'\n')
    exact=exact_short_ensemble();print('exact ensemble',exact,flush=True)
    (root/'exact_peptide_ensemble.json').write_text(json.dumps(exact,indent=2)+'\n')
    # Exploratory sampler retained for transparency, but exact automata now
    # supersede it; do not spend compute on sampling in the standard workflow.
    print('exact validation saved; seconds',time.perf_counter()-t,flush=True)

if __name__=='__main__':main()
