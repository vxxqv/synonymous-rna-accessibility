"""Exact finite automaton for uniform-peptide, both-register feasibility.

Each input letter is one of 84 aligned amino-acid pairs for which C0 is
nonempty. Determinization stores the reachable subset of reversed right
codons in Cd. A nonempty subset lies within one amino acid's <=6 codons.
Thus the number of states is finite independently of arm length.
"""
from collections import Counter,deque
from functools import lru_cache
import itertools
import json
from pathlib import Path
from bridge_transfer import CODONS,block,within,across,minimum

ALPHABET=[(a,b) for a,b in itertools.product(CODONS,repeat=2) if block(a,b)[0]]
REV=sorted({c[::-1] for cs in CODONS.values() for c in cs})
INDEX={q:i for i,q in enumerate(REV)}

@lru_cache(None)
def advance(mask,a,b,d):
    result=0
    for l,q,hist,z,num in block(a,b)[1]:
        if not within(l,q,d):continue
        if mask is None or any(mask&(1<<i) and across(prev,l,d) for i,prev in enumerate(REV)):
            result|=1<<INDEX[q]
    return result

def automaton(d):
    initial=Counter(advance(None,a,b,d) for a,b in ALPHABET)
    initial.pop(0,None)
    discovered=set(initial); todo=deque(initial); matrix={}; failures={}
    while todo:
        state=todo.popleft()
        row=Counter(advance(state,a,b,d) for a,b in ALPHABET)
        failures[state]=row.pop(0,0)
        matrix[state]=dict(row)
        for nxt in row:
            if nxt not in discovered:
                discovered.add(nxt);todo.append(nxt)
    return dict(initial),matrix,failures

def tarjan(matrix):
    index={};low={};stack=[];on=set();scc=[]
    def visit(v):
        index[v]=low[v]=len(index);stack.append(v);on.add(v)
        for w in matrix[v]:
            if w not in index:visit(w);low[v]=min(low[v],low[w])
            elif w in on:low[v]=min(low[v],index[w])
        if low[v]==index[v]:
            group=[]
            while True:
                w=stack.pop();on.remove(w);group.append(w)
                if w==v:break
            scc.append(group)
    for v in matrix:
        if v not in index:visit(v)
    return scc

def root_bounds(matrix,component,steps=2000):
    """Collatz-Wielandt numerical bounds; roundoff is NOT interval-certified."""
    inside=set(component)
    if len(component)==1 and component[0] not in matrix[component[0]]:return [0.,0.]
    x={s:1. for s in component};lo=0.;hi=float('inf')
    for _ in range(steps):
        # Iterate A+I to remove periodic oscillation; rho(A+I)=rho(A)+1.
        y={s:x[s]+sum(c*x[t] for t,c in matrix[s].items() if t in inside) for s in component}
        ratios=[y[s]/x[s] for s in component]
        lo,hi=min(ratios)-1,max(ratios)-1
        norm=max(y.values());x={s:v/norm for s,v in y.items()}
        if hi-lo<1e-12:break
    return [lo,hi]

def integer_collatz(matrix,component,steps=200):
    """Exact PF enclosure, using positive integer test vectors only."""
    inside=set(component);scale=10**24;x={s:scale for s in component}
    for _ in range(steps):
        y={s:x[s]+sum(c*x[t] for t,c in matrix[s].items() if t in inside) for s in component}
        norm=max(y.values())
        x={s:max(1,v*scale//norm) for s,v in y.items()}
    y={s:sum(c*x[t] for t,c in matrix[s].items() if t in inside) for s in component}
    # Outward rounded rational endpoints with denominator 10^15.
    den=10**15
    lower=min(y[s]*den//x[s] for s in component)
    upper=max((y[s]*den+x[s]-1)//x[s] for s in component)
    assert all(lower*x[s]<=den*y[s]<=upper*x[s] for s in component)
    return {'denominator':str(den),'lower_numerator':str(lower),'upper_numerator':str(upper),'positive_integer_vector':{str(s):str(x[s]) for s in component},'inequality':'lower*x <= denominator*A*x <= upper*x, componentwise; all exact integers'}

def build(d,lengths=100):
    initial,matrix,failures=automaton(d)
    current=Counter(initial);counts=[]
    for m in range(1,lengths+1):
        counts.append({'m':m,'C0_nonempty':str(84**m),'both_nonempty':str(sum(current.values())),'all_peptide_pairs':str(400**m)})
        nxt=Counter()
        for s,num in current.items():
            for t,c in matrix[s].items():nxt[t]+=num*c
        current=nxt
    scc=tarjan(matrix)
    spectral=[]
    for component in scc:
        lo,hi=root_bounds(matrix,component)
        spectral.append({'size':len(component),'has_selfloop':any(s in matrix[s] for s in component),'rho_numerical_bounds':[lo,hi]})
    spectral.sort(key=lambda z:z['rho_numerical_bounds'][1],reverse=True)
    assert len(scc)==1 and any(s in matrix[s] for s in scc[0])
    return {'shift':d,'alphabet_size':84,'states':len(matrix),'state_upper_bound_sum_nonempty_codon_subsets':sum(2**len(cs)-1 for cs in CODONS.values()),'initial':{str(k):v for k,v in initial.items()},'transition_matrix':{str(s):{str(t):c for t,c in row.items()} for s,row in matrix.items()},'failure_counts':{str(s):v for s,v in failures.items()},'strong_components':spectral,'exact_pf_certificate':integer_collatz(matrix,scc[0]),'exact_counts':counts}

def main():
    root=Path(__file__).parent
    results={str(d):build(d) for d in (1,2,3)}
    exact=json.loads((root/'exact_peptide_ensemble.json').read_text())
    for d in (1,2,3):
        for m in (1,2):
            assert int(results[str(d)]['exact_counts'][m-1]['both_nonempty'])==exact[f'm{m}_d{d}']['both_feasible']
        # A third exact brute-force enumeration over 84^3 input words is
        # independent of automaton matrix powering, but shares local semantics.
        # Check all words for m<=2 directly against min-plus transfer above.
        for chunks in itertools.product(ALPHABET,repeat=2):
            mask=None
            for a,b in chunks:mask=advance(mask,a,b,d)
            l=''.join(a for a,b in chunks);r=''.join(b for a,b in chunks[::-1])
            delta,_,_,_=minimum(l,r,d)
            assert bool(mask)==(delta is not None)
    (root/'exact_peptide_automata.json').write_text(json.dumps(results,indent=2)+'\n')
    for d,r in results.items():
        print('shift',d,'states',r['states'],'top SCC',r['strong_components'][:3])
        for i in [0,1,2,3,5,11,19,49,99]:print(r['exact_counts'][i])

if __name__=='__main__':main()
