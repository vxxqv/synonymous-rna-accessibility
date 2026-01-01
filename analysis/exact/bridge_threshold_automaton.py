"""Determinized min-plus automata counting peptide pairs with delta <=1 exactly.

State is two disjoint bitsets of Q-codons attainable at minimum Hamming cost
zero/one. Costs >1 cannot return below one, so dropping them is exact. Every
alphabet letter is an aligned amino-acid pair with nonempty C0 (84 total).
"""
from collections import Counter,deque
from functools import lru_cache
import itertools
import json
from pathlib import Path
from peptide_automaton import ALPHABET,REV,INDEX,tarjan,root_bounds,integer_collatz
from bridge_transfer import block,within,across,minimum

@lru_cache(None)
def advance_threshold(state,a,b,d):
    zero=0;one=0
    for l,q,hist,cost,num in block(a,b)[1]:
        if cost>1 or not within(l,q,d):continue
        bit=1<<INDEX[q]
        if state is None:
            if cost==0:zero|=bit
            else:one|=bit
            continue
        prev0,prev1=state
        yes0=any(prev0&(1<<i) and across(prev,l,d) for i,prev in enumerate(REV))
        if yes0:
            if cost==0:zero|=bit
            else:one|=bit
        if cost==0 and any(prev1&(1<<i) and across(prev,l,d) for i,prev in enumerate(REV)):
            one|=bit
    return zero,one&~zero

def build(d,lengths=100):
    initial=Counter(advance_threshold(None,a,b,d) for a,b in ALPHABET)
    initial.pop((0,0),None)
    seen=set(initial);todo=deque(initial);matrix={};failures={}
    while todo:
        s=todo.popleft();row=Counter(advance_threshold(s,a,b,d) for a,b in ALPHABET)
        failures[s]=row.pop((0,0),0);matrix[s]=dict(row)
        for t in row:
            if t not in seen:seen.add(t);todo.append(t)
    dp=Counter(initial);counts=[]
    for m in range(1,lengths+1):
        zero=sum(c for s,c in dp.items() if s[0]);one=sum(c for s,c in dp.items() if not s[0])
        counts.append({'m':m,'delta0_peptide_pairs':str(zero),'delta1_peptide_pairs':str(one),'delta_le1_peptide_pairs':str(zero+one),'C0_feasible':str(84**m),'all_peptide_pairs':str(400**m)})
        nxt=Counter()
        for s,c in dp.items():
            for t,k in matrix[s].items():nxt[t]+=c*k
        dp=nxt
    components=tarjan(matrix);spectral=[]
    for cc in components:
        if len(cc)==1 and cc[0] not in matrix[cc[0]]:continue
        spectral.append({'size':len(cc),'rho_numerical':root_bounds(matrix,cc),'exact_pf_certificate':integer_collatz(matrix,cc),'has_selfloop':any(s in matrix[s] for s in cc),'all_states_zero_cost':all(s[0] for s in cc),'all_states_positive_cost':all(not s[0] for s in cc)})
    spectral.sort(key=lambda z:z['rho_numerical'][1],reverse=True)
    # JSON string representation of tuple states; this is a transparent fixture.
    return {'shift':d,'states':len(matrix),'nontrivial_SCCs':spectral,'initial':{str(s):c for s,c in initial.items()},'matrix':{str(s):{str(t):c for t,c in row.items()} for s,row in matrix.items()},'exact_counts':counts}

def main():
    results={str(d):build(d) for d in (1,2,3)}
    exact=json.loads((Path(__file__).parent/'exact_peptide_ensemble.json').read_text())
    for d in (1,2,3):
        for m in (1,2):
            z=results[str(d)]['exact_counts'][m-1];expected=exact[f'm{m}_d{d}']['delta_counts']
            assert int(z['delta0_peptide_pairs'])==expected.get('0',0)
            assert int(z['delta1_peptide_pairs'])==expected.get('1',0)
        for chunks in itertools.product(ALPHABET,repeat=2):
            s=None
            for a,b in chunks:s=advance_threshold(s,a,b,d)
            l=''.join(a for a,b in chunks);r=''.join(b for a,b in chunks[::-1])
            delta,_,_,_=minimum(l,r,d)
            assert bool(s[0])==(delta==0)
            assert bool(s[0] or s[1])==(delta is not None and delta<=1)
    # Exact recurrence certificate: finite automaton dimension 42 bounds the
    # scalar residual recurrence by Cayley-Hamilton. Checking >=42 residual
    # terms establishes the identities for all m>=1, not just m<=100.
    r=results['1'];assert r['states']==42
    for z in r['exact_counts']:
        m=z['m']
        assert int(z['delta0_peptide_pairs'])==24*9**(m-1)
        assert int(z['delta1_peptide_pairs'])==8*(3*m-1)*9**(m-1)
    r['closed_form_certificate']={'dimension':42,'checked_terms':100,'delta0':'24*9^(m-1)','delta1':'8*(3*m-1)*9^(m-1)','delta_le1':'8*(3*m+2)*9^(m-1)','proof':'The exact 42-state linear representation implies by Cayley-Hamilton that first42 zero terms of residual (E-9)N0 or (E-9)^2 N1 force all residual terms zero. Initial terms fix the formulas. The 100 exact terms supply more than required42 residual terms.'}
    (Path(__file__).parent/'exact_bridge_threshold_automata.json').write_text(json.dumps(results,indent=2)+'\n')
    for d,r in results.items():
        print('shift',d,'states',r['states'],'SCC summaries',[{k:v for k,v in s.items() if k!='exact_pf_certificate'} for s in r['nontrivial_SCCs']])
        for i in (0,1,2,3,5,11,19,49,99):print(r['exact_counts'][i])

if __name__=='__main__':main()
