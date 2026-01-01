"""Bounded exact variable elimination for single/paired-copy RNA constraints.

This is a small reference sum-product backend, not a new treewidth algorithm.
The cap is checked before allocating each intermediate Cartesian table.
"""
from itertools import product
from math import prod
from accessibility import CODONS, PAIRS, matching

def count_compatible(peptide,pairs,fixed=None,paired_copy=False,table_cap=1000000):
    pairs=matching(pairs,3*len(peptide));fixed=fixed or {}
    domains=[]
    for k,a in enumerate(peptide):
        single=[c for c in CODONS[a] if all(3*k+i not in fixed or fixed[3*k+i]==c[i] for i in range(3))]
        states=list(product(single,repeat=2)) if paired_copy else [(c,) for c in single]
        if paired_copy and a=='S':states=[s for s in states if (s[0][:2]=='AG')==(s[1][:2]=='AG')]
        for i,j in pairs:
            if i//3==j//3==k:
                states=[s for s in states if all(c[i%3]+c[j%3] in PAIRS for c in s) and (not paired_copy or (s[0][i%3] in 'AG')==(s[1][i%3] in 'AG'))]
        domains.append(states)
    if any(not d for d in domains):return {'count':'0','paired_copy':paired_copy,'largest_table':0,'induced_width':0}
    grouped={}
    for i,j in pairs:
        a,b=i//3,j//3
        if a!=b:grouped.setdefault((a,b),[]).append((i%3,j%3))
    factors=[]
    for k,domain in enumerate(domains):
        if len(domain)>table_cap:return {'count':None,'status':'resource_limit','stage':'initial_unary_factor','next_entries':len(domain),'table_cap':table_cap,'paired_copy':paired_copy}
        factors.append(((k,),{(v,):1 for v in range(len(domain))}))
    for (a,b),edges in grouped.items():
        if len(domains[a])*len(domains[b])>table_cap:return {'count':None,'status':'resource_limit','stage':'initial_pair_factor','next_entries':len(domains[a])*len(domains[b]),'table_cap':table_cap,'paired_copy':paired_copy}
        table={}
        for va,sa in enumerate(domains[a]):
            for vb,sb in enumerate(domains[b]):
                valid=all(all(sa[c][i]+sb[c][j] in PAIRS for c in range(len(sa))) and (not paired_copy or (sa[0][i] in 'AG')==(sa[1][i] in 'AG')) for i,j in edges)
                if valid:table[va,vb]=1
        factors.append(((a,b),table))
    remaining=set(range(len(domains)));largest=0;width=0;order=[]
    while remaining:
        def scope(v):return tuple(sorted(set().union(*(set(s) for s,t in factors if v in s))))
        variable=min(remaining,key=lambda v:(prod(len(domains[i]) for i in scope(v)),v))
        joint=scope(variable);allocation=prod(len(domains[i]) for i in joint)
        if allocation>table_cap:
            return {'count':None,'status':'resource_limit','next_entries':allocation,'table_cap':table_cap,'largest_table':largest,'induced_width':width,'order':order,'paired_copy':paired_copy}
        largest=max(largest,allocation);width=max(width,len(joint)-1)
        selected=[(s,t) for s,t in factors if variable in s]
        factors=[(s,t) for s,t in factors if variable not in s]
        indices={i:k for k,i in enumerate(joint)};projection=[tuple(indices[i] for i in s) for s,t in selected]
        output_scope=tuple(i for i in joint if i!=variable);output_indices=tuple(indices[i] for i in output_scope);table={}
        for values in product(*(range(len(domains[i])) for i in joint)):
            weight=1
            for ((s,t),ix) in zip(selected,projection):
                weight*=t.get(tuple(values[k] for k in ix),0)
                if not weight:break
            if weight:
                key=tuple(values[k] for k in output_indices);table[key]=table.get(key,0)+weight
        factors.append((output_scope,table));remaining.remove(variable);order.append(variable)
    result=prod(t.get((),0) for s,t in factors)
    return {'count':str(result),'status':'exact','paired_copy':paired_copy,'largest_table':largest,'induced_width':width,'order':order}
