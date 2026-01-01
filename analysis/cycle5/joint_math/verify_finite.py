"""Exact bounded fixtures for the joint-trait proof; not a universal proof."""
import sys
if not __debug__:
    raise RuntimeError('Disabled assertions refused before reads or outputs')
sys.dont_write_bytecode=True
import argparse
from collections import deque
from fractions import Fraction as Q
import hashlib
from itertools import product, permutations
import json
from pathlib import Path
import platform
import time

CODONS={'F':'UUU UUC','L':'UUA UUG CUU CUC CUA CUG','I':'AUU AUC AUA','M':'AUG',
 'V':'GUU GUC GUA GUG','S':'UCU UCC UCA UCG AGU AGC','P':'CCU CCC CCA CCG',
 'T':'ACU ACC ACA ACG','A':'GCU GCC GCA GCG','Y':'UAU UAC','H':'CAU CAC',
 'Q':'CAA CAG','N':'AAU AAC','K':'AAA AAG','D':'GAU GAC','E':'GAA GAG',
 'C':'UGU UGC','W':'UGG','R':'CGU CGC CGA CGG AGA AGG','G':'GGU GGC GGA GGG'}
CODONS={key:tuple(value.split()) for key,value in CODONS.items()}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def encode(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x
def graph(vertices,adjacent):
    return {v:tuple(w for w in vertices if w!=v and adjacent(v,w)) for v in vertices}
def hamming(a,b):return sum(x!=y for x,y in zip(a,b))
def cube(n):return graph(tuple(product((0,1),repeat=n)),lambda a,b:hamming(a,b)==1)
def codon_graph(aa):return graph(CODONS[aa],lambda a,b:hamming(a,b)==1)
def cost(value,box):
    return max([Q(0)]+[bound-z for z,(bound,_) in zip(value,box) if bound is not None]+
               [z-bound for z,(_,bound) in zip(value,box) if bound is not None])
def reachable(g,s,t,allowed):
    if s not in allowed or t not in allowed:return None
    previous={s:None};todo=deque([s])
    while todo:
        v=todo.popleft()
        if v==t:
            result=[]
            while v is not None:result.append(v);v=previous[v]
            return list(reversed(result))
        for w in g[v]:
            if w in allowed and w not in previous:previous[w]=v;todo.append(w)
    return None
def threshold_oracle(g,s,t,costs):
    for level in sorted(set(costs.values())):
        path=reachable(g,s,t,{v for v in g if costs[v]<=level})
        if path is not None:return level,path
    return None,None
def floyd_oracle(g,s,t,costs):
    vertices=list(g); inf=None
    d={(v,w):(costs[v] if v==w else max(costs[v],costs[w]) if w in g[v] else inf) for v in vertices for w in vertices}
    for k in vertices:
        for v in vertices:
            for w in vertices:
                if d[v,k] is not None and d[k,w] is not None:
                    candidate=max(d[v,k],d[k,w])
                    if d[v,w] is None or candidate<d[v,w]:d[v,w]=candidate
    return d[s,t]
def rank(vectors):
    a=[[Q(x) for x in row] for row in vectors]
    if not a:return 0
    r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];factor=a[r][col]
        a[r]=[x/factor for x in a[r]]
        for j in range(len(a)):
            if j!=r:
                factor=a[j][col];a[j]=[x-factor*y for x,y in zip(a[j],a[r])]
        r+=1
        if r==len(a):break
    return r
def segment_distance(point,total):
    lines=[(sign*Q(p),-sign*Q(d)) for p,d in zip(point,total) for sign in (1,-1)]
    candidates={Q(0),Q(1)}
    for i,(a,b) in enumerate(lines):
        for c,d in lines[i+1:]:
            if b!=d:
                value=(c-a)/(b-d)
                if 0<=value<=1:candidates.add(value)
    return min(max(a+b*t for a,b in lines) for t in candidates)
def value_at(vectors,x):return tuple(sum(Q(v[k])*bit for v,bit in zip(vectors,x)) for k in range(len(vectors[0])))
def check_fixture(name,g,s,t,values,box,expected):
    costs={v:cost(values[v],box) for v in g}
    result,path=threshold_oracle(g,s,t,costs)
    assert result==expected,(name,result,expected)
    assert floyd_oracle(g,s,t,costs)==expected
    if path:
        assert path[0]==s and path[-1]==t and max(costs[v] for v in path)==expected
        assert len(path)==len(set(path)) and all(w in g[v] for v,w in zip(path,path[1:]))
    return {'name':name,'vertices':len(g),'edges':sum(map(len,g.values()))//2,'expected_R':expected,'observed_R':result,'path':path,'values':values,'box':box}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    out=Path(args.out)
    if out.exists():raise FileExistsError('Fresh output required')
    start=time.perf_counter(); here=Path(__file__).resolve().parent
    inputs=[Path(__file__).resolve(),here/'claim_and_proof.md',here/'prior_map.md']
    before={str(p.name):sha(p) for p in inputs}
    fixtures=[]
    # Standard-code partition and all intra-amino-acid edges from direct Hamming definition.
    assert len(CODONS)==20 and len({c for cs in CODONS.values() for c in cs})==61
    component_counts={}
    for aa in CODONS:
        g=codon_graph(aa);left=set(g);count=0
        while left:
            root=min(left);seen={v for v in g if reachable(g,root,v,set(g))};left-=seen;count+=1
        component_counts[aa]=count
        assert count==(2 if aa=='S' else 1)
    l=codon_graph('L');s,t='UUA','CUC'
    values={v:(Q(0) if v in (s,t) else Q(-3),) for v in l}
    fixtures.append(check_fixture('Leu local eta essential at rank0',l,s,t,values,[(Q(0),None)],Q(3)))
    ser=codon_graph('S');values={v:(Q(0),) for v in ser}
    fixtures.append(check_fixture('Ser different base components',ser,'UCU','AGU',values,[(Q(0),None)],None))
    g=cube(2);vectors=((1,-1),(-1,1));values={x:value_at(vectors,x) for x in g}
    fixtures.append(check_fixture('rank1 sharp lower-bound conflict',g,(0,0),(1,1),values,[(Q(0),None)]*2,Q(1)))
    # Each individual trait remains connected between these endpoints.
    assert all(reachable(g,(0,0),(1,1),{v for v in g if values[v][k]>=0}) for k in range(2))
    rational=((Q(1,2),Q(-1,3)),(Q(-1,2),Q(1,3)))
    values={x:value_at(rational,x) for x in g}
    fixtures.append(check_fixture('rational scaled effects',g,(0,0),(1,1),values,[(Q(0),None)]*2,Q(1,3)))
    g3=cube(3);vectors=((1,-1),(-1,1),(1,1));values={x:value_at(vectors,x) for x in g3}
    fixtures.append(check_fixture('background detour beats endpoint orders',g3,(0,0,0),(1,1,0),values,[(Q(0),None)]*2,Q(0)))
    assert all(max(cost(values[x],[(Q(0),None)]*2) for x in [(0,0,0),mid,(1,1,0)])==1 for mid in [(1,0,0),(0,1,0)])
    independent=[x for x in product((0,1),repeat=4) if not(any(x[:2]) and any(x[2:]))]
    coupled=graph(independent,lambda a,b:hamming(a,b)==1);values={x:(Q(sum(x)),) for x in coupled}
    fixtures.append(check_fixture('coupling invalidates product bound',coupled,(1,1,0,0),(0,0,1,1),values,[(Q(2),None)],Q(2)))
    values={v:(Q(1) if v in (s,t) else Q(0),) for v in l}
    fixtures.append(check_fixture('weak plateau versus strict',l,s,t,values,[(Q(0),None)],Q(0)))
    assert reachable(l,s,t,{v for v in l if values[v][0]>0}) is None
    values={x:(Q(sum(x)),) for x in g}
    fixtures.append(check_fixture('interval equality defeats alignment',g,(1,0),(0,1),values,[(Q(1),Q(1))],Q(1)))
    fixtures.append(check_fixture('coincident endpoints',g,(1,0),(1,0),values,[(Q(1),Q(1))],Q(0)))
    eps=Q(1,4);model={v:(-eps,) for v in l};true={v:(Q(0) if v in (s,t) else -2*eps,) for v in l}
    assert max(abs(true[v][0]-model[v][0]) for v in l)==eps
    fixtures.append(check_fixture('uniform-error coefficient2 sharp',l,s,t,true,[(Q(0),None)],2*eps))
    # Exhaustive finite vector family: all 729 ordered 2x3 matrices over {-1,0,1}.
    vectors2=list(product((-1,0,1),repeat=2)); cases=0;worst_ratio=Q(0);witness_examples=[]
    for vs in product(vectors2,repeat=3):
        total=tuple(sum(v[k] for v in vs) for k in range(2))
        bound=rank(vs)*max(abs(x) for v in vs for x in v)
        best=None;best_order=None
        for order in permutations(range(3)):
            partial=[0,0];radius=Q(0)
            for i in order:
                partial=[a+b for a,b in zip(partial,vs[i])]
                radius=max(radius,segment_distance(partial,total))
            if best is None or radius<best:best,best_order=radius,order
        assert best<=bound,(vs,best,bound)
        vals={x:value_at(vs,x) for x in g3};box=[(Q(min(0,d)),Q(max(0,d))) for d in total]
        costs={x:cost(vals[x],box) for x in g3}
        observed,path=threshold_oracle(g3,(0,0,0),(1,1,1),costs)
        assert observed<=bound and floyd_oracle(g3,(0,0,0),(1,1,1),costs)==observed
        if bound:worst_ratio=max(worst_ratio,observed/bound)
        if best>0 and any(d>0 for d in total) and any(d<0 for d in total) and not witness_examples:
            witness_examples.append({'vectors':vs,'total':total,'order':best_order,'best_segment_radius':best,'bound':bound,'exact_R':observed})
        cases+=1
    # Spectral warning: induced Laplacian is not the principal base Laplacian.
    vertices=list(g);allowed=[(0,0),(1,1)]
    principal=[[Q(len(g[v])) if v==w else Q(-1) if w in g[v] else Q(0) for w in allowed] for v in allowed]
    induced=[[Q(sum(w in allowed for w in g[v])) if v==w else Q(-1) if w in g[v] else Q(0) for w in allowed] for v in allowed]
    assert principal==[[2,0],[0,2]] and induced==[[0,0],[0,0]]
    assert {str(p.name):sha(p) for p in inputs}==before
    out.mkdir(parents=True)
    result={'status':'passed','fixtures':fixtures,'standard_code_component_counts':component_counts,
      'exhaustive_vector_contract':'All 729 ordered length3 lists of vectors in {-1,0,1}^2; all six orders tried; exact interval bottleneck checked by threshold BFS and minmax Floyd.',
      'exhaustive_vector_matrices':cases,'worst_R_over_rank_M_in_tested_range':worst_ratio,
      'mixed_sign_total_fixture':witness_examples,'spectral_counterexample':{'principal':principal,'induced':induced},
      'assertion_tested':'Finite rank/segment bound and listed counterexamples only, not proof of the universal theorem.',
      'no_biological_outcomes_read':True}
    path=out/'results.json';path.write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    manifest={'run':{'status':'completed','elapsed_seconds':time.perf_counter()-start,'argv':sys.argv},
      'software':{'python':platform.python_version(),'implementation':platform.python_implementation()},
      'arithmetic':'fractions.Fraction exact rational; no numerical tolerance or randomness',
      'inputs_before':before,'inputs_after':{str(p.name):sha(p) for p in inputs},'outputs':{'results.json':sha(path)},
      'review_provenance':'Author implementation with two distinct path algorithms; fresh independent reviewer pending.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','fixtures':len(fixtures),'exhaustive_vector_matrices':cases,'elapsed_seconds':manifest['run']['elapsed_seconds'],'results_sha256':sha(path)},indent=2))
if __name__=='__main__':main()
