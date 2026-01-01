"""Focused validation of the reusable general-pair-graph API and figure arrays."""
import itertools
import json
from pathlib import Path
import random
import sys
from graph_compiler import CODONS,PAIRS,valid,signature,endpoint_certificate,weighted_projection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'exact'))
from independent_oracle import graph,bfs

def main():
    if not __debug__:raise RuntimeError('assertions must remain enabled')
    rng=random.Random(2026102403);cases=0;endpoints=0
    for m in range(1,5):
        for _ in range(50):
            peptide=''.join(rng.choices(tuple(CODONS),k=m));initial=''.join(rng.choice(CODONS[a])for a in peptide)
            edges=[(i,j)for i in range(3*m)for j in range(i+1,3*m)if initial[i]+initial[j]in PAIRS and rng.random()<0.35]
            fixed={i:b for i,b in enumerate(initial)if rng.random()<0.2}
            vs=[''.join(cs)for cs in itertools.product(*(CODONS[a]for a in peptide))if valid(peptide,edges,''.join(cs),fixed)]
            adj=graph(vs);source=vs.index(initial);ds,counts=bfs(adj,source)
            model=weighted_projection(peptide,edges,initial,fixed)
            component=[v for v in vs if signature(v,edges,peptide)==signature(initial,edges,peptide)]
            assert int(model['volume'])==len(component)==len(ds)
            for target in rng.sample(range(len(vs)),min(3,len(vs))):
                cert=endpoint_certificate(peptide,edges,initial,vs[target],fixed)
                assert(cert['status']=='connected')==(target in ds)
                if target in ds:
                    assert cert['distance']==ds[target]
                    assert int(cert['shortest_paths'])==counts[target]
                endpoints+=1
            cases+=1
    # Genuine cyclic event graph, with independent cap statuses.
    p='KKFF';edges=[(i,j)for i in(2,5)for j in(8,11)];x='AAAAAAUUUUUU';y='AAGAAGUUCUUC'
    cert=endpoint_certificate(p,edges,x,y,subset_cap=0)
    assert cert['status']=='connected'and cert['distance']==4 and cert['count_status']=='resource_limit'and cert['shortest_paths']is None
    assert endpoint_certificate(p,edges,x,y)['shortest_paths']=='4'
    assert weighted_projection(p,edges,x,max_free=0)['status']=='resource_limit'
    assert endpoint_certificate('',[],'','')['shortest_paths']=='1'
    assert weighted_projection('',[],'')['volume']=='1'
    # Mutation-rule and composition counterexamples are direct standard-code facts.
    assert valid('S',[(0,2)],'UCG')and valid('S',[(0,2)],'AGU')
    assert signature('UCG',[(0,2)],'S')!=signature('AGU',[(0,2)],'S')
    a='AAGAAA';b='AAAAAG';assert valid('KK',[],a)and valid('KK',[],b)
    assert sum(c in'GC'for c in a)==sum(c in'GC'for c in b)==1
    # Fixed-GC=1 subspace has two vertices with no single-base edge.
    gc1=[''.join(cs)for cs in itertools.product(CODONS['K'],repeat=2)if sum(c in'GC'for c in''.join(cs))==1]
    assert set(gc1)=={a,b}and not any(graph(gc1))
    result={'seed':2026102403,'complete_genotype_domains_checked':cases,'endpoint_pairs_checked':endpoints,'status':'all assertions passed','boundary_cases':['empty peptide','cyclic dependency with count cap','separate volume cap','codon-step Ser signature change','fixed GC disconnected KK']}
    out=Path(__file__).resolve().parents[2]/'results'/'expansion';out.mkdir(parents=True,exist_ok=True)
    (out/'graph_compiler_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
