"""Exact projection/entropic-activity audit for arbitrary simultaneous pair graphs."""
from collections import Counter,defaultdict
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'exact'))
from independent_oracle import CODE,CODONS,PAIRS,RY
from arbitrary_pair_graph_oracle import signature,star_fixture,forest_path_count

def volume_certificate(peptide,edges,sample,fixed=None):
    fixed=fixed or {};paired={p for e in edges for p in e};domains={};codon_domains=[]
    for c,a in enumerate(peptide):
        cs=[z for z in CODONS[a] if all((3*c+j not in paired or RY[z[j]]==RY[sample[3*c+j]]) and
             (3*c+j not in fixed or z[j]==fixed[3*c+j])for j in range(3)) and
             (a!='S' or (z[:2]=='AG')==(sample[3*c:3*c+2]=='AG'))]
        assert cs
        codon_domains.append(cs)
        local=[3*c+j for j in range(3)if 3*c+j in paired]
        for p in local:domains[p]={z[p%3]for z in cs}
        assert set(itertools.product(*(domains[p]for p in local)))=={tuple(z[p%3]for p in local)for z in cs}
    danger=lambda p:'A' if RY[sample[p]] else 'C'
    safe=lambda p:'G' if RY[sample[p]] else 'U'
    fixed_one={p for p,ds in domains.items()if ds=={danger(p)}}
    assert all(not(i in fixed_one and j in fixed_one)for i,j in edges)
    for i,j in edges:
        if i in fixed_one:domains[j].discard(danger(j))
        if j in fixed_one:domains[i].discard(danger(i))
    assert all(domains.values())
    free=sorted(p for p,ds in domains.items()if len(ds)==2)
    baseline={p:next(iter(ds))if len(ds)==1 else safe(p)for p,ds in domains.items()}
    k=1;special=[]
    for c,(a,cs)in enumerate(zip(peptide,codon_domains)):
        local=[p for p in paired if p//3==c]
        count=lambda ass:sum(all(z[p%3]==ass[p]for p in local)for z in cs)
        k0=count(baseline);assert k0>0;k*=k0
        weights={}
        for p in local:
            if p not in free:continue
            ass={**baseline,p:danger(p)};w=Fraction(count(ass),k0);weights[p]=w
            expected=2 if a=='L' and p%3==0 and 3*c+2 not in paired and 3*c+2 not in fixed else 1
            assert w==expected,(peptide,edges,c,p,w,expected,cs)
            if w==2:special.append(p)
        # Extension multiplicities must factor over paired variables in a codon.
        for values in itertools.product(*(domains[p]for p in local)):
            ass={**baseline,**dict(zip(local,values))};pred=Fraction(k0)
            for p in local:
                if ass[p]==danger(p) and p in weights:pred*=weights[p]
            assert pred==count(ass)
    residual=[(free.index(i),free.index(j))for i,j in set(edges)if i in free and j in free]
    z=0;marginals=Counter()
    for mask in range(1<<len(free)):
        if any(mask&(1<<i)and mask&(1<<j)for i,j in residual):continue
        weight=2**sum(bool(mask&(1<<free.index(p)))for p in special);z+=weight
        for i,p in enumerate(free):
            if mask&(1<<i):marginals[p]+=weight
    return {'volume':k*z,'baseline_extension_factor':k,'weighted_independent_set_sum':z,'free_sites0':free,'leucine_activity2_sites0':special,'dangerous_marginal':{str(p):str(Fraction(num,z))for p,num in marginals.items()}}

def all_graph_volumes():
    counts=Counter();special_hist=Counter()
    for m in(1,2):
        all_edges=list(itertools.combinations(range(3*m),2))
        for peptide in itertools.product(CODONS,repeat=m):
            groups={}
            for cs in itertools.product(*(CODONS[a]for a in peptide)):
                x=''.join(cs);allowed=sum(1<<i for i,(a,b)in enumerate(all_edges)if x[a]+x[b]in PAIRS)
                mask=allowed
                while True:
                    if mask not in groups:groups[mask]=defaultdict(list)
                    edges=[e for bit,e in enumerate(all_edges)if mask&(1<<bit)]
                    groups[mask][signature(x,edges,peptide)].append(x)
                    if not mask:break
                    mask=(mask-1)&allowed
            for mask,fibres in groups.items():
                edges=[e for bit,e in enumerate(all_edges)if mask&(1<<bit)]
                for sig,vs in fibres.items():
                    result=volume_certificate(peptide,edges,vs[0])
                    assert result['volume']==len(vs),(peptide,edges,vs,result)
                    counts['nonempty_orientation_serine_fibres']+=1
                    counts['genotypes_sum']+=len(vs)
                    special_hist[len(result['leucine_activity2_sites0'])]+=1
    return {'counts':dict(counts),'activity2_site_count_histogram':dict(special_hist)}

def local_projection_audit():
    counts=Counter()
    # Every fixed-base pattern and every subset/class assignment of paired sites
    # for every amino acid, plus both Ser families. No graph assumption needed.
    for a,allcs in CODONS.items():
        for requirements in itertools.product((None,'A','C','G','U'),repeat=3):
            for classes in itertools.product((None,0,1),repeat=3):
                for family in((False,True)if a=='S'else(None,)):
                    cs=[z for z in allcs if all((requirements[j]is None or z[j]==requirements[j])and
                        (classes[j]is None or RY[z[j]]==classes[j])for j in range(3))and
                        (family is None or(z[:2]=='AG')==family)]
                    counts['local_cases']+=1
                    if not cs:continue
                    local=[j for j,c in enumerate(classes)if c is not None]
                    ds=[{z[j]for z in cs}for j in local]
                    assert set(itertools.product(*ds))=={tuple(z[j]for j in local)for z in cs}
                    counts['nonempty_local_cases']+=1
    return dict(counts)

def star_volumes():
    output=[]
    for r in range(1,11):
        peptide='L'+'K'*r;sample='CUC'+'AAG'*r;edges=[(0,3*i+2)for i in range(1,r+1)]
        cert=volume_certificate(peptide,edges,sample)
        assert cert['volume']==4+2**(r+1)
        assert Fraction(cert['dangerous_marginal']['0'])==Fraction(1,1+2**(r-1))
        cert['partners']=r;output.append(cert)
    return output

def tree_DP_validation():
    import random
    rng=random.Random(2026102402);checked=0
    for n in range(1,13):
        for _ in range(8):
            relations=set()
            for v in range(1,n):
                u=rng.randrange(v)
                relations.add((u,v)if rng.randrange(2)else(v,u))
            pre=[0]*n
            for u,v in relations:pre[v]|=1<<u
            dp={0:1}
            for mask in range(1<<n):
                if mask not in dp:continue
                for v in range(n):
                    if not mask&(1<<v)and mask&pre[v]==pre[v]:dp[mask|1<<v]=dp.get(mask|1<<v,0)+dp[mask]
            assert forest_path_count(n,relations)==dp[(1<<n)-1]
            checked+=1
    large=[]
    for r in(10,30,60):
        # path event0->event1->r interchangeable leaves.
        relations={(0,1)}|{(1,v)for v in range(2,r+2)}
        value=forest_path_count(r+2,relations);assert value==math.factorial(r)
        large.append({'events':r+2,'path_count':str(value)})
    return {'random_oriented_tree_fixtures':checked,'seed':2026102402,'large_exact_star_cases':large}

def main():
    if not __debug__:raise RuntimeError('assertions must remain enabled')
    start=time.perf_counter();import argparse;parser=argparse.ArgumentParser();parser.add_argument('--outdir');args=parser.parse_args();root=Path(args.outdir)if args.outdir else Path(__file__).resolve().parents[3]/'results'/'expansion'/'verification';root.mkdir(parents=True,exist_ok=True)
    result={'local_projection':local_projection_audit(),'global_volume':all_graph_volumes(),'leucine_star_family':star_volumes(),'forest_DP':tree_DP_validation(),'seconds':time.perf_counter()-start}
    (root/'hardcore_volume_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
