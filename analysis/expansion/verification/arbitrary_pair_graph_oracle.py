"""Exact counterexample search/generalization pilot on simultaneous pair graphs.

Independent finite universe: standard sense peptides with <=2 codons and EVERY
simple pairing graph on their nucleotide sites. Nonempty instances are covered
by enumerating every genotype and every subset of that genotype's allowed-pair
edges. This proven reduction includes graphs of arbitrary degree/odd cycles
(the latter are empty and counted by subtraction).
"""
from collections import Counter,defaultdict,deque
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'exact'))
from independent_oracle import CODE,CODONS,BASES,PAIRS,RY,graph,bfs,distance

def signature(x,edges,peptide):
    paired=sorted({p for e in edges for p in e})
    return tuple(RY[x[p]] for p in paired)+tuple(x[3*i:3*i+2]=='AG' for i,a in enumerate(peptide) if a=='S')

def dependencies(x,y,peptide,edges):
    positions=[p for p in range(len(x)) if x[p]!=y[p]]
    index={p:i for i,p in enumerate(positions)};relations=set()
    factors=[(list(range(3*c,3*c+3)), lambda z,a=a: CODE[z]==a) for c,a in enumerate(peptide)]
    factors += [([i,j],lambda z:z in PAIRS) for i,j in edges]
    for scope,valid in factors:
        changed=[p for p in scope if p in index];good=[]
        for mask in range(1<<len(changed)):
            z=''.join(y[p] if p in changed and mask&(1<<changed.index(p)) else x[p] for p in scope)
            if valid(z):good.append(mask)
        local=set()
        for a,b in itertools.permutations(range(len(changed)),2):
            if all(not(m&(1<<b)) or m&(1<<a) for m in good):
                local.add((a,b));relations.add((index[changed[a]],index[changed[b]]))
        ideals=[m for m in range(1<<len(changed)) if all(not(m&(1<<b)) or m&(1<<a) for a,b in local)]
        assert set(ideals)==set(good),(x,y,scope,'non-poset local relation')
    k=len(positions);pre=[0]*k
    for a,b in relations:pre[b]|=1<<a
    dp={0:1};height=[0]*k
    for mask in range(1<<k):
        if mask not in dp:continue
        for b in range(k):
            if mask&(1<<b) or pre[b]&mask !=pre[b]:continue
            dp[mask|(1<<b)]=dp.get(mask|(1<<b),0)+dp[mask]
    assert (1<<k)-1 in dp,(x,y,relations,'directed cycle')
    for _ in range(k):
        for b in range(k):height[b]=max([1]+[height[a]+1 for a,bb in relations if bb==b])
    assert max(height,default=0)<=3,(x,y,relations,height,'height>3')
    return positions,relations,dp[(1<<k)-1],height

def integrate(poly,upper=True):
    """Return integral_0^t f (upper=False), or integral_t^1 f (upper=True)."""
    ant=[Fraction(0)]+[c/Fraction(i+1) for i,c in enumerate(poly)]
    if not upper:return ant
    result=[sum(ant)]+[-c for c in ant[1:]]
    return result

def multiply(a,b):
    c=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def forest_path_count(k,relations):
    """Exact order-polytope polynomial DP if undirected dependencies are a forest."""
    adjacent=[set() for _ in range(k)]
    for a,b in relations:adjacent[a].add(b);adjacent[b].add(a)
    unseen=set(range(k));volume=Fraction(1)
    def recur(v,parent,visited):
        assert v not in visited,'not a forest';visited.add(v)
        poly=[Fraction(1)]
        for c in adjacent[v]-{parent}:
            child=recur(c,v,visited)
            # child before v: integrate0..t; v before child: integratet..1.
            branch=integrate(child,upper=(v,c) in relations)
            poly=multiply(poly,branch)
        return poly
    while unseen:
        r=min(unseen);visited=set();poly=recur(r,None,visited);unseen-=visited
        volume*=sum(c/Fraction(i+1) for i,c in enumerate(poly))
    answer=math.factorial(k)*volume
    assert answer.denominator==1
    return answer.numerator

def audit(max_codons=2):
    counts=Counter();max_degree=Counter();height_hist=Counter();distance_hist=Counter();ratio_hist=Counter();joint_hist=Counter();peptide_summaries=[]
    for m in range(1,max_codons+1):
        sites=3*m;all_edges=list(itertools.combinations(range(sites),2))
        counts['all_instances_including_empty']+=20**m*2**len(all_edges)
        for peptide in itertools.product(CODONS,repeat=m):
            ps={'peptide':''.join(peptide),'codons':m,'nonempty_pair_graphs':0,'multi_component_pair_graphs':0,'maximum_components':0,'maximum_vertices':0};max_effective=Fraction(0)
            instances=defaultdict(list)
            for cs in itertools.product(*(CODONS[a] for a in peptide)):
                x=''.join(cs)
                admissible=sum(1<<e for e,(i,j) in enumerate(all_edges) if x[i]+x[j] in PAIRS)
                mask=admissible
                while True:
                    instances[mask].append(x)
                    if mask==0:break
                    mask=(mask-1)&admissible
            for mask,vs in instances.items():
                edges=[e for bit,e in enumerate(all_edges) if mask&(1<<bit)]
                counts['nonempty_instances']+=1;counts['vertices_sum']+=len(vs)
                degrees=Counter(p for e in edges for p in e);max_degree[max(degrees.values(),default=0)]+=1
                adj=graph(vs);counts['edges_sum']+=sum(map(len,adj))//2
                sigs=[signature(x,edges,peptide) for x in vs]
                sizes=Counter(sigs);effective=Fraction(len(vs)**2,sum(v*v for v in sizes.values()))
                ps['nonempty_pair_graphs']+=1;ps['multi_component_pair_graphs']+=len(sizes)>1
                ps['maximum_components']=max(ps['maximum_components'],len(sizes));ps['maximum_vertices']=max(ps['maximum_vertices'],len(vs));max_effective=max(max_effective,effective)
                for a,x in enumerate(vs):
                    ds,npaths=bfs(adj,a)
                    for b in range(a+1,len(vs)):
                        counts['endpoint_pairs']+=1;same=sigs[a]==sigs[b]
                        assert same==(b in ds),(peptide,edges,x,vs[b],'component')
                        if not same:
                            counts['separated_pairs']+=1;distance_hist[('separated',distance(x,vs[b]))]+=1;continue
                        counts['connected_pairs']+=1
                        y=vs[b];assert ds[b]==distance(x,y)
                        positions,relations,num,height=dependencies(x,y,peptide,edges)
                        assert npaths[b]==num,(peptide,edges,x,y,npaths[b],num)
                        k=len(positions);fraction=Fraction(num,math.factorial(k));distance_hist[('connected',k)]+=1;ratio_hist[str(fraction)]+=1;joint_hist[(k,num,math.factorial(k),max(height,default=0))]+=1
                        height_hist[max(height,default=0)]+=1
                        # Validate every forest case with the polynomial backend;
                        # cyclic cases retain the independent subset-ideal count.
                        try:fnum=forest_path_count(len(positions),relations)
                        except AssertionError as err:
                            if str(err)=='not a forest':continue
                            raise
                        assert fnum==num
                        counts['forest_DP_pairs']+=1
            ps['maximum_effective_components_fraction']=str(max_effective);peptide_summaries.append(ps)
    counts['empty_instances']=counts['all_instances_including_empty']-counts['nonempty_instances']
    return {'counts':dict(counts),'max_pair_degree_counts':dict(max_degree),'dependency_height_counts':dict(height_hist),'distance_counts':[{'status':s,'hamming_distance':d,'endpoint_pairs':v}for(s,d),v in sorted(distance_hist.items())],'mutation_order_fraction_counts':dict(ratio_hist),'distance_order_height_counts':[{'hamming_distance':k,'valid_shortest_orders':v,'all_mutation_orders':f,'dependency_height':h,'endpoint_pairs':n}for(k,v,f,h),n in sorted(joint_hist.items())],'peptide_component_summaries':peptide_summaries}

def star_fixture(r=2):
    peptide='L'+'K'*r;x='CUC'+'AAG'*r;y='UUA'+'AAA'*r
    edges=[(0,3*i+2)for i in range(1,r+1)]
    positions,relations,num,height=dependencies(x,y,peptide,edges)
    assert num==math.factorial(r)
    return {'peptide':peptide,'x':x,'y':y,'edges':edges,'changed_coordinates0':positions,'precedences_on_event_indices':sorted(relations),'shortest_paths':num,'height':max(height),'forest_DP':forest_path_count(len(positions),relations)}

def height2_fixture(a,b,edge_mask):
    peptide='K'*a+'F'*b
    x='AAA'*a+'UUU'*b;y='AAG'*a+'UUC'*b
    edges=[(3*i+2,3*(a+j)+2)for i in range(a)for j in range(b)if edge_mask&(1<<(i*b+j))]
    positions,relations,num,height=dependencies(x,y,peptide,edges)
    target={(i,a+j)for i in range(a)for j in range(b)if edge_mask&(1<<(i*b+j))}
    assert relations==target
    return num

def main():
    if not __debug__:raise RuntimeError('assertions must remain enabled')
    start=time.perf_counter();import argparse;parser=argparse.ArgumentParser();parser.add_argument('--outdir');args=parser.parse_args();root=Path(args.outdir)if args.outdir else Path(__file__).resolve().parents[3]/'results'/'expansion'/'verification';root.mkdir(parents=True,exist_ok=True)
    results=audit();results['star_fixtures']=[star_fixture(r)for r in range(2,8)]
    # Every labelled2x3 bipartite precedence graph is realized parsimoniously.
    for mask in range(1<<6):height2_fixture(2,3,mask)
    results['height2_reduction_graphs_checked']=64
    results['seconds']=time.perf_counter()-start
    (root/'arbitrary_pair_graph_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
