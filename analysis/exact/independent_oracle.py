"""Independent finite audit: direct translation, explicit graphs, binary-subcube factors.

No code is imported from the primary implementation. RNA standard genetic code.
Run with Python >=3.10. All counts are Python arbitrary-precision integers.
"""
import collections
import itertools
import json
import math
import pathlib
import time

BASES = 'UCAG'
AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = dict(zip(map(''.join, itertools.product(BASES, repeat=3)), AAS))
CODONS = {a: tuple(c for c in CODE if CODE[c] == a) for a in sorted(set(AAS)-{'*'})}
PAIRS = set(('AU','UA','GC','CG','GU','UG'))
RY = {b: int(b in 'AG') for b in BASES}

def matchings(vertices):
    if not vertices:
        yield ()
        return
    a, *rest = vertices
    yield from matchings(rest)
    for b in rest:
        for tail in matchings([c for c in rest if c != b]):
            yield ((a,b),)+tail

def vertices(peptide, pairs):
    return [''.join(cs) for cs in itertools.product(*(CODONS[a] for a in peptide))
            if all((''.join(cs)[i]+''.join(cs)[j]) in PAIRS for i,j in pairs)]

def signature(x, pairs, peptide):
    return tuple(RY[x[i]] for i,j in pairs)+tuple(x[3*i:3*i+2]=='AG' for i,a in enumerate(peptide) if a=='S')

def distance(x,y):
    return sum(a!=b for a,b in zip(x,y))

def graph(vs):
    lookup={v:i for i,v in enumerate(vs)}
    adj=[[] for _ in vs]
    for i,x in enumerate(vs):
        for p in range(len(x)):
            for b in BASES:
                y=x[:p]+b+x[p+1:]
                if y in lookup and y!=x:
                    adj[i].append(lookup[y])
    return adj

def bfs(adj, source):
    ds={source:0}; count={source:1}; q=collections.deque([source])
    while q:
        a=q.popleft()
        for b in adj[a]:
            if b not in ds:
                ds[b]=ds[a]+1; count[b]=0; q.append(b)
            if ds[b]==ds[a]+1:
                count[b]+=count[a]
    return ds,count

def certificate(x,y,peptide,pairs):
    positions=[p for p in range(len(x)) if x[p]!=y[p]]
    index={p:i for i,p in enumerate(positions)}
    edges=set()
    factors=[(list(range(3*c,3*c+3)), lambda z,a=a: CODE[z]==a) for c,a in enumerate(peptide)]
    factors += [([i,j],lambda z: z in PAIRS) for i,j in pairs]
    for scope,valid in factors:
        changed=[p for p in scope if p in index]
        good=[]
        for mask in range(1<<len(changed)):
            z=''.join(y[p] if p in changed and mask&(1<<changed.index(p)) else x[p] for p in scope)
            if valid(z): good.append(mask)
        local=set()
        for a,b in itertools.permutations(range(len(changed)),2):
            if all(not (m&(1<<b)) or (m&(1<<a)) for m in good):
                local.add((a,b)); edges.add((index[changed[a]],index[changed[b]]))
        ideals=[m for m in range(1<<len(changed)) if all(not(m&(1<<b)) or m&(1<<a) for a,b in local)]
        assert set(ideals)==set(good), (x,y,scope,'factor not order ideals',good,ideals)
    # Remove transitive edges, then certify union of chains with <=3 vertices.
    reduction=set(edges)
    for a,b in edges:
        if any((a,c) in edges and (c,b) in edges for c in range(len(positions))):
            reduction.remove((a,b))
    indeg=collections.Counter(b for a,b in reduction); outdeg=collections.Counter(a for a,b in reduction)
    assert all(v<=1 for v in (*indeg.values(),*outdeg.values())), (x,y,reduction,'branch')
    seen=set(); sizes=[]
    for a in range(len(positions)):
        if indeg[a]: continue
        size=0
        while a not in seen:
            seen.add(a); size+=1
            nxt=[b for aa,b in reduction if aa==a]
            if not nxt: break
            a=nxt[0]
        sizes.append(size)
    assert len(seen)==len(positions) and max(sizes,default=0)<=3,(x,y,reduction,sizes,'long chain or cycle')
    return math.factorial(len(positions))//math.prod(math.factorial(s) for s in sizes),sizes

def audit_two_codons():
    counts=collections.Counter(); chain_shapes=collections.Counter()
    all_matchings=list(matchings(list(range(6))))
    assert len(all_matchings)==76
    for peptide in itertools.product(CODONS,repeat=2):
        for pairs in all_matchings:
            counts['instances']+=1
            vs=vertices(peptide,pairs); counts['vertices_sum']+=len(vs)
            if not vs:
                counts['empty_instances']+=1; continue
            counts['nonempty_instances']+=1
            adj=graph(vs)
            counts['edges_sum']+=sum(map(len,adj))//2
            for a,x in enumerate(vs):
                ds,npaths=bfs(adj,a)
                for b in range(a+1,len(vs)):
                    y=vs[b]; counts['endpoint_pairs']+=1
                    same=signature(x,pairs,peptide)==signature(y,pairs,peptide)
                    assert same==(b in ds),(peptide,pairs,x,y,'component claim')
                    if not same: counts['separated_pairs']+=1; continue
                    counts['connected_pairs']+=1
                    assert ds[b]==distance(x,y),(peptide,pairs,x,y,'geodesic claim')
                    predicted,sizes=certificate(x,y,peptide,pairs)
                    assert predicted==npaths[b],(peptide,pairs,x,y,predicted,npaths[b],'path count')
                    chain_shapes[','.join(map(str,sorted(sizes)))]+=1
    return {'counts':dict(counts),'chain_shape_counts':dict(chain_shapes)}

def register(n,d):
    return tuple((i,n+(n-1+d-i)) for i in range(d,n))

def bridge_audit():
    result={}
    rows=[]
    for d in [1,2,3]:
        counts=collections.Counter()
        for peptide in itertools.product(CODONS,repeat=2):
            a=vertices(peptide,register(3,0)); b=vertices(peptide,register(3,d))
            if not a or not b: continue
            delta=min(distance(x,y) for x in a for y in b)
            counts[delta]+=1
            inter=sorted(set(a)&set(b))
            adj=graph(inter); unseen=set(range(len(inter))); shapes=[]
            while unseen:
                j=min(unseen); ds,_=bfs(adj,j); nodes=set(ds); unseen-=nodes
                deg=[len(adj[u]) for u in nodes]
                shapes.append((len(nodes),sorted(set(deg))))
            rows.append({'peptide':''.join(peptide),'shift':d,'layer0':len(a),'layerd':len(b),'intersection':len(inter),'delta':delta,'component_sizes_and_degrees':shapes})
        result[str(d)]=dict(sorted(counts.items()))
    return result,rows

def main():
    t=time.perf_counter()
    root=pathlib.Path(__file__).resolve().parent
    results=audit_two_codons()
    results['register_delta_counts'],rows=bridge_audit()
    results['seconds']=time.perf_counter()-t
    (root/'independent_results.json').write_text(json.dumps(results,indent=2)+'\n')
    (root/'register_examples.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__': main()
