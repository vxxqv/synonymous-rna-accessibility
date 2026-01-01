"""Exact standard-code accessibility for arbitrary simultaneous pairing graphs.

General pairing graphs describe compatibility with every required target pair;
they are not a physical model allowing simultaneous multiple partners.
"""
from collections import Counter,deque
from fractions import Fraction
import itertools
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'exact'))
from independent_oracle import CODE,CODONS,PAIRS,RY

def normalize(peptide,edges):
    if any(a not in CODONS for a in peptide):raise ValueError('standard sense peptide required')
    n=3*len(peptide);result=set()
    for edge in edges:
        if len(edge)!=2:raise ValueError('pair must have two positions')
        i,j=edge
        if not isinstance(i,int)or not isinstance(j,int)or i==j or not(0<=i<n and 0<=j<n):raise ValueError('invalid pairing edge')
        result.add(tuple(sorted((i,j))))
    return tuple(sorted(result))

def valid(peptide,edges,x,fixed=None):
    return len(x)==3*len(peptide)and all(x[3*c:3*c+3]in CODONS[a]for c,a in enumerate(peptide))and all(x[i]+x[j]in PAIRS for i,j in edges)and all(x[p]==b for p,b in(fixed or{}).items())

def signature(x,edges,peptide):
    paired=sorted({p for e in edges for p in e})
    return tuple(RY[x[p]]for p in paired)+tuple(x[3*c:3*c+2]=='AG'for c,a in enumerate(peptide)if a=='S')

def dependency_poset(x,y,peptide,edges):
    """Compile local endpoint subcube factors; no global genotype search."""
    positions=[p for p in range(len(x))if x[p]!=y[p]];index={p:i for i,p in enumerate(positions)};relations=set()
    factors=[(list(range(3*c,3*c+3)),lambda z,a=a:CODE[z]==a)for c,a in enumerate(peptide)]
    factors += [([i,j],lambda z:z in PAIRS)for i,j in edges]
    for scope,test in factors:
        changed=[p for p in scope if p in index];good=[]
        for mask in range(1<<len(changed)):
            z=''.join(y[p]if p in changed and mask&(1<<changed.index(p))else x[p]for p in scope)
            if test(z):good.append(mask)
        local={(a,b)for a,b in itertools.permutations(range(len(changed)),2)if all(not(m&(1<<b))or m&(1<<a)for m in good)}
        assert set(good)=={m for m in range(1<<len(changed))if all(not(m&(1<<b))or m&(1<<a)for a,b in local)}
        relations.update((index[changed[a]],index[changed[b]])for a,b in local)
    k=len(positions);indeg=[0]*k;out=[[]for _ in range(k)];height=[1]*k
    for a,b in relations:indeg[b]+=1;out[a].append(b)
    queue=deque(i for i in range(k)if not indeg[i]);order=[]
    while queue:
        a=queue.popleft();order.append(a)
        for b in out[a]:
            height[b]=max(height[b],height[a]+1);indeg[b]-=1
            if not indeg[b]:queue.append(b)
    assert len(order)==k and max(height,default=0)<=3
    return positions,relations,height,order

def _multiply(a,b):
    out=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def _integrate(poly,upper):
    ant=[Fraction(0)]+[c/Fraction(i+1)for i,c in enumerate(poly)]
    return [sum(ant)]+[-c for c in ant[1:]]if upper else ant

def forest_extensions(k,relations):
    """Return e(P) for a forest; None if the undirected graph has a cycle."""
    adj=[set()for _ in range(k)]
    for a,b in relations:adj[a].add(b);adj[b].add(a)
    seen=set();components=[]
    for root in range(k):
        if root in seen:continue
        component={root};queue=[root];seen.add(root)
        for v in queue:
            for w in adj[v]-seen:seen.add(w);component.add(w);queue.append(w)
        if sum(len(adj[v])for v in component)//2!=len(component)-1:return None
        components.append((root,component))
    def visit(v,parent):
        f=[Fraction(1)]
        for child in adj[v]-{parent}:
            f=_multiply(f,_integrate(visit(child,v),upper=(v,child)in relations))
        return f
    volume=Fraction(1)
    for root,component in components:
        poly=visit(root,None);volume*=sum(c/Fraction(i+1)for i,c in enumerate(poly))
    value=math.factorial(k)*volume;assert value.denominator==1
    return value.numerator

def subset_extensions(k,relations,cap=24):
    if k>cap:return None
    pre=[0]*k
    for a,b in relations:pre[b]|=1<<a
    dp={0:1}
    for mask in range(1<<k):
        if mask not in dp:continue
        for b in range(k):
            if not mask&(1<<b)and mask&pre[b]==pre[b]:dp[mask|1<<b]=dp.get(mask|1<<b,0)+dp[mask]
    return dp.get((1<<k)-1,0)

def endpoint_certificate(peptide,edges,x,y,fixed=None,subset_cap=24):
    """Connectivity/distance certificate plus exact path count when permitted.

    Path-count failure is reported separately from connectivity. Mutation event
    positions and pair edges are0-based; predecessor edges index event vertices.
    """
    edges=normalize(peptide,edges)
    if not valid(peptide,edges,x,fixed)or not valid(peptide,edges,y,fixed):raise ValueError('incompatible endpoint')
    if signature(x,edges,peptide)!=signature(y,edges,peptide):return {'status':'disconnected','distance':None,'shortest_paths':'0'}
    positions,relations,height,order=dependency_poset(x,y,peptide,edges);k=len(positions)
    num=forest_extensions(k,relations);backend='forest_order_polytope'
    if num is None:num=subset_extensions(k,relations,subset_cap);backend='subset_ideals'
    witness_order=[positions[a]for a in order];z=list(x)
    for p in witness_order:z[p]=y[p];assert valid(peptide,edges,''.join(z),fixed)
    return {'status':'connected','distance':k,'changed_positions0':positions,'precedence_event_edges':sorted(relations),'dependency_height':max(height,default=0),'witness_order_positions0':witness_order,'shortest_paths':None if num is None else str(num),'count_status':'resource_limit'if num is None else'exact','count_backend':backend}

def weighted_projection(peptide,edges,sample,fixed=None,max_free=24):
    """Compile uniform component volume as K times a hard-core partition sum.

    Every residual site has activity1 except the stated Leu-first sites with
    activity2. Exact subset evaluation is capped; the projection is returned
    regardless, enabling separate width-bounded/counting backends.
    """
    edges=normalize(peptide,edges);fixed=fixed or{}
    if not valid(peptide,edges,sample,fixed):raise ValueError('incompatible sample')
    paired={p for e in edges for p in e};domains={};codon_domains=[]
    for c,a in enumerate(peptide):
        cs=[z for z in CODONS[a]if all((3*c+j not in paired or RY[z[j]]==RY[sample[3*c+j]])and(3*c+j not in fixed or z[j]==fixed[3*c+j])for j in range(3))and(a!='S'or(z[:2]=='AG')==(sample[3*c:3*c+2]=='AG'))]
        codon_domains.append(cs)
        for p in paired:
            if p//3==c:domains[p]={z[p%3]for z in cs}
    danger=lambda p:'A'if RY[sample[p]]else'C'
    safe=lambda p:'G'if RY[sample[p]]else'U'
    forced_one={p for p,ds in domains.items()if ds=={danger(p)}}
    for i,j in edges:
        if i in forced_one:domains[j].discard(danger(j))
        if j in forced_one:domains[i].discard(danger(i))
    assert all(domains.values())
    free=sorted(p for p,ds in domains.items()if len(ds)==2)
    baseline={p:next(iter(ds))if len(ds)==1 else safe(p)for p,ds in domains.items()}
    baseline_factor=1;special=[]
    for c,(a,cs)in enumerate(zip(peptide,codon_domains)):
        local=[p for p in paired if p//3==c]
        count=lambda ass:sum(all(z[p%3]==ass[p]for p in local)for z in cs)
        k0=count(baseline);assert k0>0;baseline_factor*=k0
        weights={}
        for p in local:
            if p not in free:continue
            value=Fraction(count({**baseline,p:danger(p)}),k0)
            expected=2 if a=='L'and p%3==0 and 3*c+2 not in paired and 3*c+2 not in fixed else 1
            assert value==expected;weights[p]=value
            if value==2:special.append(p)
        for values in itertools.product(*(domains[p]for p in local)):
            ass={**baseline,**dict(zip(local,values))};pred=Fraction(k0)
            for p in local:
                if ass[p]==danger(p)and p in weights:pred*=weights[p]
            assert pred==count(ass)
    residual=[(free.index(i),free.index(j))for i,j in edges if i in free and j in free]
    result={'baseline_extension_factor':str(baseline_factor),'free_sites0':free,'leucine_activity2_sites0':special,'residual_edges_on_free_indices':residual,'baseline_paired_bases':baseline}
    if len(free)>max_free:return {**result,'status':'resource_limit','volume':None}
    z=0;marginals=Counter()
    for mask in range(1<<len(free)):
        if any(mask&(1<<i)and mask&(1<<j)for i,j in residual):continue
        weight=2**sum(bool(mask&(1<<free.index(p)))for p in special);z+=weight
        for i,p in enumerate(free):
            if mask&(1<<i):marginals[p]+=weight
    return {**result,'status':'exact','volume':str(baseline_factor*z),'weighted_independent_set_sum':str(z),'dangerous_marginals':{str(p):str(Fraction(v,z))for p,v in marginals.items()}}
