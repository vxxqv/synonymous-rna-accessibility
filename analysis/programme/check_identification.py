"""Finite independent falsification checks for endpoint identification.

The reference oracle uses all simple paths, never the heap implementation.
All integer computations are exact; finite checks supplement the analytic proof.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import sys,json,hashlib,argparse,time
if not __debug__:
    raise RuntimeError('Exact verification requires enabled assertions')
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from barrier_identification import barrier,sharp_barrier_bounds,threshold_status,additive_background_barrier


def paths(adj,a,b):
    if a==b:return [(a,)]
    out=[]
    def visit(p):
        if p[-1]==b:out.append(tuple(p));return
        for v in adj[p[-1]]:
            if v not in p:visit(p+[v])
    visit([a]);return out


def reference_B(ps,f,a,b):
    return min(f[a],f[b])-max(min(f[v] for v in p) for p in ps)


def graphs(n):
    edges=list(combinations(range(n),2))
    for bits in product((0,1),repeat=len(edges)):
        adj=[[] for _ in range(n)]
        for keep,(u,v) in zip(bits,edges):
            if keep:adj[u].append(v);adj[v].append(u)
        if all(paths(adj,0,v) for v in range(n)):yield adj


def product_graph(g,h):
    n=len(h);adj=[[] for _ in range(len(g)*n)]
    for x,z in product(range(len(g)),range(n)):
        adj[x*n+z]=[y*n+z for y in g[x]]+[x*n+w for w in h[z]]
    return adj


def threshold_B(adj,f,a,b):
    # Independent connectivity oracle for product graphs with many simple paths.
    for t in sorted(set(f),reverse=True):
        if min(f[a],f[b])<t:continue
        seen={a};stack=[a]
        while stack:
            for v in adj[stack.pop()]:
                if v not in seen and f[v]>=t:seen.add(v);stack.append(v)
        if b in seen:return min(f[a],f[b])-t
    raise AssertionError('disconnected')


def check_boxes():
    counts=Counter();intervals=[(l,u) for l in (-1,0,1) for u in (-1,0,1) if l<=u]
    for n in range(1,5):
        for adj in graphs(n):
            counts['connected_labelled_graphs']+=1
            for a in range(n):assert barrier(adj,[0]*n,a,a)==0
            for a,b in combinations(range(n),2):
                ps=paths(adj,a,b);all_fields=list(product((-1,0,1),repeat=n))
                oracle={f:reference_B(ps,f,a,b) for f in all_fields}
                for f,v in oracle.items():
                    assert barrier(adj,f,a,b)==v;counts['widest_path_oracle_comparisons']+=1
                    cubed=tuple(x**3+2*x for x in f)
                    assert (barrier(adj,cubed,a,b)>0)==(v>0)
                    clipped=tuple(max(0,x) for x in f)
                    assert v>0 or barrier(adj,clipped,a,b)==0
                for box in product(intervals,repeat=n):
                    l,u=zip(*box);lo,hi=sharp_barrier_bounds(adj,l,u,a,b)
                    corners={tuple(x) for x in product(*(set(p) for p in box))}
                    vals=[oracle[x] for x in corners]
                    assert (lo,hi)==(min(vals),max(vals));counts['sharp_box_cases']+=1
                    for t in (0,1):
                        if min(l[a],l[b])<t:continue
                        status=threshold_status(adj,l,u,a,b,t)
                        decisions=[max(min(f[v] for v in p) for p in ps)>=t for f in corners]
                        expected='connected' if all(decisions) else 'separated' if not any(decisions) else 'unresolved'
                        assert status==expected;counts['threshold_completion_cases']+=1
        print('box n',n,dict(counts),flush=True)
    return dict(counts)


def check_backgrounds():
    counts=Counter();examples=[]
    for n in (2,3):
        for g in graphs(n):
            for d in (1,2):
                h=[[j for j in range(2**d) if (i^j).bit_count()==1] for i in range(2**d)]
                gh=product_graph(g,h)
                for weights in product((-1,0,1),repeat=d):
                    hv=[sum(w for k,w in enumerate(weights) if z&(1<<k)) for z in range(2**d)]
                    for gv in product((-1,0,1),repeat=n):
                        f=[v+w for v in gv for w in hv]
                        for a,b in combinations(range(n),2):
                            local=reference_B(paths(g,a,b),gv,a,b)
                            for z in range(2**d):
                                expected=additive_background_barrier(local,hv[z],max(hv))
                                actual=threshold_B(gh,f,a*len(h)+z,b*len(h)+z)
                                assert expected==actual;counts['product_cases']+=1
                                if local>0 and actual==0:counts['erased_positive_local_barriers']+=1
    # No-superlevel-access counterexample: background peak lies beyond a valley.
    g=[[1],[0,2],[1]];h=[[1],[0,2],[1]];gv=[2,0,2];hv=[0,-10,3]
    gh=product_graph(g,h);f=[x+y for x in gv for y in hv]
    actual=threshold_B(gh,f,0,6);unjustified=max(0,2-(3-0))
    assert actual==2 and unjustified==0
    examples.append({'hypothesis_removed':'background maximum reachable without dropping below initial level','g':gv,'h':hv,'actual_barrier':actual,'invalid_formula':unjustified})
    # Sharpness of the universal 2d perturbation bound.
    f=[1,0,1];perturbed=[2,-1,2]
    assert barrier(g,perturbed,0,2)-barrier(g,f,0,2)==2
    return {'counts':dict(counts),'necessary_hypothesis_counterexample':examples,'perturbation_constant_2_is_attained':True}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();t=time.perf_counter()
    result={'status':'completed','scope':'Finite exact checks, not universal proof','boxes':check_boxes(),'background':check_backgrounds(),'seconds':time.perf_counter()-t,'input_sha256':{x.relative_to(ROOT).as_posix():hashlib.sha256(x.read_bytes()).hexdigest() for x in [Path(__file__),ROOT/'src/barrier_identification.py']}}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
