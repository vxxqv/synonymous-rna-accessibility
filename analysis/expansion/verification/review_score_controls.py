"""Independent small exact audit of conditional nulls and filtration conventions."""
from collections import Counter
from fractions import Fraction
import itertools
import json
from pathlib import Path
import random
import sys
EXP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(EXP))
from conditional_vertex_null import conditional_null,degree_bin
from retention_filtration import filtration

def comps(graph,keep=None):
    unseen=set(range(len(graph)))if keep is None else set(keep);result=[]
    while unseen:
        root=min(unseen);seen={root};stack=[root];unseen.remove(root)
        while stack:
            for j in graph[stack.pop()]:
                if j in unseen:unseen.remove(j);seen.add(j);stack.append(j)
        result.append(seen)
    return result

def labels_for(graph):return {v:c for c,part in enumerate(comps(graph))for v in part}

def bfs_curve_and_bars(graph,values):
    labels=labels_for(graph);curve=[];old_peaks=set();bars=[]
    for level in sorted(set(values),reverse=True):
        parts=comps(graph,{i for i,v in enumerate(values)if v>=level});counts=Counter(labels[i]for p in parts for i in p)
        excess=1-Fraction(sum(len(p)**2 for p in parts),sum(v*v for v in counts.values()))
        curve.append((level,excess));new_peaks=set()
        for p in parts:
            elder=max(p,key=lambda i:(values[i],-i));new_peaks.add(elder)
            for dead in old_peaks&p-{elder}:
                if values[dead]>level:bars.append((values[dead],level))
        old_peaks=new_peaks
    area=sum(e*(level-(curve[i+1][0]if i+1<len(curve)else 0))for i,(level,e)in enumerate(curve))
    return curve,sorted(bars),area

def audit_filtration():
    rng=random.Random(2026102404);checked=0
    fixtures=[([[1],[0,2],[1]],[1,0,1]),([[1],[0,2],[1]],[1,.5,1]),([[1],[0,2],[1]],[1,1,1]),([[],[]],[1,1])]
    for _ in range(100):
        n=rng.randrange(1,9);g=[[]for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if rng.randrange(3)==0:g[i].append(j);g[j].append(i)
        fixtures.append((g,[rng.randrange(5)/4 for _ in range(n)]))
    for graph,values in fixtures:
        curve,bars,area=bfs_curve_and_bars(graph,values);observed=filtration(graph,values,full=True)
        assert len(observed['curve'])==len(curve)
        assert all(r['threshold']==float(v)and abs(r['excess']-float(e))<1e-12 for r,(v,e)in zip(observed['curve'],curve))
        assert sorted((b['birth'],b['death'])for b in observed['bars'])==bars
        assert abs(observed['integrated_excess']-float(area))<1e-12
        assert (observed['positive_bar_count']>0)==(max((e for v,e in curve),default=0)>0)
        checked+=1
    return checked

def audit_conditional():
    g=[[1],[0,2],[1]];lab=labels_for(g)
    a=conditional_null(g,lab,{0,2},False,42)
    assert a['exact']and a['population_size']=='3'
    assert abs(a['null_mean']-1/6)<1e-14 and a['observed_excess']==.5 and a['lower_tail_p']==1
    b=conditional_null(g,lab,{0,1},False,42)
    assert b['lower_tail_p']==2/3
    c=conditional_null(g,lab,{0,2},True,42)
    assert c['population_size']=='1'and c['null_mean']==.5 and not c['statistic_varied']
    h=[[1],[0,2],[1],[]];l=labels_for(h)
    d=conditional_null(h,l,{0,2,3},False,42)
    assert d['connection_denominator']=='5'and abs(d['null_mean']-2/15)<1e-14
    # Sampling branch reproducibility, not a new large Monte Carlo experiment.
    x=conditional_null(g,lab,{0,2},False,42,draws=10,cap=0)
    y=conditional_null(g,lab,{0,2},False,42,draws=10,cap=0)
    assert x==y and not x['exact']and x['draws']==10
    values=[Fraction(0),Fraction(1,2),Fraction(1)];edges=[(0,1),(1,2)]
    direct=sum(sum((p[i]-p[j])**2 for i,j in edges)for p in itertools.permutations(values))/6
    var=sum((v-sum(values)/3)**2 for v in values)/3
    assert direct==len(edges)*Fraction(2*3,3-1)*var
    return {'exact_hand_fixtures':4,'sampling_reproducibility_draws':10,'roughness_expectation_checked':str(direct)}

def independent_set_graph(n,edges):
    masks=[mask for mask in range(1<<n)if not any(mask&(1<<i)and mask&(1<<j)for i,j in edges)]
    graph=[[b for b,y in enumerate(masks)if(x^y).bit_count()==1]for x in masks]
    return masks,graph

def audit_additive_characterization():
    n=4;all_edges=list(itertools.combinations(range(n),2));checked=0
    for edge_mask in range(1<<len(all_edges)):
        edges=[e for bit,e in enumerate(all_edges)if edge_mask&(1<<bit)];masks,graph=independent_set_graph(n,edges)
        for weights in itertools.product((-1,0,1),repeat=n):
            raw=[sum(weights[i]for i in range(n)if mask&(1<<i))for mask in masks]
            lo,hi=min(raw),max(raw);values=[(r-lo)/(hi-lo)if hi>lo else 0.5 for r in raw]
            result=filtration(graph,values)
            positive_edge=any(weights[i]>0 and weights[j]>0 for i,j in edges)
            assert(result['positive_bar_count']>0)==positive_edge,(edges,weights,result)
            checked+=1
    return {'graph_vertices':n,'pair_graphs':64,'weight_vectors':81,'graph_weight_cases':checked}

def main():
    if not __debug__:raise RuntimeError('assertions must remain enabled')
    import argparse;parser=argparse.ArgumentParser();parser.add_argument('--outdir');args=parser.parse_args();root=Path(args.outdir)if args.outdir else Path(__file__).resolve().parents[3]/'results'/'expansion'/'verification';root.mkdir(parents=True,exist_ok=True)
    result={'filtration_fixtures':audit_filtration(),'conditional_null':audit_conditional(),'additive_superlevel_characterization':audit_additive_characterization(),'status':'all exact assertions passed'}
    (root/'score_controls_review_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
