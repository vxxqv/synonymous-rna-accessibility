"""Exact zero-dimensional superlevel persistence of native graph filters.

Finite positive bars represent additional barriers inside compatibility
components. Roots surviving to the full graph are baseline components.
No persistence threshold, embedding or geometric interpolation is required.
"""
from pathlib import Path
import sys,json,itertools,argparse,hashlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from provenance import signature as cache_signature,cache_valid
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from accessibility import components

def filtration(graph, values, baseline=None, full=False):
    if baseline is None:baseline,_=components(graph)
    n=len(graph);parent=list(range(n));sizes=[1]*n;birth=list(map(float,values));active=[False]*n
    originals={};square_sum=0;denominator=0;bars=[];curve=[]
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    order=sorted(range(n),key=lambda i:(-values[i],i))
    groups=[(float(v),list(ix)) for v,ix in itertools.groupby(order,key=lambda i:values[i])]
    for index,(level,vertices) in enumerate(groups):
        for i in vertices:
            active[i]=True;square_sum+=1;c=baseline[i];old=originals.get(c,0);originals[c]=old+1;denominator+=2*old+1
        for i in vertices:
            for j in graph[i]:
                if not active[j]:continue
                a,b=find(i),find(j)
                if a==b:continue
                if (birth[a],-a)<(birth[b],-b):a,b=b,a
                if birth[b]>level:bars.append({'birth':birth[b],'death':level,'lifetime':birth[b]-level,'peak_vertex':b,'merge_edge':[i,j],'baseline_component':baseline[i]})
                square_sum+=2*sizes[a]*sizes[b];sizes[a]+=sizes[b];parent[b]=a
        curve.append({'threshold':level,'excess':1-square_sum/denominator,'retained':sum(originals.values()),'filtered_components':sum(active[i] and parent[i]==i for i in range(n))})
    area=sum(r['excess']*(r['threshold']-(curve[i+1]['threshold'] if i+1<len(curve) else 0)) for i,r in enumerate(curve))
    maximum=max((r['excess'] for r in curve),default=0)
    result={'positive_bar_count':len(bars),'longest_lifetime':max((b['lifetime'] for b in bars),default=0.),
      'integrated_excess':area,'maximum_excess':maximum,'baseline_components':len(set(baseline.values())),
      'threshold_at_maximum':next((r['threshold'] for r in curve if r['excess']==maximum),None)}
    if full:result.update({'bars':sorted(bars,key=lambda b:-b['lifetime']),'curve':curve})
    return result

def evaluate(task):
    block,condition,draws=task;r=json.loads((ROOT/'results/native'/f"{block['id']}_{condition}.json").read_text())
    graph=[[] for _ in r['genotypes']]
    for a,b in r['edges']:graph[a].append(b);graph[b].append(a)
    labels,_=components(graph);values=np.array(r['pair_retention']);observed=filtration(graph,values,labels,True)
    # Verify every distinct critical value against a separate BFS oracle.
    for row in observed['curve']:
        keep={i for i,v in enumerate(values) if v>=row['threshold']};_,sizes=components(graph,keep)
        counts=[sum(i in keep for i in range(len(graph)) if labels[i]==c) for c in set(labels.values())]
        assert abs(row['excess']-(1-sum(s*s for s in sizes)/sum(c*c for c in counts)))<1e-12
    groups=[[i for i in range(len(graph)) if labels[i]==c] for c in sorted(set(labels.values()))]
    rng=np.random.default_rng(int(hashlib.sha256(f"filtration:1729:{block['id']}:{condition}".encode()).hexdigest()[:16],16))
    null=[]
    for _ in range(draws):
        permuted=values.copy()
        for ix in groups:permuted[ix]=rng.permutation(values[ix])
        null.append(filtration(graph,permuted,labels))
    out={k:r[k] for k in ('id','transcript_id','gene','locus_tag','condition','n_compatible')};out['observed']=observed;out['null_draws']=draws
    out['run_provenance']=cache_signature('retention_filtration.py',block,condition,{'draws':draws,'seed':1729})
    out['null']={}
    for key in ('positive_bar_count','longest_lifetime','integrated_excess','maximum_excess'):
        x=np.array([q[key] for q in null]);o=observed[key]
        out['null'][key]={'mean':float(x.mean()),'sd':float(x.std()),'quantiles':np.quantile(x,[.025,.5,.975]).tolist(),'lower_tail_p':float((1+sum(x<=o))/(draws+1)),'varied':len(set(x))>1}
    return out

def summary(rows):
    out={}
    for c in ('incell','cellfree'):
        rr=[r for r in rows if r['condition']==c]
        out[c]={'graphs':len(rr),'graphs_with_positive_bar':sum(r['observed']['positive_bar_count']>0 for r in rr),
         'graphs_with_lifetime_above_005':sum(r['observed']['longest_lifetime']>.05 for r in rr),
         'graphs_with_lifetime_above_010':sum(r['observed']['longest_lifetime']>.1 for r in rr)}
        for key in ('longest_lifetime','integrated_excess','maximum_excess'):
            out[c][key]={'observed_mean':float(np.mean([r['observed'][key] for r in rr])),
             'null_mean':float(np.mean([r['null'][key]['mean'] for r in rr]))}
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',type=int,default=0);p.add_argument('--workers',type=int,default=8);p.add_argument('--draws',type=int,default=199);a=p.parse_args()
    blocks=json.loads((ROOT/'data/native_windows.json').read_text());blocks=blocks[:a.pilot] if a.pilot else blocks
    out=ROOT/'results/expansion'/('filtration_pilot' if a.pilot else 'filtration');out.mkdir(parents=True,exist_ok=True)
    tasks=[(b,c,a.draws) for b in blocks for c in ('incell','cellfree') if not cache_valid(out/f"{b['id']}_{c}.json",cache_signature('retention_filtration.py',b,c,{'draws':a.draws,'seed':1729}))]
    print('pending',len(tasks),flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for i,f in enumerate(as_completed([pool.submit(evaluate,t) for t in tasks]),1):
            r=f.result();(out/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False))
            if i%50==0:print('completed',i,flush=True)
    rows=[json.loads(p.read_text()) for p in out.glob('*.json') if p.name!='summary.json'];s=summary(rows);(out/'summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
