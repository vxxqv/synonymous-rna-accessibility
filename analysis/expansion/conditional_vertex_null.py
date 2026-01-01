"""Conditional exact/sampled site-percolation nulls on frozen native graphs.

The number of retained vertices in every compatibility component is fixed.
The optional degree stratification also fixes coarse local degree composition.
These are model-geometry controls, not tests of natural selection.
"""
from pathlib import Path
import sys, json, math, itertools, argparse, hashlib, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from provenance import signature as cache_signature,cache_valid
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from accessibility import components

def degree_bin(value):
    return 0 if value==0 else 1 if value==1 else 2 if value<4 else 3 if value<8 else 4

def numerator(graph, keep):
    _, sizes=components(graph,keep)
    return sum(s*s for s in sizes)

def conditional_null(graph, labels, keep, degree_matched, seed, draws=999, cap=2000):
    groups={}
    for i in range(len(graph)):
        key=(labels[i],degree_bin(len(graph[i])) if degree_matched else 0)
        groups.setdefault(key,[]).append(i)
    strata=[(vertices,sum(i in keep for i in vertices)) for vertices in groups.values()]
    population=math.prod(math.comb(len(v),k) for v,k in strata)
    exact=population<=cap
    fixed={i for vertices,k in strata if k==len(vertices) for i in vertices}
    variable=[(vertices,k) for vertices,k in strata if 0<k<len(vertices)]
    denominator=sum(sum(i in keep for i in range(len(graph)) if labels[i]==c)**2 for c in set(labels.values()))
    observed=numerator(graph,keep)
    if exact:
        subsets=itertools.product(*(itertools.combinations(v,k) for v,k in variable))
        nums=[numerator(graph,fixed.union(*(set(s) for s in choice))) for choice in subsets]
    else:
        rng=np.random.default_rng(seed)
        nums=[]
        for _ in range(draws):
            sample=fixed.copy()
            for vertices,k in variable:sample.update(rng.choice(vertices,k,replace=False).tolist())
            nums.append(numerator(graph,sample))
    values=1-np.array(nums,dtype=float)/denominator
    observed_excess=1-observed/denominator
    tail=sum(n>=observed for n in nums)
    p=tail/len(nums) if exact else (tail+1)/(len(nums)+1)
    return {'exact':exact,'population_size':str(population),'draws':len(nums),
       'observed_excess':observed_excess,'null_mean':float(values.mean()),
       'null_sd':float(values.std(ddof=0)),'null_min':float(values.min()),'null_max':float(values.max()),
       'null_quantiles':np.quantile(values,[.025,.5,.975]).tolist(),
       'statistic_varied':len(set(nums))>1,'lower_tail_p':p,
       'mc_se_mean':None if exact else float(values.std(ddof=1)/math.sqrt(len(nums))),
       'connection_denominator':str(denominator),'observed_connection_numerator':str(observed)}

def evaluate(task):
    block, condition,draws=task
    r=json.loads((ROOT/'results/native'/f"{block['id']}_{condition}.json").read_text())
    graph=[[] for _ in r['genotypes']]
    for a,b in r['edges']:graph[a].append(b);graph[b].append(a)
    labels,_=components(graph);keep={i for i,value in enumerate(r['pair_retention']) if value>=.5}
    out={k:r[k] for k in ('id','transcript_id','gene','locus_tag','condition','n_compatible')}
    out['n_retained']=len(keep)
    out['run_provenance']=cache_signature('conditional_vertex_null.py',block,condition,{'draws':draws,'cap':2000,'seed':1729,'cutoff':.5})
    if not keep:out['status']='empty';return out
    out['status']='evaluated'
    base_seed=int(hashlib.sha256(f"1729:{block['id']}:{condition}".encode()).hexdigest()[:16],16)
    for j,name in enumerate(('uniform','degree')):
        out[name]=conditional_null(graph,labels,keep,j==1,base_seed+j,draws=draws)
        assert abs(out[name]['observed_excess']-r['thresholds']['0.5']['excess'])<1e-12
    values=np.array(r['pair_retention'])
    observed=sum((values[a]-values[b])**2 for a,b in r['edges'])
    expected=0.
    for c in set(labels.values()):
        ix=[i for i in range(len(graph)) if labels[i]==c]
        edges=sum(len(graph[i]) for i in ix)/2
        if len(ix)>1:expected+=edges*2*len(ix)/(len(ix)-1)*np.var(values[ix])
    out['edge_retention_roughness']=observed
    out['permuted_expected_roughness']=float(expected)
    out['roughness_ratio']=float(observed/expected) if expected else None
    return out

def summarize(rows):
    out={}
    for condition in ('incell','cellfree'):
        use=[r for r in rows if r['condition']==condition and r['status']=='evaluated']
        out[condition]={'nonempty':len(use),'empty':sum(r['condition']==condition and r['status']=='empty' for r in rows)}
        for method in ('uniform','degree'):
            rr=[r[method] for r in use];informative=[r for r in rr if r['statistic_varied']]
            out[condition][method]={'statistic_varied':len(informative),'exact':sum(r['exact'] for r in rr),
              'mean_observed_excess':float(np.mean([r['observed_excess'] for r in rr])) if rr else None,
              'mean_null_excess':float(np.mean([r['null_mean'] for r in rr])) if rr else None,
              'mean_null_minus_observed':float(np.mean([r['null_mean']-r['observed_excess'] for r in rr])) if rr else None,
              'informative_mean_contrast':float(np.mean([r['null_mean']-r['observed_excess'] for r in informative])) if informative else None}
        rough=[r['roughness_ratio'] for r in use if r['roughness_ratio'] is not None]
        out[condition]['median_roughness_ratio']=float(np.median(rough)) if rough else None
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',type=int,default=0);p.add_argument('--workers',type=int,default=8);p.add_argument('--draws',type=int,default=999);args=p.parse_args()
    blocks=json.loads((ROOT/'data/native_windows.json').read_text());blocks=blocks[:args.pilot] if args.pilot else blocks
    target=ROOT/'results/expansion'/('conditional_null_pilot' if args.pilot else 'conditional_null');target.mkdir(parents=True,exist_ok=True)
    tasks=[(b,c,args.draws) for b in blocks for c in ('incell','cellfree') if not cache_valid(target/f"{b['id']}_{c}.json",cache_signature('conditional_vertex_null.py',b,c,{'draws':args.draws,'cap':2000,'seed':1729,'cutoff':.5}))]
    print('pending',len(tasks),flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,f in enumerate(as_completed([pool.submit(evaluate,t) for t in tasks]),1):
            r=f.result();(target/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,indent=2,allow_nan=False))
            if i%50==0:print('completed',i,flush=True)
    rows=[json.loads(p.read_text()) for p in sorted(target.glob('*.json')) if p.name!='summary.json']
    summary=summarize(rows);(target/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print(json.dumps(summary,indent=2))
