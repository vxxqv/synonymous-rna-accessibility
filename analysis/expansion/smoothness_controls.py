"""Generic correlation-aware reference fields and exact conditional controls.

No reference model here establishes evolutionary selection or validates a
folding engine. The heat-rank reference is not a uniform conditional null.
"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,itertools,math,hashlib,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from provenance import signature as cache_signature,cache_valid
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from accessibility import components
from retention_filtration import filtration

def roughness(edges, values):
    return sum((values[i]-values[j])**2 for i,j in edges)

def rank_assign(fields, sorted_values):
    order=np.argsort(fields,axis=0,kind='stable');out=np.empty_like(fields)
    out[order,np.arange(fields.shape[1])]=np.asarray(sorted_values)[:,None]
    return out

def smooth_component(graph, indices, values, rng_fit, rng_test, draws):
    n=len(indices);ix={v:i for i,v in enumerate(indices)}
    edges=[(ix[i],ix[j]) for i in indices for j in graph[i] if j>i]
    observed=np.array(values)[indices];target=roughness(edges,observed);sorted_values=np.sort(observed)
    if n==1 or np.ptp(observed)==0:
        return np.repeat(observed[:,None],draws,axis=1),{'n':n,'status':'constant','target_roughness':float(target),'time':0.,'evaluation_roughness':float(target)}
    lap=np.zeros((n,n))
    for i,j in edges:lap[i,i]+=1;lap[j,j]+=1;lap[i,j]-=1;lap[j,i]-=1
    lap/=2*len(edges)/n
    eig,q=np.linalg.eigh(lap);good=eig>1e-9;eig=eig[good];q=q[:,good]
    if len(eig)!=n-1:raise AssertionError('Expected connected component')
    # A common positive rescaling leaves ranks unchanged and prevents underflow.
    eig=eig-eig[0];zfit=rng_fit.standard_normal((len(eig),64));ztest=rng_test.standard_normal((len(eig),draws))
    edge_array=np.array(edges)
    def generate(t,z):return rank_assign(q@(np.exp(-t*eig/2)[:,None]*z),sorted_values)
    def mean_rough(t):
        f=generate(t,zfit);return float(np.mean(np.sum((f[edge_array[:,0]]-f[edge_array[:,1]])**2,axis=0)))
    times=[0.,.125,.25,.5,1.,2.,4.,8.,16.,32.,64.,128.,256.]
    fits=[mean_rough(t) for t in times]
    best=min(range(len(times)),key=lambda i:abs(fits[i]-target))
    # Refine the nearest bracket using calibration draws only.
    lo=times[max(0,best-1)];hi=times[min(len(times)-1,best+1)]
    candidates=list(zip(times,fits))
    for _ in range(8):
        mid=(lo+hi)/2;value=mean_rough(mid);candidates.append((mid,value))
        if value>target:lo=mid
        else:hi=mid
    t,fit=min(candidates,key=lambda r:abs(r[1]-target));sample=generate(t,ztest)
    rr=np.sum((sample[edge_array[:,0]]-sample[edge_array[:,1]])**2,axis=0)
    return sample,{'n':n,'status':'fitted','time':t,'target_roughness':float(target),'calibration_roughness':fit,
       'evaluation_roughness':float(rr.mean()),'evaluation_roughness_se':float(rr.std(ddof=1)/math.sqrt(draws)),
       'relative_evaluation_mismatch':float((rr.mean()-target)/target) if target else None,
       'evaluation_roughness_quantiles':np.quantile(rr,[.025,.5,.975]).tolist()}

def exact_conditional(graph,labels,values,cap=40320):
    groups=[[i for i in range(len(graph)) if labels[i]==c] for c in sorted(set(labels.values()))]
    population=math.prod(math.factorial(len(ix)) for ix in groups)
    if population>cap:return {'status':'resource_gate','label_permutations':str(population)}
    edges=[(i,j) for i,nb in enumerate(graph) for j in nb if i<j];target=roughness(edges,values)
    results={str(t):[] for t in (.05,.1,.2)}
    for choices in itertools.product(*(itertools.permutations(values[ix]) for ix in groups)):
        permuted=np.array(values)
        for ix,v in zip(groups,choices):permuted[ix]=v
        r=roughness(edges,permuted);difference=abs(r-target)
        accepted=[t for t in (.05,.1,.2) if difference<=t*max(target,1e-14)+1e-14]
        if not accepted:continue
        f=filtration(graph,permuted,labels)
        for t in accepted:results[str(t)].append((f['integrated_excess'],f['longest_lifetime'],r))
    out={'status':'exact','label_permutations':str(population),'target_roughness':float(target),'tolerances':{}}
    for t,rr in results.items():
        if not rr:raise AssertionError('Original assignment must be in conditional population')
        x=np.array(rr);out['tolerances'][t]={'accepted':len(rr),'integrated_excess_mean':float(x[:,0].mean()),'longest_lifetime_mean':float(x[:,1].mean()),'integrated_excess_varied':len(set(x[:,0]))>1,
          'integrated_excess_range':[float(x[:,0].min()),float(x[:,0].max())]}
    return out

def evaluate(task):
    block,c,draws=task;r=json.loads((ROOT/'results/native'/f"{block['id']}_{c}.json").read_text());graph=[[] for _ in r['genotypes']]
    for i,j in r['edges']:graph[i].append(j);graph[j].append(i)
    labels,_=components(graph);values=np.array(r['pair_retention']);groups=[[i for i in range(len(graph)) if labels[i]==v] for v in sorted(set(labels.values()))]
    seed=int(hashlib.sha256(f"smooth:1729:{block['id']}:{c}".encode()).hexdigest()[:16],16)
    rng_fit=np.random.default_rng(seed);rng_test=np.random.default_rng(seed+1);fields=np.zeros((len(graph),draws));fits=[]
    for indices in groups:
        fields[indices,:],fit=smooth_component(graph,indices,values,rng_fit,rng_test,draws);fits.append(fit)
    out={k:r[k] for k in ('id','transcript_id','gene','locus_tag','condition','n_compatible')};out['observed']=filtration(graph,values,labels,True);out['heat_fits']=fits;out['heat_draws']=draws
    null=[filtration(graph,fields[:,i],labels) for i in range(draws)]
    out['heat']={k:{'mean':float(np.mean([v[k] for v in null])),'quantiles':np.quantile([v[k] for v in null],[.025,.5,.975]).tolist()} for k in ('integrated_excess','longest_lifetime','positive_bar_count')}
    target=sum(f['target_roughness'] for f in fits);sample=sum(f['evaluation_roughness'] for f in fits)
    out['roughness_target']=target;out['roughness_sample_mean']=sample;out['relative_roughness_mismatch']=(sample-target)/target if target else 0.
    out['exact_conditional']=exact_conditional(graph,labels,values)
    out['run_provenance']=cache_signature('smoothness_controls.py',block,c,{'draws':draws,'fit_draws':64,'cap':40320,'seed':1729})
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',type=int,default=0);p.add_argument('--workers',type=int,default=8);p.add_argument('--draws',type=int,default=199);a=p.parse_args();blocks=json.loads((ROOT/'data/native_windows.json').read_text());blocks=blocks[:a.pilot] if a.pilot else blocks
    out=ROOT/'results/expansion'/('smoothness_pilot' if a.pilot else 'smoothness');out.mkdir(parents=True,exist_ok=True);tasks=[(b,c,a.draws) for b in blocks for c in ('incell','cellfree') if not cache_valid(out/f"{b['id']}_{c}.json",cache_signature('smoothness_controls.py',b,c,{'draws':a.draws,'fit_draws':64,'cap':40320,'seed':1729}))]
    print('pending',len(tasks),flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for i,f in enumerate(as_completed([pool.submit(evaluate,t) for t in tasks]),1):
            r=f.result();(out/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False))
            if i%25==0:print('completed',i,flush=True)
    rows=[json.loads(p.read_text()) for p in out.glob('*.json') if p.name!='summary.json'];s={}
    for c in ('incell','cellfree'):
        rr=[r for r in rows if r['condition']==c];exact=[r for r in rr if r['exact_conditional']['status']=='exact'];s[c]={'graphs':len(rr),'roughness_mismatch_gt10percent':sum(abs(r['relative_roughness_mismatch'])>.1 for r in rr),'exact_conditional_graphs':len(exact)}
        for k in ('integrated_excess','longest_lifetime'):
            s[c][k]={'observed_mean':float(np.mean([r['observed'][k] for r in rr])),'heat_reference_mean':float(np.mean([r['heat'][k]['mean'] for r in rr]))}
        s[c]['exact_conditional_tolerances']={t:{'graphs':len(exact),'varied':sum(r['exact_conditional']['tolerances'][t]['integrated_excess_varied'] for r in exact),'observed_mean_integrated_excess':float(np.mean([r['observed']['integrated_excess'] for r in exact])) if exact else None,'reference_mean_integrated_excess':float(np.mean([r['exact_conditional']['tolerances'][t]['integrated_excess_mean'] for r in exact])) if exact else None} for t in ('0.05','0.1','0.2')}
    (out/'summary.json').write_text(json.dumps(s,indent=2,allow_nan=False));print(json.dumps(s,indent=2))
