"""Descriptive additive approximations, not held-out biological predictors.

OLS scores reveal how much local field variation is additive by nucleotide or
codon. Rank reassignment restores the observed component score multiset.
Random tie orders are explicitly summarized instead of chosen favorably.
"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,argparse,hashlib,itertools
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from provenance import signature as cache_signature,cache_valid
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from accessibility import components
from retention_filtration import filtration
from smoothness_controls import roughness

def design(sequences, codon=False):
    length=len(sequences[0]);width=3 if codon else 1;columns=[np.ones(len(sequences))]
    for pos in range(0,length,width):
        states=sorted(set(s[pos:pos+width] for s in sequences))
        columns.extend(np.array([s[pos:pos+width]==v for s in sequences],float) for v in states[1:])
    return np.array(columns).T

def evaluate(task):
    b,c=task;r=json.loads((ROOT/'results/native'/f"{b['id']}_{c}.json").read_text());seq=r['genotypes'];graph=[[] for _ in seq]
    for i,j in r['edges']:graph[i].append(j);graph[j].append(i)
    labels,_=components(graph);values=np.array(r['pair_retention']);groups=[[i for i in range(len(seq)) if labels[i]==v] for v in sorted(set(labels.values()))]
    seed=int(hashlib.sha256(f"additive:1729:{b['id']}:{c}".encode()).hexdigest()[:16],16);rng=np.random.default_rng(seed)
    out={k:r[k] for k in ('id','transcript_id','gene','locus_tag','condition','n_compatible')};out['observed']=filtration(graph,values,labels)
    out['approximations']={}
    for codon,name in ((False,'nucleotide'),(True,'codon')):
        prediction=np.zeros(len(seq));fits=[];rank_fields=np.zeros((len(seq),99))
        for ix in groups:
            x=design([seq[i] for i in ix],codon);y=values[ix];coeff,_,rank,_=np.linalg.lstsq(x,y,rcond=None);p=x@coeff;prediction[ix]=p
            total=float(np.sum((y-y.mean())**2));residual=float(np.sum((y-p)**2));fits.append({'n':len(ix),'rank':int(rank),'residual_df':int(len(ix)-rank),'sum_squares':total,'residual_sum_squares':residual,'R_squared':1-residual/total if total>1e-24 else None})
            tolerance=1e-10*max(1.,float(np.ptp(p)));rounded=np.rint(p/tolerance)
            for draw in range(99):
                order=np.lexsort((rng.random(len(ix)),rounded));f=np.empty(len(ix));f[order]=np.sort(y);rank_fields[ix,draw]=f
        null=[filtration(graph,rank_fields[:,i],labels) for i in range(99)]
        total=sum(f['sum_squares'] for f in fits);residual=sum(f['residual_sum_squares'] for f in fits)
        out['approximations'][name]={'in_sample_R_squared':1-residual/total if total>1e-24 else None,'component_fits':fits,
          'total_residual_df':sum(f['residual_df'] for f in fits),'raw_prediction_filtration':filtration(graph,prediction,labels),
          'rank_restored':{k:{'mean':float(np.mean([f[k] for f in null])),'minimum':float(min(f[k] for f in null)),'maximum':float(max(f[k] for f in null))} for k in ('integrated_excess','longest_lifetime')},
          'observed_roughness':float(roughness(r['edges'],values)),'rank_restored_mean_roughness':float(np.mean([roughness(r['edges'],rank_fields[:,i]) for i in range(99)]))}
    out['run_provenance']=cache_signature('additive_score_control.py',b,c,{'tie_draws':99,'seed':1729})
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=8);a=p.parse_args();blocks=json.loads((ROOT/'data/native_windows.json').read_text());out=ROOT/'results/expansion/additive';out.mkdir(parents=True,exist_ok=True)
    tasks=[(b,c) for b in blocks for c in ('incell','cellfree') if not cache_valid(out/f"{b['id']}_{c}.json",cache_signature('additive_score_control.py',b,c,{'tie_draws':99,'seed':1729}))]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for i,f in enumerate(as_completed([pool.submit(evaluate,t) for t in tasks]),1):
            r=f.result();(out/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False))
            if i%100==0:print('completed',i,flush=True)
    rows=[json.loads(p.read_text()) for p in out.glob('*.json') if p.name!='summary.json'];s={}
    for c in ('incell','cellfree'):
        rr=[r for r in rows if r['condition']==c];s[c]={'graphs':len(rr)}
        for method in ('nucleotide','codon'):
            fits=[r['approximations'][method] for r in rr];r2=[f['in_sample_R_squared'] for f in fits if f['in_sample_R_squared'] is not None]
            s[c][method]={'nonconstant_graphs':len(r2),'median_in_sample_R_squared':float(np.median(r2)),
              'saturated_fits':sum(f['total_residual_df']==0 for f in fits),
              'observed_mean_integrated_excess':float(np.mean([r['observed']['integrated_excess'] for r in rr])),
              'rank_restored_mean_integrated_excess':float(np.mean([f['rank_restored']['integrated_excess']['mean'] for f in fits])),
              'rank_restored_mean_longest_lifetime':float(np.mean([f['rank_restored']['longest_lifetime']['mean'] for f in fits]))}
    (out/'summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
