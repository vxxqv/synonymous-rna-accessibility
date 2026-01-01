"""Transcript-cluster uncertainty for exploratory model-geometry contrasts."""
from pathlib import Path
import json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

def interval(rows, value, seed=1729):
    groups={}
    for r in rows:groups.setdefault(r['transcript_id'],[]).append(value(r))
    sums=np.array([sum(v) for v in groups.values()]);counts=np.array([len(v) for v in groups.values()]);n=len(groups)
    rng=np.random.default_rng(seed);weights=rng.multinomial(n,np.full(n,1/n),size=10000)
    estimates=weights@sums/(weights@counts)
    return {'estimate':float(sums.sum()/counts.sum()),'lower':float(np.quantile(estimates,.025)),
       'upper':float(np.quantile(estimates,.975)),'rows':len(rows),'transcript_clusters':n,'resamples':10000}

if __name__=='__main__':
    out={}
    for source in ('conditional_null','filtration','smoothness','additive'):
        rows=[json.loads(p.read_text()) for p in (ROOT/'results/expansion'/source).glob('*.json') if p.name!='summary.json']
        for condition in ('incell','cellfree'):
            use=[r for r in rows if r['condition']==condition and r.get('status','evaluated')=='evaluated']
            for name in (('uniform','degree') if source=='conditional_null' else (('nucleotide','codon') if source=='additive' else ('longest_lifetime','integrated_excess'))):
                value=(lambda r:r[name]['null_mean']-r[name]['observed_excess']) if source=='conditional_null' else ((lambda r:r['heat'][name]['mean']-r['observed'][name]) if source=='smoothness' else ((lambda r:r['approximations'][name]['rank_restored']['integrated_excess']['mean']-r['observed']['integrated_excess']) if source=='additive' else (lambda r:r['null'][name]['mean']-r['observed'][name])))
                out[f'{source}/{condition}/{name}/null_minus_observed']=interval(use,value)
                if source=='filtration':out[f'{source}/{condition}/{name}/observed']=interval(use,lambda r:r['observed'][name])
                groups=sorted(set(r['transcript_id'] for r in use));loo=[]
                for g in groups:
                    retained=[r for r in use if r['transcript_id']!=g]
                    loo.append({'excluded_transcript':g,'estimate':float(np.mean([value(r) for r in retained]))})
                original=float(np.mean([value(r) for r in use]))
                out[f'{source}/{condition}/{name}/leave_one_cluster_out']={'minimum':min(r['estimate'] for r in loo),'maximum':max(r['estimate'] for r in loo),'sign_changes':sum(r['estimate']*original<0 for r in loo),'zero_estimates':sum(r['estimate']==0 for r in loo),'details':loo}
            # A sensitivity chooses first SHA-ordered context per locus before aggregation.
            loci={}
            for r in sorted(use,key=lambda r:r['id']):loci.setdefault(r['locus_tag'],r)
            if source=='conditional_null':out[f'{source}/{condition}/degree/deduplicated']=interval(list(loci.values()),lambda r:r['degree']['null_mean']-r['degree']['observed_excess'])
            elif source=='filtration':out[f'{source}/{condition}/integrated_excess/deduplicated']=interval(list(loci.values()),lambda r:r['null']['integrated_excess']['mean']-r['observed']['integrated_excess'])
            elif source=='smoothness':
                out[f'{source}/{condition}/integrated_excess/deduplicated']=interval(list(loci.values()),lambda r:r['heat']['integrated_excess']['mean']-r['observed']['integrated_excess'])
                good=[r for r in use if abs(r['relative_roughness_mismatch'])<=.1]
                out[f'{source}/{condition}/integrated_excess/roughness_match_within10percent']=interval(good,lambda r:r['heat']['integrated_excess']['mean']-r['observed']['integrated_excess'])
                exact=[r for r in use if r['exact_conditional']['status']=='exact']
                for tol in ('0.05','0.1','0.2'):
                    out[f'{source}/{condition}/exact_conditional/{tol}']=interval(exact,lambda r:r['exact_conditional']['tolerances'][tol]['integrated_excess_mean']-r['observed']['integrated_excess'])
    (ROOT/'results/expansion/inference.json').write_text(json.dumps(out,indent=2,allow_nan=False));print(json.dumps({k:v for k,v in out.items() if k.startswith(('smoothness/','additive/')) and '/leave_one' not in k},indent=2))
