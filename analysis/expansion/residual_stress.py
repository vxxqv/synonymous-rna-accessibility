"""Focused retrospective stress tests; no discovery result becomes validation."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from extension_inference import interval
from additive_score_control import design
from retention_filtration import filtration
from accessibility import components
def family(name):return [json.loads(p.read_text()) for p in sorted((ROOT/'results/expansion'/name).glob('*.json')) if p.name!='summary.json']
def influence(rows,value):
    groups={}
    for r in rows:groups.setdefault(r['transcript_id'],[]).append(value(r))
    ranked=sorted(groups,key=lambda g:-abs(sum(groups[g])))
    return {f'exclude_top{k}_absolute_cluster_contributions':{'excluded':ranked[:k],**interval([r for r in rows if r['transcript_id'] not in ranked[:k]],value)} for k in (1,5)}
def main():
    sm=family('smoothness');ad=family('additive');result={'scope':'retrospective discovery diagnostics; nominal cluster intervals, not independent validation','conditions':{}}
    for c in ('incell','cellfree'):
        s=[r for r in sm if r['condition']==c];a=[r for r in ad if r['condition']==c];d={};value=lambda r:r['heat']['integrated_excess']['mean']-r['observed']['integrated_excess']
        d['heat_influence']=influence(s,value);d['heat_roughness_sensitivities']={}
        for t in (.01,.05,.1):
            rows=[r for r in s if abs(r['relative_roughness_mismatch'])<=t]
            d['heat_roughness_sensitivities'][str(t)]=interval(rows,value) if rows else {'rows':0}
        for method in ('nucleotide','codon'):
            field=lambda r:r['approximations'][method]
            value=lambda r:field(r)['rank_restored']['integrated_excess']['mean']-r['observed']['integrated_excess']
            good=[r for r in a if field(r)['total_residual_df']>0]
            away=[r for r in a if field(r)['total_residual_df']/r['n_compatible']>=.1]
            entry={'exclude_saturated':interval(good,value),'residual_df_at_least10percent_n':interval(away,value),'influence':influence(a,value)}
            entry['99_sampled_tie_orders_lower_envelope']=interval(a,lambda r:field(r)['rank_restored']['integrated_excess']['minimum']-r['observed']['integrated_excess'])
            entry['99_sampled_tie_orders_upper_envelope']=interval(a,lambda r:field(r)['rank_restored']['integrated_excess']['maximum']-r['observed']['integrated_excess'])
            entry['graphs_with_sampled_tie_variation']=sum(field(r)['rank_restored']['integrated_excess']['maximum']-field(r)['rank_restored']['integrated_excess']['minimum']>1e-12 for r in a)
            d[method]=entry
        # Refit only to assess numerical tie-scale sensitivity, no new fitting claim.
        ties=[]
        for r in a:
            n=json.loads((ROOT/'results/native'/f"{r['id']}_{c}.json").read_text());values=np.array(n['pair_retention']);graph=[[] for _ in values]
            for i,j in n['edges']:graph[i].append(j);graph[j].append(i)
            labels,_=components(graph);groups=[[i for i in range(len(values)) if labels[i]==g] for g in sorted(set(labels.values()))];fields=np.empty((len(values),3));differences=0
            for ix in groups:
                x=design([n['genotypes'][i] for i in ix],True);p=x@np.linalg.lstsq(x,values[ix],rcond=None)[0];rounded=[np.rint(p/(t*max(1.,float(np.ptp(p))))) for t in (1e-12,1e-10,1e-8)]
                # Compare equality partitions, not numerically different scale labels.
                partitions=[tuple(tuple(np.flatnonzero(v==v[j])) for j in range(len(ix))) for v in rounded]
                differences+=int(partitions[0]!=partitions[1] or partitions[1]!=partitions[2])
                for k,v in enumerate(rounded):
                    order=np.lexsort((np.arange(len(ix)),v));local=np.empty(len(ix));local[order]=np.sort(values[ix]);fields[ix,k]=local
            ties.append({'id':r['id'],'transcript_id':r['transcript_id'],'changed_components':differences,'values':[filtration(graph,fields[:,k],labels)['integrated_excess']-r['observed']['integrated_excess'] for k in range(3)]})
        d['tie_scale']={'method':'codon OLS, deterministic genotype-index order inside ties; sensitivity only','graphs_with_changed_tie_partition':sum(r['changed_components']>0 for r in ties),'contrasts':{str(t):interval(ties,lambda r,k=k:r['values'][k]) for k,t in enumerate((1e-12,1e-10,1e-8))}}
        result['conditions'][c]=d
    source=Path(__file__);result['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();(ROOT/'results/expansion/residual_stress.json').write_text(json.dumps(result,indent=2,allow_nan=False));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
