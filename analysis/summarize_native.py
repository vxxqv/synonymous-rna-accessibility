"""Export exact native-model results and covariates for independent R inference."""
from pathlib import Path
import json,csv,math
ROOT=Path(__file__).resolve().parents[1]
blocks=json.loads((ROOT/'data/native_windows.json').read_text())
rows=[];paired=[]
for b in blocks:
    rr={}
    for condition in ('incell','cellfree','kasugamycin'):
        r=json.loads((ROOT/'results/native'/f"{b['id']}_{condition}.json").read_text())
        assert r['id']==b['id'] and r['gene']==b['gene']
        pairs=b['touching_pairs'][condition]
        f=r['thresholds']['0.5']
        row={k:b[k] for k in ('id','transcript_id','gene','locus_tag','peptide','cartesian_size')}
        row.update(condition=condition,n_pairs=len(pairs),mean_span=sum(j-i for i,j in pairs)/len(pairs),
            n_compatible=r['n_compatible'],a_compatible=r['a_compatible'],
            effective_components=1/r['a_compatible'],log_effective_components=-math.log(r['a_compatible']),
            native_retention=r['native_retention'],n_filtered=f['n'],excess=f['excess'],
            n_filtered_components=f.get('components'),mean_retention=sum(r['pair_retention'])/len(r['pair_retention']),
            coverage=b['coverage'][condition],n_internal_pairs=r['n_internal_pairs'],
            pair_phase00=sum((i-b['variable_start0'])%3==0 and (j-b['variable_start0'])%3==0 for i,j in pairs))
        rows.append(row);rr[condition]=row
    a,z=rr['incell'],rr['cellfree']
    pair={k:b[k] for k in ('id','transcript_id','gene','locus_tag')}
    pair.update(d_log_effective=a['log_effective_components']-z['log_effective_components'],
        d_log_abundance=math.log(a['n_compatible'])-math.log(z['n_compatible']),
        d_pairs=a['n_pairs']-z['n_pairs'],d_span=a['mean_span']-z['mean_span'],
        d_internal=a['n_internal_pairs']-z['n_internal_pairs'],
        d_excess=None if a['excess'] is None or z['excess'] is None else a['excess']-z['excess'],
        d_retention=a['native_retention']-z['native_retention'],
        incell_nonempty=a['n_filtered']>0,cellfree_nonempty=z['n_filtered']>0)
    paired.append(pair)
for name,data in [('native_metrics.csv',rows),('native_paired.csv',paired)]:
    with (ROOT/'results'/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
summary={}
for condition in ('incell','cellfree','kasugamycin'):
    subset=[r for r in rows if r['condition']==condition]
    nonempty=[r for r in subset if r['n_filtered']]
    summary[condition]={'n_blocks':len(subset),'n_compatible_genotypes':sum(r['n_compatible'] for r in subset),
        'nonempty_filtered_blocks':len(nonempty),'blocks_with_excess':sum(r['excess']>1e-12 for r in nonempty),
        'mean_excess_nonempty':sum(r['excess'] for r in nonempty)/len(nonempty) if nonempty else None,
        'total_filtered_genotypes':sum(r['n_filtered'] for r in subset)}
summary['cohort']={'blocks':len(blocks),'transcripts':len(set(r['transcript_id'] for r in blocks)),
    'unique_loci':len(set(r['locus_tag'] for r in blocks)),
    'cartesian_per_condition':sum(r['cartesian_size'] for r in blocks)}
(ROOT/'results/native_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
