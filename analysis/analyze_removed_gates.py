"""Counterfactual restoration of certified local GU bridge vertices.

The endpoint population remains the filtered vertices, including after adding
intermediates. These certificates identify a sufficient graph mechanism; they
do not identify the energetic cause of a genotype's removal.
"""
from pathlib import Path
import sys,json,itertools
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from accessibility import components
blocks={b['id']:b for b in json.loads((ROOT/'data/native_windows.json').read_text())}
results=[]
for path in sorted((ROOT/'results/native').glob('*.json')):
    r=json.loads(path.read_text());f=r['thresholds']['0.5']
    if f['excess'] is None or f['excess']<=1e-12:continue
    b=blocks[r['id']];sequences=r['genotypes'];keep={i for i,x in enumerate(r['pair_retention']) if x>=.5}
    graph=[[] for s in sequences]
    for i,j in r['edges']:graph[i].append(j);graph[j].append(i)
    labels,_=components(graph,keep);v=b['variable_start0']
    pairs={tuple((i-v,j-v)) for i,j in b['touching_pairs'][r['condition']] if v<=i<j<v+15}
    restored=set();coding=set();generic=set();certificates=[]
    for z in set(range(len(sequences)))-keep:
        neighbors=[i for i in graph[z] if i in keep]
        for i,j in itertools.combinations(neighbors,2):
            if labels[i]==labels[j]:continue
            changed=tuple(k for k,(a,c) in enumerate(zip(sequences[i],sequences[j])) if a!=c)
            generic.add(z)
            if len(changed)==2 and changed[0]//3==changed[1]//3 and changed[0]%3==0 and changed[1]%3==2:
                coding.add(z)
            if changed in pairs and ''.join(sequences[z][k] for k in changed) in ('GU','UG'):
                restored.add(z);certificates.append({'removed':sequences[z],'endpoints':[sequences[i],sequences[j]],'pair_positions0':changed})
    newlabels,_=components(graph,keep|restored)
    counts=Counter(newlabels[i] for i in keep)
    connected=sum(n*n for n in counts.values())
    jointlabels,_=components(graph,keep|restored|coding)
    jointcounts=Counter(jointlabels[i] for i in keep)
    jointconnected=sum(n*n for n in jointcounts.values())
    genericlabels,_=components(graph,keep|generic)
    genericcounts=Counter(genericlabels[i] for i in keep)
    genericconnected=sum(n*n for n in genericcounts.values())
    original=int(f['same_component_pairs_filtered']);baseline=int(f['same_component_pairs_compatibility'])
    results.append({'id':r['id'],'condition':r['condition'],'gene':r['gene'],
        'n_endpoints':len(keep),'original_excess':f['excess'],
        'gu_bridge_vertices':len(restored),'certificates':certificates,
        'restored_connected_pairs':str(connected),'baseline_connected_pairs':str(baseline),
        'fraction_excess_recovered':(connected-original)/(baseline-original),
        'fully_explained_by_single_GU_addbacks':connected==baseline,
        'coding_gate_bridge_vertices':len(coding),'any_single_addback_vertices':len(generic),
        'fully_restored_by_GU_and_coding':jointconnected==baseline,
        'fully_restored_by_all_single_addbacks':genericconnected==baseline})
(ROOT/'results/local_gate_restoration.json').write_text(json.dumps(results,indent=2))
print(json.dumps({'cases':len(results),'any_local_gu_certificate':sum(r['gu_bridge_vertices']>0 for r in results),
    'fully_restored':sum(r['fully_explained_by_single_GU_addbacks'] for r in results)},indent=2))
