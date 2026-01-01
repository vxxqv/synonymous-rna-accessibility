"""Replay whole-threshold geometry on every completed context sensitivity."""
from pathlib import Path
import json,hashlib
from retention_filtration import filtration
ROOT=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    variants={(b['id'],b['sensitivity']):b for b in json.loads((ROOT/'data/sensitivity_windows.json').read_text())};rows=[]
    for p in sorted((ROOT/'results/sensitivity').glob('*.json')):
        r=json.loads(p.read_text());variant=variants[r['id'],r['sensitivity']];graph=[[] for _ in r['genotypes']]
        for i,j in r['edges']:graph[i].append(j);graph[j].append(i)
        rows.append({'id':r['id'],'transcript_id':r['transcript_id'],'gene':r['gene'],'condition':r['condition'],'sensitivity':r['sensitivity'],
          'selection_group':variant.get('selection_group','hash_selected'),'filtration':filtration(graph,r['pair_retention'],full=True),
          'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    out={'rows':rows,'filtration_source_sha256':hashlib.sha256((Path(__file__).parent/'retention_filtration.py').read_bytes()).hexdigest(),'summary':{}}
    for group in ('hash_selected','posthoc_exception'):
        for name in ('context210','context300','Turner1999'):
            rr=[r for r in rows if r['selection_group']==group and r['sensitivity']==name];out['summary'][group+'/'+name]={'graphs':len(rr),'positive_bar_graphs':sum(r['filtration']['positive_bar_count']>0 for r in rr),'mean_integrated_excess':sum(r['filtration']['integrated_excess'] for r in rr)/len(rr) if rr else None}
    (ROOT/'results/expansion/context_filtration.json').write_text(json.dumps(out,allow_nan=False));print(json.dumps(out['summary'],indent=2))
