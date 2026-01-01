from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text())
variants=load(ROOT/'data/sensitivity_windows.json');cases=load(ROOT/'results/local_gate_restoration.json')
original={(r['id'],r['condition']) for r in cases}
out={'hash_selected':{},'exceptions':{}}
for name in ('context210','context300','Turner1999'):
    rows=[];extra=[]
    for b in variants:
        if b['sensitivity']!=name:continue
        for c in ('incell','cellfree'):
            r=load(ROOT/'results/sensitivity'/f"{b['id']}_{c}_{name}.json")
            base=load(ROOT/'results/native'/f"{b['id']}_{c}.json")
            assert r['n_compatible']==base['n_compatible'] and r['component_sizes']==base['component_sizes']
            f=r['thresholds']['0.5']
            if b.get('selection_group','hash_selected')=='hash_selected':rows.append(f)
            if (b['id'],c) in original:extra.append(f)
    non=[r for r in rows if r['n']]
    out['hash_selected'][name]={'condition_windows':len(rows),'nonempty':len(non),'with_extra':sum(r['excess']>1e-12 for r in non),'mean_excess_nonempty':sum(r['excess'] for r in non)/len(non) if non else None}
    non=[r for r in extra if r['n']]
    out['exceptions'][name]={'tested_original_cases':len(extra),'nonempty':len(non),'fragmentation_persists':sum(r['excess']>1e-12 for r in non),'empty':len(extra)-len(non)}
(ROOT/'results/sensitivity_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
