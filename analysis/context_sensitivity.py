"""Prespecified context and parameter sensitivities on hash-selected blocks."""
from pathlib import Path
import json,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from native_models import ROOT,read_dbn,evaluate_one,CONDITIONS

def prepare(data_dir):
    records=read_dbn(Path(data_dir)/'weeks194_orig/incell.dbn')
    blocks=json.loads((ROOT/'data/native_windows.json').read_text())
    chosen=[]
    for block in blocks:
        pos=block['transcript_pos0'];sequence=records[block['transcript_id']]['sequence']
        if pos-142>=0 and pos-142+300<=len(sequence):chosen.append(block)
        if len(chosen)==24:break
    variants=[]
    for block in chosen:
        sequence=records[block['transcript_id']]['sequence']
        for length in (210,300):
            start=block['transcript_pos0']-(length-15)//2
            item=dict(block);item['context']=sequence[start:start+length];item['context_start0']=start
            item['variable_start0']=block['transcript_pos0']-start
            item['touching_pairs']={c:[(a+block['context_start0']-start,b+block['context_start0']-start) for a,b in block['touching_pairs'][c]] for c in CONDITIONS}
            item['sensitivity']=f'context{length}';item['selection_group']='hash_selected';variants.append(item)
        if len([v for v in variants if v.get('model')=='Turner1999'])<16:
            item=dict(block);item['model']='Turner1999';item['sensitivity']='Turner1999';item['selection_group']='hash_selected';variants.append(item)
    (ROOT/'data/sensitivity_windows.json').write_text(json.dumps(variants,indent=2))
    print('variants',len(variants))

def add_exceptions(data_dir):
    records=read_dbn(Path(data_dir)/'weeks194_orig/incell.dbn')
    blocks={b['id']:b for b in json.loads((ROOT/'data/native_windows.json').read_text())}
    variants=json.loads((ROOT/'data/sensitivity_windows.json').read_text())
    existing={(b['id'],b['sensitivity']) for b in variants}
    cases=json.loads((ROOT/'results/local_gate_restoration.json').read_text())
    for case in cases:
        block=blocks[case['id']];sequence=records[block['transcript_id']]['sequence']
        for length in (210,300):
            name=f'context{length}'
            if (block['id'],name) in existing:continue
            start=block['transcript_pos0']-(length-15)//2
            if start<0 or start+length>len(sequence):continue
            item=dict(block);item['context']=sequence[start:start+length];item['context_start0']=start
            item['variable_start0']=block['transcript_pos0']-start
            item['touching_pairs']={c:[(a+block['context_start0']-start,b+block['context_start0']-start) for a,b in block['touching_pairs'][c]] for c in CONDITIONS}
            item['sensitivity']=name;item['selection_group']='posthoc_exception';variants.append(item);existing.add((block['id'],name))
        name='Turner1999'
        if (block['id'],name) not in existing:
            item=dict(block);item['model']=name;item['sensitivity']=name;item['selection_group']='posthoc_exception';variants.append(item);existing.add((block['id'],name))
    (ROOT/'data/sensitivity_windows.json').write_text(json.dumps(variants,indent=2))
    print('variants after exception audit',len(variants))

def task(args):
    block,condition=args
    result=evaluate_one((block,condition));result['sensitivity']=block['sensitivity'];return result

def run(workers):
    variants=json.loads((ROOT/'data/sensitivity_windows.json').read_text())
    out=ROOT/'results/sensitivity';out.mkdir(exist_ok=True)
    tasks=[(b,c) for b in variants for c in CONDITIONS[:2] if not(out/f"{b['id']}_{c}_{b['sensitivity']}.json").exists()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,f in enumerate(as_completed([pool.submit(task,t) for t in tasks]),1):
            r=f.result();(out/f"{r['id']}_{r['condition']}_{r['sensitivity']}.json").write_text(json.dumps(r))
            if i%10==0:print('sensitivity',i,'of',len(tasks),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','exceptions','run']);parser.add_argument('--data-dir');parser.add_argument('--workers',type=int,default=10);a=parser.parse_args()
    if a.action=='prepare':prepare(a.data_dir)
    elif a.action=='exceptions':add_exceptions(a.data_dir)
    else:run(a.workers)
