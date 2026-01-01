"""Pilot longer COMPLETE synonymous spaces using codon-CSP pruning."""
from pathlib import Path
import sys,json,csv,math,argparse,time
from collections import defaultdict,Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'analysis'))
from accessibility import *
from constrained_enumeration import enumerate_context
from native_models import read_dbn,read_react
from retention_filtration import filtration

def prepare(data_dir):
    data_dir=Path(data_dir);records={c:read_dbn(data_dir/'weeks194_orig'/f'{c}.dbn') for c in ('incell','cellfree')};reacts={c:read_react(data_dir/'weeks194_orig'/f'{c}.react') for c in records}
    overlaps=defaultdict(set)
    for row in csv.DictReader((data_dir/'coding_position_map.csv').open()):
        if row['overlaps_other_CDS']=='1':overlaps[row['transcript_id']].add(int(row['transcript_pos'])-1)
    cds={(r['transcript_id'],r['locus_tag']):r for r in csv.DictReader((data_dir/'cds_mapping.csv').open())};blocks=json.loads((ROOT/'data/native_windows.json').read_text());variants=[];gates=[]
    for codons,number in ((7,12),(9,8),(11,4)):
        counts=Counter();tested=[]
        for block in blocks:
            tid=block['transcript_id'];pos=block['transcript_pos0']-3*((codons-5)//2);row=cds[tid,block['locus_tag']];end=pos+3*codons;seq=records['incell'][tid]['sequence']
            if pos<int(row['transcript_start'])-1+3 or end>int(row['transcript_end'])-3:counts['CDS_boundary']+=1;continue
            if any(i in overlaps[tid] for i in range(pos,end)):counts['overlap']+=1;continue
            if any(sum(v is not None for v in reacts[c][tid][pos:end])<math.ceil(.8*3*codons) for c in records):counts['coverage']+=1;continue
            # Same 150-nt context as the frozen central block. New target pairs
            # are all pairs touching the expanded interval, not the old subset.
            cs=block['context_start0'];ce=cs+len(block['context']);pp={};bad=False
            for c in records:
                pairs=sorted({(min(i,records[c][tid]['partners'][i]),max(i,records[c][tid]['partners'][i])) for i in range(pos,end) if i in records[c][tid]['partners']})
                if len(pairs)<3 or any(a<cs or b>=ce or b-a<4 for a,b in pairs):bad=True;break
                pp[c]=[(a-cs,b-cs) for a,b in pairs]
            if bad:counts['pair_containment']+=1;continue
            item={**block,'native':seq[pos:end],'peptide':translate(seq[pos:end]),'variable_start0':pos-cs,'transcript_pos0':pos,'codons':codons,'touching_pairs':pp}
            item['cartesian_size']=math.prod(len(CODONS[a]) for a in item['peptide']);item['enumeration']={};spaces={}
            try:
                for c in records:spaces[c],item['enumeration'][c]=enumerate_context(item['peptide'],item['context'],item['variable_start0'],pp[c])
            except MemoryError as e:counts['resource_limit']+=1;tested.append({'id':block['id'],'status':'resource_limit','reason':str(e)});continue
            if any(item['native'] not in spaces[c] for c in spaces):raise AssertionError('Native is not compatible')
            item['genotypes']=spaces;variants.append(item);counts['selected']+=1
            if counts['selected']==number:break
        gates.append({'codons':codons,'requested':number,'selected':counts['selected'],'sequential_rejections':dict(counts),'resource_failures':tested})
    (ROOT/'data/longer_windows.json').write_text(json.dumps(variants,indent=2));(ROOT/'results/expansion/longer_selection.json').write_text(json.dumps(gates,indent=2));print(json.dumps(gates,indent=2))

def evaluate(task):
    import RNA
    b,c=task
    if not RNA.params_load_RNA_Turner2004():raise RuntimeError('Parameters not loaded')
    md=RNA.md();md.temperature=37;md.dangles=2;md.noLP=0
    seqs=b['genotypes'][c];v=b['variable_start0'];values=[];start=time.time()
    for s in seqs:
        full=b['context'][:v]+s+b['context'][v+len(s):];fc=RNA.fold_compound(full,md);_,e=fc.mfe();fc.exp_params_rescale(e);fc.pf();p=fc.bpp();values.append(sum(p[a+1][d+1] for a,d in b['touching_pairs'][c])/len(b['touching_pairs'][c]))
    graph=adjacency(seqs);_,sizes=components(graph);internal=[(a-v,d-v) for a,d in b['touching_pairs'][c] if v<=a<d<v+len(b['native'])]
    assert sorted(Counter(signature(s,internal) for s in seqs).values())==sorted(sizes)
    return {'id':b['id'],'gene':b['gene'],'transcript_id':b['transcript_id'],'codons':b['codons'],'condition':c,'cartesian_size':b['cartesian_size'],
     'n_compatible':len(seqs),'genotypes':seqs,'edges':[(i,j) for i,nb in enumerate(graph) for j in nb if i<j],'pair_retention':values,'component_sizes':sorted(sizes),
     'enumeration':b['enumeration'][c],'thresholds':{str(t):excess_fragmentation(seqs,internal,graph,[i for i,x in enumerate(values) if x>=t]) for t in (.1,.25,.5,.75)},
     'filtration':filtration(graph,values,full=True),'seconds':time.time()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','verify']);p.add_argument('--data-dir');p.add_argument('--workers',type=int,default=8);a=p.parse_args();(ROOT/'results/expansion').mkdir(exist_ok=True)
    if a.action=='prepare':prepare(a.data_dir)
    elif a.action=='verify':
        n=0
        for b in json.loads((ROOT/'data/native_windows.json').read_text()):
            for c in ('incell','cellfree'):
                seqs,stats=enumerate_context(b['peptide'],b['context'],b['variable_start0'],b['touching_pairs'][c]);old=json.loads((ROOT/'results/native'/f"{b['id']}_{c}.json").read_text());assert seqs==sorted(old['genotypes']);n+=1
        (ROOT/'results/expansion/constrained_enumerator_check.json').write_text(json.dumps({'complete_graphs_checked':n,'all_genotype_sets_equal':True}));print(n,'sets match')
    else:
        out=ROOT/'results/expansion/longer';out.mkdir(exist_ok=True);blocks=json.loads((ROOT/'data/longer_windows.json').read_text());tasks=[(b,c) for b in blocks for c in ('incell','cellfree') if not (out/f"{b['id']}_{b['codons']}_{c}.json").exists()]
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for f in as_completed([pool.submit(evaluate,t) for t in tasks]):
                r=f.result();(out/f"{r['id']}_{r['codons']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False));print('completed',r['codons'],r['id'],r['condition'],r['n_compatible'],flush=True)
