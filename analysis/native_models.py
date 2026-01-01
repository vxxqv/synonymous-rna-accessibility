"""Select complete native coding blocks and evaluate unrestrained folding models."""
from pathlib import Path
import sys,json,csv,hashlib,math,itertools,time,os,argparse
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from accessibility import *

CONDITIONS=('incell','cellfree','kasugamycin')

def read_dbn(path):
    lines=[l.strip() for l in path.read_text().splitlines() if l.strip()]
    assert len(lines)%3==0
    out={}
    for i in range(0,len(lines),3):
        name,seq,structure=lines[i][1:],lines[i+1],lines[i+2]
        assert len(seq)==len(structure)
        stack=[];pairs=[];partners={}
        for p,c in enumerate(structure):
            if c=='(':stack.append(p)
            elif c==')':
                q=stack.pop();pairs.append((q,p));partners[q]=p;partners[p]=q
            else:assert c=='.'
        assert not stack and all(seq[a]+seq[b] in PAIRS for a,b in pairs)
        out[name]={'sequence':seq,'structure':structure,'pairs':sorted(pairs),'partners':partners}
    return out

def read_react(path):
    out={};name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):name=line[1:];out[name]=[]
        else:
            i,v=line.split();assert int(i)==len(out[name])+1
            out[name].append(None if v=='NA' else float(v))
    return out

def select(data_dir):
    data_dir=Path(data_dir)
    records={c:read_dbn(data_dir/'weeks194_orig'/f'{c}.dbn') for c in CONDITIONS}
    reactivities={c:read_react(data_dir/'weeks194_orig'/f'{c}.react') for c in CONDITIONS}
    rows=list(csv.DictReader((data_dir/'cds_mapping.csv').open()))
    overlaps=defaultdict(set)
    with (data_dir/'coding_position_map.csv').open() as handle:
        for row in csv.DictReader(handle):
            if row['overlaps_other_CDS']=='1':overlaps[row['transcript_id']].add(int(row['transcript_pos'])-1)
    selected=[];audit=[]
    for row in rows:
        tid=row['transcript_id'];seq=records['incell'][tid]['sequence']
        assert all(records[c][tid]['sequence']==seq for c in CONDITIONS)
        start,end=int(row['transcript_start'])-1,int(row['transcript_end'])
        cds=seq[start:end];counts=Counter();candidates=[]
        if '*' in translate(cds[3:-3]):
            audit.append({**row,'eligible_windows':0,'reason':'internal_stop_or_recoding'});continue
        for pos in range(start+3,end-3-14,3):
            variable=seq[pos:pos+15];pep=translate(variable)
            n=math.prod(len(CODONS[a]) for a in pep)
            if n>1024:counts['cartesian_cap']+=1;continue
            if any(p in overlaps[tid] for p in range(pos,pos+15)):
                counts['overlapping_CDS']+=1;continue
            if any(sum(reactivities[c][tid][p] is not None and reactivities[c][tid][p]>=0 for p in range(pos,pos+15))<12 for c in CONDITIONS[:2]):
                counts['coverage']+=1;continue
            cs=pos-67;ce=cs+150
            if cs<0 or ce>len(seq):counts['context_boundary']+=1;continue
            pp={};bad=False
            for c in CONDITIONS:
                partners=records[c][tid]['partners']
                pairs=sorted(set((min(p,partners[p]),max(p,partners[p])) for p in range(pos,pos+15) if p in partners))
                if any(a<cs or b>=ce or b-a<4 for a,b in pairs):bad=True;break
                if c in CONDITIONS[:2] and len(pairs)<3:bad=True;break
                pp[c]=[(a-cs,b-cs) for a,b in pairs]
            if bad:counts['pair_context_or_count']+=1;continue
            key=f"{tid}:{row['locus_tag']}:{pos}:1729"
            candidates.append({'id':hashlib.sha256(key.encode()).hexdigest()[:16],
                'transcript_id':tid,'gene':row['gene'],'locus_tag':row['locus_tag'],
                'reference':'U00096.2','strand':row['strand'],'transcript_pos0':pos,
                'codon_index1':(pos-start)//3+1,'peptide':pep,'native':variable,
                'context_start0':cs,'variable_start0':pos-cs,'context':seq[cs:ce],
                'touching_pairs':pp,'cartesian_size':n,
                'coverage':{c:sum(reactivities[c][tid][p] is not None and reactivities[c][tid][p]>=0 for p in range(pos,pos+15))/15 for c in CONDITIONS},
                'reactivity':{c:reactivities[c][tid][cs:ce] for c in CONDITIONS}})
        if candidates:selected.append(min(candidates,key=lambda r:r['id']))
        audit.append({**row,'eligible_windows':len(candidates),'reason':'eligible' if candidates else 'no_eligible_window','sequential_rejections':dict(counts)})
    selected.sort(key=lambda r:r['id'])
    (ROOT/'data/native_windows.json').write_text(json.dumps(selected,indent=2,allow_nan=False))
    (ROOT/'results/native_selection.json').write_text(json.dumps({'selected':len(selected),'transcripts':len(set(r['transcript_id'] for r in selected)),'audit':audit},indent=2))
    print(json.dumps({'selected':len(selected),'transcripts':len(set(r['transcript_id'] for r in selected)),'cartesian_total':sum(r['cartesian_size'] for r in selected)}))

def evaluate_one(task):
    import RNA
    block,condition=task;started=time.time();seq=block['native'];vstart=block['variable_start0'];context=block['context'];pairs=block['touching_pairs'][condition]
    internal=[(a-vstart,b-vstart) for a,b in pairs if vstart<=a<b<vstart+len(seq)]
    allseq=enumerate_space(block['peptide'],cap=1024)
    compatible_seq=[];fullseq=[]
    for s in allseq:
        full=context[:vstart]+s+context[vstart+len(s):]
        if compatible(full,pairs):compatible_seq.append(s);fullseq.append(full)
    graph=adjacency(compatible_seq);labels,sizes=components(graph)
    sig=Counter(signature(s,internal) for s in compatible_seq)
    assert sorted(sig.values())==sorted(sizes)
    occupancy=[];mfe=[];energy=[]
    md=RNA.md();md.temperature=37;md.dangles=2;md.noLP=0
    for full in fullseq:
        fc=RNA.fold_compound(full,md)
        structure,e=fc.mfe();fc.exp_params_rescale(e);fc.pf();bpp=fc.bpp()
        occupancy.append(sum(bpp[a+1][b+1] for a,b in pairs)/len(pairs) if pairs else None)
        mfe.append(structure);energy.append(e)
    thresholds={}
    for cutoff in (0.1,0.25,0.5,0.75):
        keep=[i for i,value in enumerate(occupancy) if value is not None and value>=cutoff]
        thresholds[str(cutoff)]=excess_fragmentation(compatible_seq,internal,graph,keep)
    result={'id':block['id'],'transcript_id':block['transcript_id'],'gene':block['gene'],'locus_tag':block['locus_tag'],
        'condition':condition,'context_nt':len(context),'n_pairs':len(pairs),'n_internal_pairs':len(internal),
        'cartesian_size':len(allseq),'n_compatible':len(compatible_seq),'component_sizes':sorted(sizes,reverse=True),
        'a_compatible':connected_probability(sizes),'thresholds':thresholds,
        'genotypes':compatible_seq,'edges':[(i,j) for i,neighbors in enumerate(graph) for j in neighbors if i<j],
        'pair_retention':occupancy,'mfe_structures':mfe,'mfe_energy_kcal':energy,
        'native_retention':occupancy[compatible_seq.index(seq)],'seconds':time.time()-started,'engine':RNA.__version__}
    return result

def run(workers,limit):
    blocks=json.loads((ROOT/'data/native_windows.json').read_text())
    if limit:blocks=blocks[:limit]
    outdir=ROOT/'results/native';outdir.mkdir(exist_ok=True)
    tasks=[(r,c) for r in blocks for c in CONDITIONS if not (outdir/f"{r['id']}_{c}.json").exists()]
    print('pending',len(tasks),flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(evaluate_one,t):t for t in tasks}
        for i,f in enumerate(as_completed(futures),1):
            r=f.result();(outdir/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False))
            if i%10==0:print('completed',i,'of',len(tasks),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['select','run']);parser.add_argument('--data-dir');parser.add_argument('--workers',type=int,default=8);parser.add_argument('--limit',type=int,default=0);args=parser.parse_args()
    if args.action=='select':select(args.data_dir)
    else:run(args.workers,args.limit)
