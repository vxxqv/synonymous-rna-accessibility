"""Independent folding-engine sensitivity on complete, hash-selected domains.

Requires an externally installed RNAstructure command-line distribution.
The distribution itself and transient PFS files are not part of the archive.
"""
from pathlib import Path
import sys,json,os,argparse,subprocess,tempfile,hashlib,time,math
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).parent))
from accessibility import components,excess_fragmentation
from retention_filtration import filtration

def run_block(task):
    block,install,scratch=task;install=Path(install);scratch=Path(scratch)
    env=os.environ.copy();env['DATAPATH']=str(install/'data_tables');env['OMP_NUM_THREADS']='1'
    extension='.exe' if os.name=='nt' else ''
    records={c:json.loads((ROOT/'results/native'/f"{block['id']}_{c}.json").read_text()) for c in ('incell','cellfree')}
    sequences=sorted(set(s for r in records.values() for s in r['genotypes']))
    pairs=sorted(set(tuple(p) for c in records for p in block['touching_pairs'][c]));v=block['variable_start0'];context=block['context']
    cache={};start=time.time()
    with tempfile.TemporaryDirectory(dir=scratch,prefix='rs_') as tmp:
        tmp=Path(tmp);fa=tmp/'s.fa';pfs=tmp/'s.pfs';dp=tmp/'p.txt'
        for seq in sequences:
            full=context[:v]+seq+context[v+len(seq):];fa.write_text('>sequence\n'+full+'\n')
            cmd=[str(install/'exe'/('partition'+extension)),str(fa),str(pfs),'-q','-i','-t','310.15']
            r=subprocess.run(cmd,env=env,capture_output=True,text=True,timeout=120)
            if r.returncode or not pfs.exists():raise RuntimeError(f'Partition failed: {r.stderr[:500]} {r.stdout[:500]}')
            r=subprocess.run([str(install/'exe'/('ProbabilityPlot'+extension)),str(pfs),str(dp),'-t','-min','0','-max','1000'],env=env,capture_output=True,text=True,timeout=120)
            if r.returncode or not dp.exists():raise RuntimeError(f'Probability export failed: {r.stderr[:500]}')
            lines=dp.read_text().splitlines()
            if int(lines[0])!=len(full) or 'log10' not in lines[1]:raise ValueError('Unexpected probability format')
            probs={};bounds={}
            for line in lines[2:]:
                a,b,z=line.split();a,b=int(a)-1,int(b)-1
                if (a,b) not in pairs:continue
                z=float(z);p=10**(-z);unit=10**(math.floor(math.log10(abs(z)))-5) if z else 0.
                probs[a,b]=p;bounds[a,b]=(10**(-z-.5*unit),min(1.,10**(-z+.5*unit)))
            cache[seq]={'mean':{c:sum(probs.get(tuple(p),0.) for p in block['touching_pairs'][c])/len(block['touching_pairs'][c]) for c in records},
               'bounds':{c:[sum(bounds.get(tuple(p),(0.,0.))[j] for p in block['touching_pairs'][c])/len(block['touching_pairs'][c]) for j in (0,1)] for c in records}}
            pfs.unlink();dp.unlink()
    outputs=[]
    for c,r in records.items():
        vals=[cache[s]['mean'][c] for s in r['genotypes']];bounds=[cache[s]['bounds'][c] for s in r['genotypes']]
        graph=[[] for _ in vals]
        for a,b in r['edges']:graph[a].append(b);graph[b].append(a)
        labels,_=components(graph);internal=[(a-v,b-v) for a,b in block['touching_pairs'][c] if v<=a<b<v+len(block['native'])]
        thresholds={}
        for tau in (.1,.25,.5,.75):thresholds[str(tau)]=excess_fragmentation(r['genotypes'],internal,graph,[i for i,x in enumerate(vals) if x>=tau])
        outputs.append({'id':block['id'],'condition':c,'transcript_id':block['transcript_id'],'gene':block['gene'],'locus_tag':block['locus_tag'],
          'engine':'RNAstructure 6.6','temperature_K':310.15,'isolated_pairs':True,'coaxial_stacking':'default enabled',
          'genotypes':r['genotypes'],'pair_retention':vals,'retention_rounding_bounds':bounds,'thresholds':thresholds,
          'cutoff_ambiguities':{str(t):sum(lo<t<=hi for lo,hi in bounds) for t in (.1,.25,.5,.75)},
          'filtration':filtration(graph,vals,labels,True),'vienna_filtration':json.loads((ROOT/'results/expansion/filtration'/f"{block['id']}_{c}.json").read_text())['observed'],
          'seconds_shared_block':time.time()-start,'unique_context_genotypes_evaluated':len(sequences)})
    return outputs

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--install',required=True);p.add_argument('--scratch',required=True);p.add_argument('--limit',type=int,default=24);p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    blocks=json.loads((ROOT/'data/native_windows.json').read_text())[:a.limit];out=ROOT/'results/expansion/rnastructure';out.mkdir(parents=True,exist_ok=True);Path(a.scratch).mkdir(exist_ok=True,parents=True)
    tasks=[(b,str(Path(a.install).resolve()),str(Path(a.scratch).resolve())) for b in blocks if not all((out/f"{b['id']}_{c}.json").exists() for c in ('incell','cellfree'))]
    print('pending blocks',len(tasks),flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(run_block,t) for t in tasks]):
            rr=f.result()
            for r in rr:(out/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False))
            print('completed block',rr[0]['id'],'genotypes',rr[0]['unique_context_genotypes_evaluated'],flush=True)
    rows=[json.loads(p.read_text()) for p in out.glob('*.json') if p.name!='summary.json']
    s={'selection':'first '+str(a.limit)+' frozen SHA-ordered blocks; no energetic-outcome selection','graphs':len(rows),'primary_nonempty':sum(r['thresholds']['0.5']['n']>0 for r in rows),
       'primary_extra':sum((r['thresholds']['0.5']['excess'] or 0)>1e-12 for r in rows),'any_threshold_extra':sum(r['filtration']['positive_bar_count']>0 for r in rows),
       'longest_lifetime_mean':sum(r['filtration']['longest_lifetime'] for r in rows)/len(rows),'integrated_excess_mean':sum(r['filtration']['integrated_excess'] for r in rows)/len(rows),
       'cutoff_ambiguities':sum(sum(r['cutoff_ambiguities'].values()) for r in rows)}
    (out/'summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
