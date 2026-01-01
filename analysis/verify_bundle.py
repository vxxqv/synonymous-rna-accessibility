"""Cross-check archived native graphs, control fingerprints and exact receipts."""
from pathlib import Path
import sys,json,hashlib,subprocess,math,argparse,os
ROOT=Path(__file__).resolve().parents[1]
# Deep, byte-pinned provenance paths can exceed Windows' ordinary path limit.
# Keep every tracked file in the audit instead of treating those paths as absent.
if os.name=='nt' and not str(ROOT).startswith('\\\\?\\'):
    ROOT=Path('\\\\?\\UNC\\'+str(ROOT)[2:] if str(ROOT).startswith('\\\\')
              else '\\\\?\\'+str(ROOT))
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'analysis/expansion'))
from accessibility import translate,compatible,adjacency,components,signature,connected_probability,excess_fragmentation
from provenance import signature as cache_signature,cache_valid

def require(value,message):
    if not value:raise ValueError(message)

def verify():
    blocks={b['id']:b for b in json.loads((ROOT/'data/native_windows.json').read_text())};graphs=0;genotypes=0
    for path in sorted((ROOT/'results/native').glob('*.json')):
        r=json.loads(path.read_text());b=blocks[r['id']];seq=r['genotypes'];v=b['variable_start0'];pairs=b['touching_pairs'][r['condition']];internal=[(a-v,d-v) for a,d in pairs if v<=a<d<v+len(b['native'])]
        require(r['n_compatible']==len(seq),str(path)+': cardinality');require(len(set(seq))==len(seq),str(path)+': duplicate genotype')
        for s in seq:
            require(translate(s)==b['peptide'],str(path)+': translation');require(compatible(b['context'][:v]+s+b['context'][v+len(s):],pairs),str(path)+': pair compatibility')
        graph=adjacency(seq);edges=[(i,j) for i,neighbors in enumerate(graph) for j in neighbors if i<j];require(edges==[tuple(e) for e in r['edges']],str(path)+': edge completeness')
        _,sizes=components(graph);require(sorted(sizes,reverse=True)==r['component_sizes'],str(path)+': component sizes');require(abs(connected_probability(sizes)-r['a_compatible'])<1e-12,str(path)+': connection probability')
        require(r['native_retention']==r['pair_retention'][seq.index(b['native'])],str(path)+': native score')
        for t in (.1,.25,.5,.75):
            calculated=excess_fragmentation(seq,internal,graph,[i for i,value in enumerate(r['pair_retention']) if value>=t]);old=r['thresholds'][str(t)]
            require(calculated==old,str(path)+': filtered statistic')
        graphs+=1;genotypes+=len(seq)
    require(graphs==1212,'expected complete frozen cohort')
    families=[('conditional_null','conditional_vertex_null.py',{'draws':999,'cap':2000,'seed':1729,'cutoff':.5}),('filtration','retention_filtration.py',{'draws':199,'seed':1729}),('smoothness','smoothness_controls.py',{'draws':199,'fit_draws':64,'cap':40320,'seed':1729}),('additive','additive_score_control.py',{'tie_draws':99,'seed':1729})]
    cached=0
    for folder,method,parameters in families:
        for b in blocks.values():
            for c in ('incell','cellfree'):
                path=ROOT/'results/expansion'/folder/f"{b['id']}_{c}.json";require(cache_valid(path,cache_signature(method,b,c,parameters)),str(path)+': stale cache');cached+=1
    receipts=[]
    for name in ('verification_release_v2','feasibility_verification_release_v2'):
        manifest=json.loads((ROOT/'results/expansion'/name/'manifest.json').read_text())
        for item in manifest['input_artifacts']:
            require(hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'],name+': input hash changed '+item['path'])
        for item in manifest['outputs']:
            require(hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'],name+': output hash changed '+item['path'])
        receipts.append(name)
    # Newline handling must not change any tracked byte-level source receipt.
    tracked=subprocess.check_output(['git','ls-tree','-rz','HEAD'],cwd=ROOT).split(b'\0');different=[]
    for item in filter(None,tracked):
        header,name=item.split(b'\t',1);name=name.decode();expected=header.split()[2].decode();path=ROOT/name
        if not path.is_file():different.append(name);continue
        raw=path.read_bytes();actual=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if actual!=expected:different.append(name)
    require(not different,'tracked working files differ from HEAD: '+', '.join(different))
    return {'status':'passed','native_complete_graphs':graphs,'native_genotype_instances':genotypes,'current_control_fingerprints':cached,'verified_exact_receipts':receipts,'tracked_worktree_differences':different,'interpretation':'Integrity and computational checks, not proof-assistant certification or prospective biological validation'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/bundle_verification.json');args=p.parse_args();result=verify();path=ROOT/args.out;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
