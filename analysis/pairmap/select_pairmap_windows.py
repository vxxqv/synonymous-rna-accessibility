import sys,json,csv,math,itertools,collections
from pathlib import Path
root=(Path(__file__).resolve().parents[2]/'data/pairmap');sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'));from accessibility import *
models=json.loads((root/'ct_models.json').read_text());maps=list(csv.DictReader((root/'coding_position_map.csv').open()));out=[];audit=[]
for gene in ['rpsM','rpsB']:
 rows=[r for r in maps if r['gene']==gene and r['is_CDS']=='True'];start=min(int(r['dataset_pos']) for r in rows)-1;end=max(int(r['dataset_pos']) for r in rows);allprofiles={}
 for c in ['incell','cellfree']:
  p=root/'source'/'processed_data'/'shapemapper_profiles'/f'{c}_{gene}_profile.txt';allprofiles[c]=list(csv.DictReader(p.open(),delimiter='\t'))
 seq=models[f'{gene}/processed_data/structure_models/pairmap-dms/incell_{gene}.all.f.ct']['sequence'];assert all(len(v)==len(seq) for v in allprofiles.values())
 for pos in range(start+3,end-14,3):
  pep=translate(seq[pos:pos+15]);n=math.prod(len(CODONS[a]) for a in pep);reason=[]
  if n>1024:reason.append('cartesian_cap')
  coverage={c:sum(math.isfinite(float(allprofiles[c][p]['HQ_profile'])) for p in range(pos,pos+15))/15 for c in allprofiles}
  if min(coverage.values())<0.8:reason.append('HQ_coverage')
  pp={}
  for c in allprofiles:
   rec=models[f'{gene}/processed_data/structure_models/pairmap-dms/{c}_{gene}.all.f.ct'];pp[c]=[(a-1,b-1) for a,b in rec['pairs'] if pos<=a-1<pos+15 or pos<=b-1<pos+15]
   if len(pp[c])<3:reason.append(c+'_fewer3touchingpairs')
  audit.append(dict(gene=gene,pos1=pos+1,peptide=pep,cartesian_size=n,reasons=reason,coverage=coverage))
  if reason:continue
  out.append(dict(id=f'PAIRMAP_{gene}_{pos+1}',transcript_id=f'PAIRMAP2019_{gene}',gene=gene,locus_tag=rows[0]['locus_tag'],reference='U00096.2',strand=rows[0]['strand'],transcript_pos0=pos,codon_index1=(pos-start)//3+1,peptide=pep,native=seq[pos:pos+15],context_start0=0,variable_start0=pos,context=seq,touching_pairs=pp,cartesian_size=n,coverage=coverage,source='original_Figshare2019_CT_pairmap+DMS',reactivity={c:[None if not math.isfinite(float(r['HQ_profile'])) else float(r['HQ_profile']) for r in allprofiles[c]] for c in allprofiles}))
(root/'eligible_windows.json').write_text(json.dumps(out,indent=2));(root/'window_selection_audit.json').write_text(json.dumps(audit,indent=2));print('eligible',len(out),collections.Counter(r['gene'] for r in out),'cartesian total',sum(r['cartesian_size'] for r in out));print([(r['id'],r['cartesian_size']) for r in out])
