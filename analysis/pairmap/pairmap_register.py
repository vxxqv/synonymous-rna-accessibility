from pathlib import Path
import json,csv,collections,itertools
root=(Path(__file__).resolve().parents[2]/'data/pairmap');models=json.loads((root/'ct_models.json').read_text());mapping=list(csv.DictReader((root/'coding_position_map.csv').open()));allrows=[];runs=[]
for gene in ['rpsM','rpsB']:
 cds={int(r['dataset_pos']) for r in mapping if r['gene']==gene and r['is_CDS']=='True'}
 for condition in ['incell','cellfree']:
  p=root/'source'/'processed_data'/'pairmap_correlations'/f'{condition}_{gene}-pairmap.txt';lines=p.read_text().splitlines();rs=list(csv.DictReader(lines[1:],delimiter='\t'))
  for r in rs:
   i,j=int(r['i']),int(r['j']);r.update(gene=gene,condition=condition,first_window_CDS_n=sum(i+k in cds for k in range(3)),second_window_CDS_n=sum(j+k in cds for k in range(3)));allrows.append(r)
 a=models[f'{gene}/processed_data/structure_models/pairmap-dms/cellfree_{gene}.all.f.ct'];b=models[f'{gene}/processed_data/structure_models/pairmap-dms/incell_{gene}.all.f.ct'];pa=dict(a['pairs']);pb=dict(b['pairs']);N=a['length']
 for delta in [1,2,3,-1,-2,-3]:
  positions=sorted(i for i,j in pa.items() if pb.get(i)==j+delta)
  groups=[]
  for i in positions:
   if groups and i==groups[-1][-1]+1 and pa[i]==pa[groups[-1][-1]]-1:groups[-1].append(i)
   else:groups.append([i])
  for run in groups:
   if len(run)>=3:
    pairs1=[(i,pa[i]) for i in run];pairs2=[(i,pb[i]) for i in run];arms_cds=all(i in cds and j in cds for i,j in pairs1+pairs2)
    runs.append(dict(gene=gene,delta=delta,length=len(run),left_start=run[0],left_end=run[-1],pairs_cellfree=pairs1,pairs_incell=pairs2,both_arms_CDS=arms_cds))
print('register candidate runs',runs)
print('evidence summaries',collections.Counter((r['gene'],r['condition'],r['Class'],bool(r['first_window_CDS_n'] or r['second_window_CDS_n'])) for r in allrows))
with (root/'pairmap_windows.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
(root/'register_shift_audit.json').write_text(json.dumps(dict(definition='maximal >=3bp antiparallel runs sharing left positions and shifted right partner by constant +/-1..3 between deposited pairmap+DMS condition CTs; no outcome selection',runs=runs,eligible_supported_coding_register_examples=0 if not any(r['both_arms_CDS'] for r in runs) else 'requires_window_support_audit'),indent=2))
