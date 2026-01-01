from pathlib import Path
from Bio import SeqIO
import json,csv,collections,itertools
root=(Path(__file__).resolve().parents[2]/'data/pairmap');data=root/'source';gb=SeqIO.read(root/'U00096.2.gb','genbank');genome=str(gb.seq);comp=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(comp)[::-1]
def readct(p):
 lines=p.read_text().splitlines();n=int(lines[0].split()[0]);rows=[x.split() for x in lines[1:n+1]];assert len(rows)==n
 sequence=''.join(r[1] for r in rows).upper().replace('T','U');pair=[];table=[int(r[4]) for r in rows]
 for k,r in enumerate(rows,1):
  assert int(r[0])==k
  if table[k-1]:assert table[table[k-1]-1]==k
  if table[k-1]>k:pair.append((k,table[k-1]))
 return dict(sequence=sequence,pairs=pair,partners=table,length=n,crossings=sum(a<c<b<d for (a,b),(c,d) in itertools.combinations(pair,2)))
summary=[];allpairs=[];mapping=[];ct_records={}
for gene in ['rpsM','rpsB']:
 fa=SeqIO.read(data/'alignment_sequences'/f'{gene}.fa','fasta');s=str(fa.seq).upper();hits=[]
 for strand,ss in [(1,s),(-1,rc(s))]:
  p=genome.find(ss)
  while p>=0:hits.append((strand,p+1,p+len(s)));p=genome.find(ss,p+1)
 assert len(hits)==1,hits;st,a,b=hits[0];feat=next(f for f in gb.features if f.type=='CDS' and f.qualifiers.get('gene')==[gene]);cdspos=list(feat.location);posmap={pos:i for i,pos in enumerate(cdspos)}
 for i in range(1,len(s)+1):
  gp=a+i-1 if st==1 else b-i+1;idx=posmap.get(gp-1)
  mapping.append(dict(gene=gene,reference='U00096.2',dataset_pos=i,genomic_pos=gp,strand='+' if st==1 else '-',is_CDS=idx is not None,codon_index=idx//3+1 if idx is not None else '',codon_phase=idx%3+1 if idx is not None else '',base=s[i-1],protein_id=feat.qualifiers.get('protein_id',[''])[0],locus_tag=feat.qualifiers.get('locus_tag',[''])[0]))
 info=dict(gene=gene,length=len(s),reference='U00096.2',strand=st,start=a,end=b,CDS_dataset_positions=[r['dataset_pos'] for r in mapping if r['gene']==gene and r['is_CDS']],source_full_sequence_exact=True,models=[])
 for p in list((data/'reference_structures').glob(f'{gene}*.ct'))+list((data/'processed_data'/'structure_models').rglob(f'*{gene}*.ct')):
  rec=readct(p);assert rec['sequence']==s.replace('T','U');key=gene+'/'+str(p.relative_to(data)).replace('\\','/');ct_records[key]=rec
  row=dict(file=str(p.relative_to(data)).replace('\\','/'),n_pairs=len(rec['pairs']),crossings=rec['crossings']);info['models'].append(row)
  for i,j in rec['pairs']:allpairs.append(dict(gene=gene,model_file=str(p.relative_to(data)),i=i,j=j,i_in_CDS=(a+i-2 if st==1 else b-i) in posmap,j_in_CDS=(a+j-2 if st==1 else b-j) in posmap))
 summary.append(info)
for fn,rows in [('coding_position_map.csv',mapping),('structure_pairs.csv',allpairs)]:
 with (root/fn).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(root/'mapping_summary.json').write_text(json.dumps(summary,indent=2));(root/'ct_models.json').write_text(json.dumps(ct_records,indent=2));print([{k:v for k,v in s.items() if k!='CDS_dataset_positions'}|dict(CDS_range=[min(s['CDS_dataset_positions']),max(s['CDS_dataset_positions'])]) for s in summary])
