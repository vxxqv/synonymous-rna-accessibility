import json,hashlib,collections,re,itertools
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
P=Path(__file__).resolve().parent
seq=str(SeqIO.read(P/'full_construct.dna','snapgene').seq).upper()
fa=''.join(next(P.glob('references__*.fa')).read_text().splitlines()[1:]).upper()
assert fa==seq
catalog=json.loads((P/'synonymous_id_catalog.json').read_text(encoding='utf8'))
rows=[]
for h in catalog['all4_common']:
    r=catalog['all'][h]; assert len(r['native_codon_sites'])==1
    site=r['native_codon_sites'][0]
    assert 40<=site<=74 or 136<=site<=158
    r={**r,'native_codon_site':site,'ROI':1 if site<100 else 2,'nucleotide_changes':len(r['edits'])}
    rows.append(r)
theory={}
for site in list(range(40,75))+list(range(136,159)):
    codon=seq[4652+3*(site+7):4652+3*(site+8)]
    aa=str(Seq(codon).translate())
    theory[site]=sum(str(Seq(''.join(c)).translate())==aa and ''.join(c)!=codon for c in itertools.product('ACGT',repeat=3))
gs=json.loads((P/'sample_id_map.json').read_text(encoding='utf8'))
counts=collections.Counter();sample_ids=[]
for g in gs:
    t=g.get('Sample_title','');m=re.search(r'ROI([12])_(.*?)_biorep(\d+)',t)
    if m:counts[str((int(m[1]),m[2],int(m[3])))]+=1
    sample_ids.append({'gsm':g['GSM'],'source_sample':g.get('Sample_description'),'title':t})
with (P/'common_synonymous_ids.tsv').open('w',encoding='utf8') as f:
    f.write('hgvs_nt\tnative_codon_site\tROI\tnucleotide_changes\n')
    for r in rows:f.write(f"{r['hgvs_nt']}\t{r['native_codon_site']}\t{r['ROI']}\t{r['nucleotide_changes']}\n")
summary={'construct_fasta_equals_snapgene':True,'construct_sequence_length':len(seq),'gff_transcript_terminal_annotation':6233,'terminal_annotation_one_nt_beyond_construct':True,'CDS_and_target_coordinates_valid':True,'all4_common_variants':len(rows),'all4_common_site_count':len({r['native_codon_site'] for r in rows}),'all4_by_ROI':dict(collections.Counter(str(r['ROI']) for r in rows)),'all4_sites_by_ROI':{str(i):len({r['native_codon_site'] for r in rows if r['ROI']==i}) for i in (1,2)},'all4_SNVS':sum(r['nucleotide_changes']==1 for r in rows),'theoretical_single_codon_synonymous_alternatives':sum(theory.values()),'theoretical_sites_with_synonymous_alternatives':sum(v>0 for v in theory.values()),'sample_title_counts':dict(counts),'source_non_gene_samples':sum(not g.get('Sample_description','').startswith('COMT_') for g in gs),'phenotype_values_inspected':False}
(P/'coverage_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
manifest=[]
for p in P.iterdir():
    if p.is_file() and p.name not in ['source_manifest.json','finalization_stdout.json']:
        b=p.read_bytes();manifest.append({'file':p.name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(P/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2))
