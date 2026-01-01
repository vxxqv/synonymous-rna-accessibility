"""Metadata gate. Never decodes phenotype, count, confidence or FDR values."""
import csv,hashlib,json,re,tarfile,collections,itertools
from pathlib import Path
from Bio.Seq import Seq
P=Path(__file__).resolve().parent
def selected_csv(path,selected):
    # Only allowlisted identifier bytes are retained or decoded. Quotes handled explicitly.
    with path.open('rb') as f:
        header=f.readline().rstrip(b'\r\n').decode('utf8').split(',')
        cols={header.index(k):k for k in selected}
        for line in f:
            if not line.strip():continue
            result={};i=0;quoted=False;buf=bytearray();pos=0
            while pos<len(line):
                c=line[pos]
                if c==34:
                    if quoted and pos+1<len(line) and line[pos+1]==34:
                        if i in cols:buf.append(c)
                        pos+=2;continue
                    quoted=not quoted
                elif c==44 and not quoted:
                    if i in cols:result[cols[i]]=bytes(buf).decode('utf8')
                    i+=1;buf.clear()
                elif c not in (10,13) and i in cols:buf.append(c)
                pos+=1
            if quoted:raise ValueError('Multiline fields require separate audit')
            if i in cols:result[cols[i]]=bytes(buf).decode('utf8')
            yield result
seq=''.join(next(P.glob('references__*.fa')).read_text().splitlines()[1:]).upper()
target=json.loads((P/'mavedb_urn_mavedb_00000661-c-1.json').read_text(encoding='utf8'))['targetGenes'][0]['targetSequence']['sequence']
assert len(seq)==6232 and seq[4793:5150]==target
assert all(json.loads(p.read_text(encoding='utf8'))['targetGenes'][0]['targetSequence']['sequence']==target for p in P.glob('mavedb_urn*.json'))
reference_protein=str(Seq(seq[4652:6203]).translate())
assert len(reference_protein)==517 and '*' not in reference_protein and seq[6203:6206] in ('TAA','TAG','TGA')
allsets={};summaries={};unparsed=[];catalog={}
for assay in ['rna','mid','heavy','protein']:
    rows=list(selected_csv(P/f'mavedb_{assay}_scores.csv',['accession','hgvs_nt','hgvs_pro','hgvs_splice']))
    ids=set();syn=set();types=collections.Counter();sites=set()
    for row in rows:
        h=row['hgvs_nt'];ids.add(h)
        edits=[]
        for atom in h.removeprefix('c.').strip('[]').split(';'):
            m=re.fullmatch(r'(\d+)([ACGT])>([ACGT])',atom)
            if m:edits.append((int(m[1]),m[2],m[3]))
            else:
                m=re.fullmatch(r'(\d+)_(\d+)delins([ACGT]+)',atom)
                if not m or int(m[2])-int(m[1])+1!=len(m[3]):break
                edits.extend((pos,target[pos-1],alt) for pos,alt in zip(range(int(m[1]),int(m[2])+1),m[3]) if target[pos-1]!=alt)
        else:
            mutated=list(target)
            for pos,a,b in edits:
                assert target[pos-1]==a,(assay,h,pos)
                mutated[pos-1]=b
            mutated=''.join(mutated)
            full=seq[:4793]+mutated+seq[5150:]
            same=str(Seq(full[4652:6203]).translate())==reference_protein
            site=sorted({39+(x[0]-1)//3+1 for x in edits})
            # target begins native codon40; 39+floor((p-1)/3)+1.
            types[(len(edits),'synonymous' if same else 'other')]+=1
            if same:
                syn.add(h);sites.update(site)
                catalog[h]={'hgvs_nt':h,'native_codon_sites':site,'edits':edits,'full_translation_identical':True,'hgvs_pro':row['hgvs_pro']}
            continue
        unparsed.append({'assay':assay,**row})
    assert len(rows)==len(ids),'Duplicate nucleotide IDs'
    allsets[assay]={'all':ids,'syn':syn}
    summaries[assay]={'rows':len(rows),'synonymous_rows':len(syn),'synonymous_sites':len(sites),'types':{str(k):v for k,v in types.items()}}
common=set.intersection(*(x['syn'] for x in allsets.values()))
common_rp=allsets['rna']['syn']&allsets['protein']['syn']
(P/'synonymous_id_catalog.json').write_text(json.dumps({'all':catalog,'all4_common':sorted(common),'RNA_protein_common':sorted(common_rp)},indent=2),encoding='utf8')
# Header-only archive inspection. No data values are read in this pass.
with tarfile.open(P/'author_archive.tar.gz','r:gz') as t:
    members=[m for m in t.getmembers() if m.name.endswith('.var.cand.vcf.summary.txt')]
    headers={m.name:t.extractfile(m).readline().decode('utf8').rstrip().split('\t') for m in members}
(P/'raw_headers.json').write_text(json.dumps(headers,indent=2),encoding='utf8')
geo=[];entry=None
for line in (P/'geo_metadata.txt').read_text(encoding='utf8').splitlines():
    if line.startswith('^SAMPLE = '):entry={'GSM':line.split(' = ',1)[1]};geo.append(entry)
    elif entry is not None and line.startswith(('!Sample_title = ','!Sample_description = ')):
        a,b=line.split(' = ',1);entry[a[1:]]=b
(P/'sample_id_map.json').write_text(json.dumps(geo,indent=2),encoding='utf8')
out={'reference_length_nt':len(seq),'fused_CDS_length_without_stop':1551,'fusion_protein_length':len(reference_protein),'target_length_nt':len(target),'target_native_codons':[40,158],'mutagenesis_native_codons':[[40,74],[136,158]],'all_target_references_identical':True,'assays':summaries,'all4_common_synonymous':len(common),'all4_common_sites':len({s for x in common for s in catalog[x]['native_codon_sites']}),'rna_protein_common_synonymous':len(common_rp),'rna_protein_common_sites':len({s for x in common_rp for s in catalog[x]['native_codon_sites']}),'unparsed_identifiers':unparsed,'raw_variant_summary_members':len(headers),'GEO_samples':len(geo),'score_values_decoded':0,'uncertainty_values_decoded':0,'count_values_decoded':0}
(P/'sequence_schema_audit.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps({**out,'unparsed_identifiers':len(unparsed)},indent=2))
