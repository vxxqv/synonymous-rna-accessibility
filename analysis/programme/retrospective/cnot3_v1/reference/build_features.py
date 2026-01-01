"""Outcome-blind, pinned reference proxy. Standard code; no folds or motifs."""
import sys
if not __debug__: raise RuntimeError('Optimization is refused before input access')
sys.dont_write_bytecode = True
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re

CODON_ORDER = tuple(''.join(x) for x in itertools.product('ACGT', repeat=3))
# NCBI standard code in T,C,A,G nested order, used independently of installed tools.
_BASES = 'TCAG'
_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = dict(zip((''.join(x) for x in itertools.product(_BASES,repeat=3)),_AA))
SENSE = tuple(c for c in CODON_ORDER if CODE[c] != '*')
AMINO = tuple(sorted(set(CODE.values())-{'*'}))
CHROMS = {'chr'+str(i) for i in range(1,23)} | {'chrX','chrY'}
PINNED = {'gtf': '37a298c2b57e41988cfe12667134106523c69c24a1807347679d19e88aa1f4a9',
          'fasta': 'd29486ce127cf4ab0f681fdcfa9bcd9fe2cfd18e7331b75257eb12c4cbb385f5',
          'ids': '85fb60552f3c12477d28dbeeafbd171d52828ecb66e30591cefa60f77f482c52'}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk:=stream.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def text_sha(text): return hashlib.sha256(text.encode('ascii')).hexdigest()

def stable(gene):
    match=re.fullmatch(r'(ENSG\d+)\.\d+',gene)
    return match[1] if match else gene

def attrs(text):
    result=defaultdict(list)
    for key,value in re.findall(r'(\w+) "([^"]*)";',text): result[key].append(value)
    return dict(result)

def merge(intervals):
    result=[]
    for start,end in sorted(intervals):
        if start>=end: raise ValueError('Nonpositive interval')
        if result and start<=result[-1][1]: result[-1]=(result[-1][0],max(end,result[-1][1]))
        else: result.append((start,end))
    return result

def subtract(intervals,cut):
    result=[]
    for a,b in intervals:
        if a<min(b,cut[0]): result.append((a,min(b,cut[0])))
        if max(a,cut[1])<b: result.append((max(a,cut[1]),b))
    return merge(result)

def project(exons,strand,regions):
    """1-based closed genomic intervals -> 0-based half-open spliced intervals."""
    ordered=sorted(exons,reverse=(strand=='-'))
    if strand not in ('+','-'): raise ValueError('Invalid strand')
    genomic=sorted(exons)
    if any(a[1]>=b[0] for a,b in zip(genomic,genomic[1:])): raise ValueError('Overlapping exons')
    answer=[]
    for left,right in regions:
        matches=[]; offset=0
        for a,b in ordered:
            if a<=left<=right<=b:
                matches.append((offset+left-a,offset+right-a+1) if strand=='+' else (offset+b-right,offset+b-left+1))
            offset+=b-a+1
        if len(matches)!=1: raise ValueError('Feature is not contained in exactly one exon')
        answer.extend(matches)
    return merge(answer)

def regions(header):
    fields=header.split('|'); result={}
    if len(fields)<7: raise ValueError('Invalid GENCODE FASTA header')
    for field in fields[7:]:
        if not field: continue
        match=re.fullmatch(r'(UTR5|CDS|UTR3):(\d+)-(\d+)',field)
        if match:
            if match[1] in result: raise ValueError('Duplicate FASTA region')
            result[match[1]]=(int(match[2])-1,int(match[3]))
    return fields,result

def numeric(cds,utr5,utr3):
    if len(cds)<6 or len(cds)%3: raise ValueError('cds_length')
    if set(cds)-set('ACGT'): raise ValueError('ambiguous_cds')
    codons=[cds[i:i+3] for i in range(0,len(cds),3)]
    if codons[0]!='ATG': raise ValueError('non_ATG_start')
    if CODE[codons[-1]]!='*': raise ValueError('missing_terminal_standard_stop')
    if any(CODE[c]=='*' for c in codons[:-1]): raise ValueError('internal_standard_stop')
    sense=codons[:-1]; counts=Counter(sense); aa=''.join(CODE[c] for c in sense); aa_counts=Counter(aa)
    result={**{'codon_'+c:counts[c]/len(sense) for c in SENSE},
            **{'aa_'+a:aa_counts[a]/len(aa) for a in AMINO},
            'CDS_length_nt':len(cds),'CDS_GC':(cds.count('G')+cds.count('C'))/len(cds)}
    for name,sequence in [('UTR5',utr5),('UTR3',utr3)]:
        if set(sequence)-set('ACGT'): raise ValueError('ambiguous_'+name)
        result[name+'_length_nt']=len(sequence)
        result[name+'_GC']=(sequence.count('G')+sequence.count('C'))/len(sequence) if sequence else None
    return result,aa

def validate(record,header,sequence):
    fields,declared=regions(header)
    if fields[0]!=record['transcript_id'] or fields[1]!=record['gene_id']: raise ValueError('header_identity_mismatch')
    if int(fields[6])!=len(sequence): raise ValueError('header_length_mismatch')
    tags=set(record['tags'])
    if tags & {'cds_start_NF','cds_end_NF'}: raise ValueError('annotated_incomplete_CDS')
    if record['features'].get('Selenocysteine'): raise ValueError('selenocysteine_standard_code_incompatible')
    exons=record['features'].get('exon',[])
    if not exons or sum(b-a+1 for a,b in exons)!=len(sequence): raise ValueError('exon_transcript_length_mismatch')
    projected={name:project(exons,record['strand'],record['features'].get(name,[])) for name in ('CDS','start_codon','stop_codon','UTR')}
    full=merge(projected['CDS']+projected['stop_codon'])
    if len(full)!=1: raise ValueError('noncontiguous_spliced_CDS')
    begin,end=full[0]
    if projected['start_codon']!=[(begin,begin+3)] or projected['stop_codon']!=[(end-3,end)]: raise ValueError('incomplete_start_stop_annotation')
    if declared.get('CDS')!=(begin,end): raise ValueError('GTF_FASTA_CDS_coordinate_mismatch')
    expected_utr=[x for x in [(0,begin),(end,len(sequence))] if x[0]<x[1]]
    if subtract(projected['UTR'],(begin,end))!=expected_utr: raise ValueError('GTF_UTR_coordinate_mismatch')
    for label,expected in [('UTR5',(0,begin)),('UTR3',(end,len(sequence)))]:
        if expected[0]<expected[1] and declared.get(label)!=expected: raise ValueError('GTF_FASTA_'+label+'_coordinate_mismatch')
        if expected[0]==expected[1] and label in declared: raise ValueError('Unexpected_zero_length_'+label)
    cds,utr5,utr3=sequence[begin:end],sequence[:begin],sequence[end:]
    # Validate all bases and translate now. Fractions are produced only for selected rows.
    result,protein=numeric(cds,utr5,utr3)
    return result,protein,{'CDS_spliced_0based_halfopen':[begin,end],
                          'UTR5_spliced_0based_halfopen':[0,begin],
                          'UTR3_spliced_0based_halfopen':[end,len(sequence)],
                          'transcript_sha256':text_sha(sequence),'CDS_sha256':text_sha(cds),
                          'protein_sha256':text_sha(protein),
                          'UTR5_sha256':text_sha(utr5),'UTR3_sha256':text_sha(utr3)}

def load_gtf(path,targets):
    genes={}; transcripts={}; counts=Counter()
    with gzip.open(path,'rt',encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#'): continue
            columns=line.rstrip('\n').split('\t')
            if len(columns)!=9: raise ValueError('GTF column count')
            match=re.search(r'gene_id "([^"]+)";',columns[8])
            if not match: raise ValueError('GTF gene_id missing')
            gene=stable(match[1])
            if gene not in targets: continue
            a=attrs(columns[8]); kind=columns[2]; chromosome=columns[0]
            if kind=='gene':
                if gene in genes: raise ValueError('Stable gene maps to multiple GTF gene records')
                genes[gene]={'gene_id':a['gene_id'][0],'gene_name':a['gene_name'][0],
                             'gene_type':a['gene_type'][0],'chromosome':chromosome,
                             'strand':columns[6],'hgnc_id':a.get('hgnc_id',[None])[0]}
            if a['gene_type'][0]!='protein_coding' or chromosome not in CHROMS: continue
            if 'transcript_id' not in a: continue
            tid=a['transcript_id'][0]
            if kind=='transcript':
                if tid in transcripts: raise ValueError('Duplicate GTF transcript')
                transcripts[tid]={'gene_id':a['gene_id'][0],'stable_gene_id':gene,'transcript_id':tid,
                                  'transcript_type':a['transcript_type'][0],'chromosome':chromosome,
                                  'strand':columns[6],'tags':a.get('tag',[]),'features':defaultdict(list)}
            elif kind in ('exon','CDS','start_codon','stop_codon','UTR','Selenocysteine'):
                if tid not in transcripts: raise ValueError('GTF transcript does not precede features')
                r=transcripts[tid]
                if r['chromosome']!=chromosome or r['strand']!=columns[6]: raise ValueError('Transcript location mismatch')
                r['features'][kind].append((int(columns[3]),int(columns[4])))
    return genes,transcripts

def fasta(path):
    header=None; chunks=[]
    with gzip.open(path,'rt',encoding='ascii') as stream:
        for line in stream:
            if line.startswith('>'):
                if header is not None: yield header,''.join(chunks)
                header=line[1:].strip(); chunks=[]
            else: chunks.append(line.strip())
        if header is not None: yield header,''.join(chunks)

def run(args):
    paths={name:Path(getattr(args,name)).resolve() for name in PINNED}
    before={name:sha(path) for name,path in paths.items()}
    if before!=PINNED: raise ValueError('Pinned reference/GeneID-only input hashes differ')
    out=Path(args.out).resolve()
    if out.exists(): raise FileExistsError(out)
    ids=json.loads(paths['ids'].read_text(encoding='utf-8'))
    targets={row['id'] for row in ids['entries']}
    if len(targets)!=len(ids['entries']) or any(not re.fullmatch(r'ENSG\d+',x) for x in targets): raise ValueError('Unexpected ID inventory')
    out.mkdir(parents=True)
    genes,records=load_gtf(paths['gtf'],targets)
    print(json.dumps({'GTF_target_genes':len(genes),'candidate_transcripts':len(records)}),flush=True)
    best={}; seen=set(); valid=Counter(); invalid=Counter(); per_gene_invalid=defaultdict(Counter); lengths=defaultdict(list)
    with gzip.open(out/'candidate_audit.jsonl.gz','wt',encoding='utf-8',newline='\n') as audit:
        for header,sequence in fasta(paths['fasta']):
            tid=header.split('|',1)[0]
            if tid not in records: continue
            if tid in seen: raise ValueError('Duplicate FASTA target transcript')
            seen.add(tid); r=records[tid]; gene=r['stable_gene_id']
            record={'gene_id':gene,'transcript_id':tid,'transcript_type':r['transcript_type']}
            try:
                features,protein,coordinates=validate(r,header,sequence)
            except ValueError as exc:
                reason=str(exc); invalid[reason]+=1; per_gene_invalid[gene][reason]+=1
                record.update({'status':'excluded','reason':reason})
            else:
                valid[gene]+=1; lengths[gene].append((features['CDS_length_nt'],tid))
                record.update({'status':'valid','CDS_length_nt':features['CDS_length_nt']})
                key=(-features['CDS_length_nt'],tid)
                if gene not in best or key<best[gene]['selection_key']:
                    metadata={**genes[gene],**coordinates,'stable_gene_id':gene,'transcript_id':tid,
                              'transcript_type':r['transcript_type'],'tags':r['tags'],
                              'sense_codon_count':len(protein),
                              'UTR5_status': 'annotated_incomplete' if 'mRNA_start_NF' in r['tags'] else ('annotated_nonempty' if features['UTR5_length_nt'] else 'no_annotated_bases'),
                              'UTR3_status': 'annotated_incomplete' if 'mRNA_end_NF' in r['tags'] else ('annotated_nonempty' if features['UTR3_length_nt'] else 'no_annotated_bases')}
                    best[gene]={'selection_key':key,'features':features,'metadata':metadata}
            audit.write(json.dumps(record,separators=(',',':'))+'\n')
    missing=Counter()
    for tid,r in records.items():
        if tid not in seen and r['features'].get('CDS'):
            missing[r['stable_gene_id']]+=1
    # Source omission is explicit; it is not silently treated as an invalid sequence.
    numeric_names=['codon_'+x for x in SENSE]+['aa_'+x for x in AMINO]+['CDS_length_nt','CDS_GC','UTR5_length_nt','UTR5_GC','UTR3_length_nt','UTR3_GC']
    with gzip.open(out/'features.csv.gz','wt',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['GeneID','transcript_id']+numeric_names)
        writer.writeheader()
        for gene in sorted(best): writer.writerow({'GeneID':gene,'transcript_id':best[gene]['metadata']['transcript_id'],**best[gene]['features']})
    with gzip.open(out/'selected_transcripts.jsonl.gz','wt',encoding='utf-8',newline='\n') as stream:
        for gene in sorted(best):
            selected=best[gene]; selected['metadata']['valid_candidate_count']=valid[gene]
            longest=-selected['selection_key'][0]
            selected['metadata']['longest_CDS_tie_transcripts']=sorted(tid for length,tid in lengths[gene] if length==longest)
            stream.write(json.dumps(selected['metadata'],separators=(',',':'))+'\n')
    status_counts=Counter()
    with (out/'gene_inventory.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['GeneID','status','gene_type','chromosome','selected_transcript','valid_candidates','missing_coding_FASTA_transcripts','invalid_reason_counts'])
        writer.writeheader()
        for gene in sorted(targets):
            g=genes.get(gene,{})
            status=('selected' if gene in best else 'not_in_CHR_annotation' if not g else 'non_nuclear_reference' if g['chromosome'] not in CHROMS else 'not_protein_coding_gene' if g['gene_type']!='protein_coding' else 'no_valid_complete_standard_CDS')
            status_counts[status]+=1
            writer.writerow({'GeneID':gene,'status':status,'gene_type':g.get('gene_type',''),'chromosome':g.get('chromosome',''),
                             'selected_transcript':best[gene]['metadata']['transcript_id'] if gene in best else '',
                             'valid_candidates':valid[gene],'missing_coding_FASTA_transcripts':missing[gene],
                             'invalid_reason_counts':json.dumps(per_gene_invalid[gene],sort_keys=True)})
    grouping=defaultdict(list)
    for gene,row in best.items(): grouping[row['metadata']['protein_sha256']].append(gene)
    groups=[{'protein_sha256':h,'GeneIDs':sorted(genes)} for h,genes in sorted(grouping.items()) if len(genes)>1]
    (out/'exact_protein_group_option.json').write_text(json.dumps({'partition_applied':False,'grouping':'Identical complete standard-code translated representative proteins only','nontrivial_groups':groups,'limitation':'Broader homologs and curated families are not guaranteed to share a group; distinct genes are not independent biological units.'},indent=2)+'\n',encoding='utf-8')
    after={name:sha(path) for name,path in paths.items()}
    if after!=before: raise RuntimeError('Input mutation')
    summary={'status':'reference_proxy_features_complete','target_gene_count':len(targets),'annotation_target_genes':len(genes),
             'candidate_transcripts':len(records),'candidate_FASTA_records_seen':len(seen),'coding_records_missing_FASTA':sum(missing.values()),
             'valid_transcripts':sum(valid.values()),'selected_genes':len(best),'gene_status_counts':dict(status_counts),
             'transcript_exclusions':dict(invalid),'features':numeric_names,'numeric_feature_count':len(numeric_names),
             'longest_CDS_tied_genes':sum(len(r['metadata']['longest_CDS_tie_transcripts'])>1 for r in best.values()),
             'selected_transcript_types':dict(Counter(r['metadata']['transcript_type'] for r in best.values())),
             'UTR5_status':dict(Counter(r['metadata']['UTR5_status'] for r in best.values())),
             'UTR3_status':dict(Counter(r['metadata']['UTR3_status'] for r in best.values())),
             'exact_protein_nontrivial_groups':len(groups),'exact_protein_grouped_genes':sum(len(g['GeneIDs']) for g in groups),
             'input_sha256_before':before,'input_sha256_after':after,'inputs_immutable':before==after,
             'code_sha256':sha(Path(__file__)),'outcome_fields_read':[],'folds':0,'partitions_applied':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('status','selected_genes','gene_status_counts','coding_records_missing_FASTA','numeric_feature_count')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('gtf','fasta','ids','out'): p.add_argument('--'+name,required=True)
    run(p.parse_args())
