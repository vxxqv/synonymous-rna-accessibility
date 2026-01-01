"""Independent cohort audit. Never imports author code or decodes KO targets."""
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import time
from collections import Counter, defaultdict

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASE=ROOT/'work/decay_cycle4'
RUN=BASE/'execution/operational_cohort_v2'
CODE=dict(zip((''.join(x) for x in itertools.product('TCAG',repeat=3)),
              'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
CODONS=sorted(k for k,v in CODE.items() if v!='*')
AMINO=sorted(set(CODE.values())-{'*'})
FEATURES={'codon_'+c for c in CODONS}|{'aa_'+a for a in AMINO}|{'CDS_length_nt','CDS_GC','UTR5_length_nt','UTR5_GC','UTR3_length_nt','UTR3_GC'}
FIELDS=FEATURES|{'GeneID','group_ID','HL_WT','LOG2_FC_IP','SE_IP'}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def require(condition,message):
    if not condition:raise ValueError(message)
def finite(value):return type(value) in (int,float) and math.isfinite(value)


def selective_rows(text):
    """Parse flat row structure; skip the KO literal without numeric decoding."""
    decoder=json.JSONDecoder();pos=0;rows=[];skipped=0
    literal=re.compile(r'"(?:[^"\\]|\\.)*"|-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?|true|false|null')
    def ws():
        nonlocal pos
        while pos<len(text) and text[pos].isspace():pos+=1
    def char(expected):
        nonlocal pos
        ws();require(text[pos:pos+1]==expected,'Malformed source JSON structure');pos+=1
    char('[')
    while True:
        ws()
        if text[pos:pos+1]==']':pos+=1;break
        char('{');row={};keys=set()
        while True:
            ws();key,end=decoder.raw_decode(text,pos);pos=end
            require(isinstance(key,str) and key not in keys,'Duplicate/malformed source key');keys.add(key)
            char(':');ws();match=literal.match(text,pos);require(match is not None,'Nonflat source cell')
            if key=='HL_CNOT3KO':skipped+=1
            else:row[key]=json.loads(match.group())
            pos=match.end();ws()
            if text[pos:pos+1]=='}':pos+=1;break
            char(',')
        require(keys=={'GeneID','HL_WT','HL_CNOT3KO','LOG2_FC_IP','SE_IP'},'Source row schema')
        rows.append(row);ws()
        if text[pos:pos+1]==']':continue
        char(',')
    ws();require(pos==len(text),'Trailing source content')
    return rows,skipped


def distance(a,b):
    """Elementary full Wagner-Fischer recurrence; separate from native alignment."""
    if len(a)>len(b):a,b=b,a
    previous=list(range(len(a)+1))
    for i,right in enumerate(b,1):
        current=[i]
        for j,left in enumerate(a,1):
            current.append(min(current[-1]+1,previous[j]+1,previous[j-1]+(left!=right)))
        previous=current
    return previous[-1]


def run():
    started=time.monotonic()
    # Poisoned target literal proves this reader does not decode that field.
    poison='[{"GeneID":"SYN","HL_WT":1,"HL_CNOT3KO":"UNREAD_TARGET","LOG2_FC_IP":0,"SE_IP":1}]'
    toy,skipped=selective_rows(poison);require(skipped==1 and 'HL_CNOT3KO' not in toy[0],'Selective reader boundary')
    for a,b,d in [('', '',0),('A','',1),('AC','CA',2),('AAAA','AA',2),('ABC','AXC',1),('ACGT','ACGTT',1)]:
        require(distance(a,b)==d,'Independent edit-distance fixture')
    receipt=load(RUN/'receipt.json')
    require(receipt['status']=='cohort_group_partition_ready','Cohort phase incomplete')
    sources={'rows':BASE/'execution/reproduction_v1/cohort_rows.json',
             'original_ids':BASE/'data_eligibility/original_gene_ids_only.json',
             'features':BASE/'reference_features/reference_v1/features.csv.gz',
             'metadata':BASE/'reference_features/reference_v1/selected_transcripts.jsonl.gz',
             'fasta':BASE/'reference_features/gencode.v48.pc_transcripts.fa.gz',
             'protocol':BASE/'operational_successor/protocol.md',
             'original_gate':BASE/'execution/principal_gate.json'}
    checked_inputs={key:digest(path) for key,path in sources.items()}
    require(receipt['input_sha256_before']==receipt['input_sha256_after'],'Input mutation receipt')
    for key,value in checked_inputs.items():require(value==receipt['input_sha256_before'][key],'Input hash '+key)
    selected_outputs=['eligibility.json','grouping.json','assignments.json','design_guards.json',
                      'model_inputs/training_predictors.json','model_inputs/holdout_predictors.json']
    checked_outputs={name:digest(RUN/name) for name in selected_outputs}
    for name,value in checked_outputs.items():require(value==receipt['output_sha256'][name],'Output hash '+name)
    rows,skipped=selective_rows(sources['rows'].read_text(encoding='utf-8'))
    row_by_id={r['GeneID']:r for r in rows};require(len(row_by_id)==len(rows)==5480,'Source IDs')
    original_ids=[r['id'] for r in load(sources['original_ids'])['entries']]
    source_ids={}
    for raw in original_ids:
        m=re.fullmatch(r'(ENSG[0-9]+)(?:\.[0-9]+)?',raw)
        if m:
            require(m[1] not in source_ids,'Normalized source ID collision');source_ids[m[1]]=raw
    meta={}
    with gzip.open(sources['metadata'],'rt',encoding='utf-8') as f:
        for line in f:
            record=json.loads(line);gene=record['stable_gene_id'];require(gene not in meta,'Metadata duplicate');meta[gene]=record
    features={}
    with gzip.open(sources['features'],'rt',encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f);require(set(reader.fieldnames)==FEATURES|{'GeneID','transcript_id'},'Feature header')
        for row in reader:
            gene=row['GeneID'];require(gene not in features,'Feature duplicate')
            require(gene in meta and row['transcript_id']==meta[gene]['transcript_id'],'Reference transcript join')
            features[gene]={k:None if row[k]=='' else float(row[k]) for k in FEATURES}
    eligibility=load(RUN/'eligibility.json');expected_excluded={};kept=[]
    for gene,row in sorted(row_by_id.items()):
        reasons=[]
        if gene not in source_ids:reasons.append('absent_original_SLAM_ID_inventory')
        if gene not in features or gene not in meta:reasons.append('absent_selected_reference')
        else:
            f,m=features[gene],meta[gene];tags=set(m['tags'])
            if m['transcript_type']!='protein_coding':reasons.append('selected_transcript_not_protein_coding')
            if tags&{'mRNA_start_NF','mRNA_end_NF'}:reasons.append('incomplete_UTR_annotation')
            if tags&{'cds_start_NF','cds_end_NF'}:reasons.append('incomplete_CDS_annotation')
            n=f['CDS_length_nt'];sense=m['sense_codon_count']
            if not finite(n) or n<6 or n%3 or type(sense)!=int or sense<1 or n!=3*(sense+1):reasons.append('invalid_complete_CDS_length')
            invalid=False
            for key,value in f.items():
                if key in ('UTR5_GC','UTR3_GC') and f[key[:4]+'_length_nt']==0:
                    invalid|=value is not None
                elif not finite(value):invalid=True
                elif key.endswith('_length_nt'):invalid|=value<0 or int(value)!=value
                else:invalid|=not 0<=value<=1
            if invalid:reasons.append('nonfinite_or_invalid_reference_feature')
            else:
                inconsistent=(abs(math.fsum(f['codon_'+c] for c in CODONS)-1)>1e-12 or abs(math.fsum(f['aa_'+a] for a in AMINO)-1)>1e-12)
                inconsistent|=any(abs(math.fsum(f['codon_'+c] for c in CODONS if CODE[c]==a)-f['aa_'+a])>1e-12 for a in AMINO)
                if inconsistent:reasons.append('reference_composition_inconsistent')
        if not finite(row['LOG2_FC_IP']):reasons.append('nonfinite_enrichment')
        if not finite(row['SE_IP']) or row['SE_IP']<=0:reasons.append('nonpositive_or_nonfinite_SE_IP')
        if reasons:expected_excluded[gene]=reasons
        else:kept.append(gene)
    require(kept==eligibility['eligible_GeneIDs'],'Eligible gene set mismatch')
    actual_excluded={r['GeneID']:r['all_reasons'] for r in eligibility['excluded']}
    require(expected_excluded==actual_excluded,'Exclusion reasons mismatch')
    S=len(set(row_by_id)&set(source_ids));A=len(kept)
    require((eligibility['P'],eligibility['S'],eligibility['A'])==(5480,S,A),'Coverage count mismatch')
    require(eligibility['A_over_S']==A/S and eligibility['S_over_P']==S/5480,'Coverage arithmetic')
    first=Counter(v[0] for v in expected_excluded.values());remaining=5480
    for stage in eligibility['sequential_counts']:
        require(stage['excluded_at_this_stage']==first[stage['reason']],'Sequential exclusion count')
        remaining-=first[stage['reason']];require(stage['remaining']==remaining,'Sequential remaining count')
    assignments=load(RUN/'assignments.json');by_id={r['GeneID']:r for r in assignments}
    require(set(by_id)==set(kept) and len(by_id)==len(assignments),'Assignment gene universe')
    grouping=load(RUN/'grouping.json');require(grouping['status']=='complete','Incomplete final grouping')
    component_by_id={};counts=Counter();folds=Counter();groups=defaultdict(set)
    for component in grouping['components']:
        members=component['GeneIDs'];group=component['group_ID']
        require(members==sorted(set(members)) and group==min(members),'Canonical component')
        outer=hashlib.sha256(('cnot3-focus-v1|outer|'+group).encode()).digest()
        test=int.from_bytes(outer,'big')%5==0;kind='holdout' if test else 'training'
        inner=None if test else hashlib.sha256(('cnot3-focus-v1|inner|'+group).encode()).digest()
        for gene in members:
            require(gene not in component_by_id,'Overlapping components');component_by_id[gene]=group
            a=by_id[gene]
            expected=dict(GeneID=gene,group_ID=group,outer=kind,inner_fold=None if test else int.from_bytes(inner,'big')%5,
                          outer_SHA256=outer.hex(),inner_SHA256=None if test else inner.hex())
            require(a==expected,'Hash partition mismatch '+gene);counts[kind]+=1;groups[kind].add(group)
            if not test:folds[a['inner_fold']]+=1
    require(set(component_by_id)==set(kept),'Component universe')
    for kind in ('training','holdout'):
        predictors=load(RUN/'model_inputs'/f'{kind}_predictors.json')
        ids=[r['GeneID'] for r in predictors]
        require(ids==sorted(g for g in kept if by_id[g]['outer']==kind),'Predictor partition membership')
        for row in predictors:
            require(set(row)==FIELDS and len(row)==92,'Predictor schema or outcome leakage')
            gene=row['GeneID'];expected={'GeneID':gene,'group_ID':component_by_id[gene],**{k:row_by_id[gene][k] for k in ['HL_WT','LOG2_FC_IP','SE_IP']},**features[gene]}
            require(row==expected,'Predictor/reference/source value mismatch '+gene)
    print('Independent joins, all predictor fields and partitions passed.',flush=True)
    # Reconstruct only eligible proteins from the pinned transcript FASTA.
    chosen={meta[g]['transcript_id']:g for g in kept};proteins={}
    def consume(header,chunks):
        if header is None:return
        tid=header.split('|')[0]
        if tid not in chosen:return
        gene=chosen[tid];require(gene not in proteins,'Repeated selected transcript')
        seq=''.join(chunks);m=meta[gene];lo,hi=m['CDS_spliced_0based_halfopen'];cds=seq[lo:hi]
        require(len(cds)%3==0 and cds.startswith('ATG'),'Complete CDS')
        triplets=[cds[i:i+3] for i in range(0,len(cds),3)]
        require(CODE[triplets[-1]]=='*' and all(CODE[c]!='*' for c in triplets[:-1]),'Stop convention')
        protein=''.join(CODE[c] for c in triplets[:-1]);proteins[gene]=protein
        for label,string in [('transcript',seq),('CDS',cds),('UTR5',seq[:lo]),('UTR3',seq[hi:]),('protein',protein)]:
            require(hashlib.sha256(string.encode('ascii')).hexdigest()==m[label+'_sha256'],'Reference sequence hash '+gene)
    with gzip.open(sources['fasta'],'rt',encoding='ascii') as f:
        header=None;chunks=[]
        for line in f:
            if line.startswith('>'):consume(header,chunks);header=line[1:].strip();chunks=[]
            elif header is not None and header.split('|')[0] in chosen:chunks.append(line.strip())
        consume(header,chunks)
    require(set(proteins)==set(kept),'Protein selected universe')
    identical=defaultdict(list)
    for gene in kept:identical[proteins[gene]].append(gene)
    identity=[{'representative':members[0],'members':members,'protein_sha256':hashlib.sha256(seq.encode()).hexdigest()}
              for seq,members in identical.items() if len(members)>1]
    require(sorted(identity,key=lambda x:x['representative'])==grouping['identity_groups'],'Identity groups')
    entries=sorted((members[0],seq) for seq,members in identical.items())
    kmers=[Counter(seq[i:i+5] for i in range(len(seq)-4)) for _,seq in entries]
    posting=defaultdict(list)
    for index,counts5 in enumerate(kmers):
        for word,count in counts5.items():posting[word].append((index,count))
    expected_postings=sum(len(x)*(len(x)-1)//2 for x in posting.values())
    require(expected_postings==grouping['counts']['posting_pair_events'],'Posting count mismatch')
    candidates=[];edges=[]
    for i,(gene,left) in enumerate(entries):
        shared=[0]*len(entries)
        for word,count in kmers[i].items():
            for j,other_count in posting[word]:
                if j>i:shared[j]+=min(count,other_count)
        for j,common in enumerate(shared):
            if j<=i:continue
            other,right=entries[j];maximum=max(len(left),len(right));budget=maximum//10
            if budget==0 or abs(len(left)-len(right))>budget:continue
            # For maximum>=10, maximum-4-5*floor(maximum/10)>0.
            # Therefore every passing pair has a common5-mer; no missing postings.
            require(maximum-4-5*budget>0,'Lossless candidate positivity')
            if common<maximum-4-5*budget:continue
            d=distance(left,right);candidates.append(dict(left=gene,right=other,distance=d))
            if 10*d<=maximum:edges.append({'left_representative':gene,'right_representative':other,'distance':d,'left_length':len(left),'right_length':len(right)})
    require(edges==grouping['near_protein_edges'],'Exact near edge set/distance mismatch')
    require(len(candidates)==grouping['counts']['exact_pairs_completed'],'Candidate completeness count')
    cells=sum(len(proteins[x['left']])*len(proteins[x['right']]) for x in candidates)
    require(cells==grouping['counts']['alignment_cells_completed'],'Exact pair cell count')
    adjacency={g:set() for g in kept}
    for members in identical.values():
        for gene in members[1:]:adjacency[members[0]].add(gene);adjacency[gene].add(members[0])
    for edge in edges:
        a,b=edge['left_representative'],edge['right_representative'];adjacency[a].add(b);adjacency[b].add(a)
    visited=set();recovered=[]
    for gene in kept:
        if gene in visited:continue
        pending=[gene];members=[];visited.add(gene)
        while pending:
            node=pending.pop();members.append(node)
            for neighbor in adjacency[node]-visited:visited.add(neighbor);pending.append(neighbor)
        recovered.append({'group_ID':min(members),'GeneIDs':sorted(members)})
    require(recovered==grouping['components'],'Connected components mismatch')
    guards=load(RUN/'design_guards.json')
    require(guards['status']=='design_coverage_passed' and all(guards['checks'].values()),'Stored design guards')
    require(A*5>=4*S and counts['training']>=2000 and counts['holdout']>=300,'Independent gene coverage')
    require(len(groups['training'])>=100 and len(groups['holdout'])>=30,'Independent group coverage')
    for fold in range(5):
        fitting_groups={r['group_ID'] for r in assignments if r['outer']=='training' and r['inner_fold']!=fold}
        require(folds[fold]>=100 and len(fitting_groups)>=50,'Independent inner fold coverage')
    report=dict(status='passed',independent_logic='stdlib only; no author module imports',source_sha256=checked_inputs,
                output_sha256=checked_outputs,source_KO_literals_skipped=skipped,source_KO_numeric_decodes=0,
                target_files_opened=0,model_fits=0,P=5480,S=S,A=A,reference_coverage=A/S,
                exclusions_by_first_reason=dict(first),training_genes=counts['training'],holdout_genes=counts['holdout'],
                training_groups=len(groups['training']),holdout_groups=len(groups['holdout']),inner_assessment_genes=dict(folds),
                predictor_fields=92,protein_sequence_hashes_checked=5*A,unique_proteins=len(entries),
                posting_pair_events=expected_postings,all_candidate_pairs_independently_checked=len(candidates),
                exact_DP_cells=cells,near_edges=len(edges),identity_groups=len(identity),components=len(recovered),
                candidate_scope='Complete lossless5-mer candidate set independently regenerated; all candidate global edit distances independently recomputed',
                limits='No biological-replicate, broader paralog/domain independence, or outcome accuracy claim; paired-HL validity relies on separate source audit',
                elapsed_seconds=time.monotonic()-started,review_source_sha256=digest(__file__))
    (HERE/'result.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (HERE/'independent_candidate_distances.json').write_text(json.dumps(candidates,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    try:run()
    except Exception as error:
        (HERE/'failure.json').write_text(json.dumps({'type':type(error).__name__,'message':str(error)},indent=2)+'\n',encoding='utf-8')
        raise
