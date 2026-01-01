"""Frozen-rule cohort primitives. No XLSX reader, fitting, tuning or outcome search."""
import sys
if not __debug__: raise RuntimeError('Python -O/-OO is refused before input access')
sys.dont_write_bytecode=True
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import hashlib
import itertools
import math
import re
import time
from Bio.Align import PairwiseAligner

BASES='TCAG'
AA64='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE=dict(zip((''.join(x) for x in itertools.product(BASES,repeat=3)),AA64))
SENSE=tuple(sorted(c for c,a in CODE.items() if a!='*'))
AMINO=tuple(sorted(set(CODE.values())-{'*'}))
FEATURE_NAMES=tuple('codon_'+c for c in SENSE)+tuple('aa_'+a for a in AMINO)+('CDS_length_nt','CDS_GC','UTR5_length_nt','UTR5_GC','UTR3_length_nt','UTR3_GC')
ROW_FIELDS={'GeneID','HL_WT','HL_CNOT3KO','LOG2_FC_IP','SE_IP'}

class IntegrityError(ValueError): pass

class ResourceLimit(RuntimeError):
    def __init__(self,cap,receipt):
        super().__init__('Grouping resource cap: '+cap)
        self.receipt=receipt

@dataclass(frozen=True)
class Limits:
    posting_pair_events:int=20_000_000
    candidate_pairs:int=50_000
    alignment_cells:int=3_000_000_000
    elapsed_seconds:float=300.0

LIMITS=Limits()

def stable_original(value):
    if not isinstance(value,str): raise IntegrityError('Original ID must be a string')
    match=re.fullmatch(r'(ENSG\d+)(?:\.\d+)?',value)
    return match[1] if match else None

def source_membership(values):
    """Strict original-ID normalization. No aliases or collision resolution."""
    sources={}; excluded=[]; raw_seen=set()
    for value in values:
        if not isinstance(value,str) or value in raw_seen: raise IntegrityError('Malformed/duplicate original source ID')
        raw_seen.add(value); stable=stable_original(value)
        if stable is None:
            excluded.append({'source_id':value,'reason':'not_ordinary_ENSG_or_terminal_version'})
            continue
        if stable in sources: raise IntegrityError('Distinct source IDs collide after terminal version stripping: '+stable)
        sources[stable]=value
    return frozenset(sources),{'mapped_source_IDs':dict(sorted(sources.items())),'unmappable_source_IDs':sorted(excluded,key=lambda r:r['source_id'])}

def finite_number(value):
    if isinstance(value,bool) or not isinstance(value,(int,float)): return False
    try: return math.isfinite(value)
    except (ValueError,OverflowError): return False

def require_reference_schema(features,meta,gene):
    if not isinstance(features,dict) or set(features)!=set(FEATURE_NAMES): raise IntegrityError('Reference feature schema mismatch: '+gene)
    needed={'stable_gene_id','transcript_id','transcript_type','tags','sense_codon_count','CDS_sha256','protein_sha256','transcript_sha256','UTR5_sha256','UTR3_sha256','CDS_spliced_0based_halfopen'}
    if not isinstance(meta,dict) or not needed.issubset(meta): raise IntegrityError('Reference metadata schema mismatch: '+gene)
    if meta['stable_gene_id']!=gene or not isinstance(meta['tags'],list) or any(not isinstance(x,str) for x in meta['tags']): raise IntegrityError('Reference metadata identity/tags mismatch')
    if not isinstance(meta['transcript_id'],str) or not re.fullmatch(r'ENST\d+\.\d+',meta['transcript_id']): raise IntegrityError('Expected versioned transcript ID')
    for key in ('CDS_sha256','protein_sha256','transcript_sha256','UTR5_sha256','UTR3_sha256'):
        if not isinstance(meta[key],str) or not re.fullmatch('[0-9a-f]{64}',meta[key]): raise IntegrityError('Malformed reference sequence hash')

def reference_reasons(features,meta):
    reasons=[]
    tags=set(meta['tags'])
    if meta['transcript_type']!='protein_coding': reasons.append('selected_transcript_not_protein_coding')
    if tags & {'mRNA_start_NF','mRNA_end_NF'}: reasons.append('incomplete_UTR_annotation')
    if tags & {'cds_start_NF','cds_end_NF'}: reasons.append('incomplete_CDS_annotation')
    bad=False
    for name,value in features.items():
        if name in ('UTR5_GC','UTR3_GC') and features[name[:4]+'_length_nt']==0:
            if value is not None: bad=True
            continue
        if not finite_number(value): bad=True; continue
        if name.endswith('_length_nt'):
            if int(value)!=value or value<0: bad=True
        elif not 0<=value<=1: bad=True
    length=features['CDS_length_nt']; sense=meta['sense_codon_count']
    if not finite_number(length) or length<6 or length%3 or not isinstance(sense,int) or isinstance(sense,bool) or sense<1 or length!=3*(sense+1):
        reasons.append('invalid_complete_CDS_length')
    if bad: reasons.append('nonfinite_or_invalid_reference_feature')
    if not bad:
        if abs(math.fsum(features['codon_'+c] for c in SENSE)-1)>1e-12 or abs(math.fsum(features['aa_'+a] for a in AMINO)-1)>1e-12: reasons.append('reference_composition_inconsistent')
        if any(abs(math.fsum(features['codon_'+c] for c in SENSE if CODE[c]==a)-features['aa_'+a])>1e-12 for a in AMINO): reasons.append('reference_composition_inconsistent')
    return list(dict.fromkeys(reasons))

def eligibility(rows,original_ids,reference_features,reference_metadata):
    """Pure join and declared exclusions; no grouping and no replacement transcript."""
    source_set,source_report=source_membership(original_ids)
    seen=set(); normalized=[]
    for row in rows:
        if not isinstance(row,dict) or set(row)!=ROW_FIELDS: raise IntegrityError('Malformed reproduction row schema')
        gene=row['GeneID']
        if not isinstance(gene,str) or not re.fullmatch(r'ENSG\d+',gene) or gene in seen: raise IntegrityError('Malformed/duplicate reproduction GeneID')
        seen.add(gene)
        for name in ('HL_WT','HL_CNOT3KO'):
            if not finite_number(row[name]) or row[name]<=0: raise IntegrityError('Reproduction cohort must retain only positive finite half-lives')
        for name in ('LOG2_FC_IP','SE_IP'):
            if row[name] is not None and (isinstance(row[name],bool) or not isinstance(row[name],(int,float))): raise IntegrityError('Malformed numeric reproduction field')
        normalized.append(dict(row))
    kept=[]; excluded=[]; first_reason=Counter(); denominator=0
    for row in sorted(normalized,key=lambda r:r['GeneID']):
        gene=row['GeneID']; reasons=[]
        if gene not in source_set: reasons.append('absent_original_SLAM_ID_inventory')
        else: denominator+=1
        if gene not in reference_features or gene not in reference_metadata:
            reasons.append('absent_selected_reference')
        else:
            require_reference_schema(reference_features[gene],reference_metadata[gene],gene)
            reasons.extend(reference_reasons(reference_features[gene],reference_metadata[gene]))
        if not finite_number(row['LOG2_FC_IP']): reasons.append('nonfinite_enrichment')
        if not finite_number(row['SE_IP']) or row['SE_IP']<=0: reasons.append('nonpositive_or_nonfinite_SE_IP')
        if reasons:
            excluded.append({'GeneID':gene,'primary_reason':reasons[0],'all_reasons':reasons})
            first_reason[reasons[0]]+=1
        else: kept.append(row)
    # coverage denominator is Figure2 intersect original source inventory, before reference exclusions.
    coverage=len(kept)/denominator if denominator else None
    sequential=[]; remaining=len(normalized)
    reason_order=('absent_original_SLAM_ID_inventory','absent_selected_reference','selected_transcript_not_protein_coding','incomplete_UTR_annotation','incomplete_CDS_annotation','invalid_complete_CDS_length','nonfinite_or_invalid_reference_feature','reference_composition_inconsistent','nonfinite_enrichment','nonpositive_or_nonfinite_SE_IP')
    for reason in reason_order:
        remaining-=first_reason[reason]
        sequential.append({'reason':reason,'excluded_at_this_stage':first_reason[reason],'remaining':remaining})
    if remaining!=len(kept): raise IntegrityError('Sequential exclusion accounting mismatch')
    report={'input_reproduction_genes':len(rows),'original_membership_denominator':denominator,'eligible_genes':len(kept),'eligible_GeneIDs':[r['GeneID'] for r in kept],
            'P':len(rows),'S':denominator,'A':len(kept),'S_over_P':denominator/len(rows) if rows else None,'A_over_S':coverage,
            'reference_coverage':coverage,'reference_coverage_pass':denominator>0 and 5*len(kept)>=4*denominator,
            'excluded':excluded,'sequential_counts':sequential,'source_membership':source_report,
            'UTR_imputation_performed':False,'outcome_based_trimming':False,'replacement_transcripts_selected':False}
    return kept,report

def protein_digest(sequence): return hashlib.sha256(sequence.encode('ascii')).hexdigest()

def exact_distance(a,b):
    aligner=PairwiseAligner(mode='global',match_score=0,mismatch_score=-1,open_gap_score=-1,extend_gap_score=-1)
    score=-aligner.score(a,b)
    if score!=int(score): raise IntegrityError('Nonintegral C-backed edit score')
    return int(score)

def group_proteins(proteins_by_gene,*,limits=LIMITS,clock=time.monotonic,digest=protein_digest,progress=None):
    """Lossless5-mer filter and exact score. Any cap returns NO usable partition."""
    start=clock()
    counts={'genes':len(proteins_by_gene),'unique_proteins':0,'posting_pair_events_upper_bound':0,
            'posting_pair_events':0,'candidate_pairs_seen':0,'exact_pairs_completed':0,'alignment_cells_completed':0}
    def receipt(status,cap=None):
        return {'status':status,'cap':cap,'counts':dict(counts),'elapsed_seconds':clock()-start,
                'limits':asdict(limits),'components_complete':status=='complete',
                'elapsed_limit_semantics':'Acceptance cap checked throughout and after every native score, not hard preemption of a C call',
                'partial_components_must_not_be_used':status!='complete'}
    def stop(cap): raise ResourceLimit(cap,receipt('resource_inconclusive',cap))
    def check_time():
        if clock()-start>limits.elapsed_seconds: stop('elapsed_seconds')
    def checkpoint():
        if progress is not None: progress(receipt('in_progress'))
    genes=sorted(proteins_by_gene)
    if any(not isinstance(g,str) or not re.fullmatch(r'ENSG\d+',g) for g in genes): raise IntegrityError('Protein gene IDs must be exact stable ENSG IDs')
    groups={}; hashes={}
    for gene in genes:
        check_time(); sequence=proteins_by_gene[gene]
        if not isinstance(sequence,str) or not sequence or set(sequence)-set(AMINO): raise IntegrityError('Protein must be nonempty standard20-AA sequence')
        h=digest(sequence)
        if h in hashes and hashes[h]!=sequence: raise IntegrityError('Protein hash collision')
        hashes[h]=sequence; groups.setdefault(sequence,[]).append(gene)
    entries=sorted(((members[0],sequence,members) for sequence,members in groups.items()))
    counts['unique_proteins']=len(entries)
    checkpoint()
    parent={g:g for g in genes}
    def find(g):
        while parent[g]!=g: parent[g]=parent[parent[g]]; g=parent[g]
        return g
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b: parent[max(a,b)]=min(a,b)
    identity=[]
    for representative,sequence,members in entries:
        if len(members)>1: identity.append({'representative':representative,'members':members,'protein_sha256':digest(sequence)})
        for g in members[1:]: union(representative,g)
    counters=[]; index=defaultdict(list)
    for i,(_,sequence,_) in enumerate(entries):
        check_time(); counter=Counter(sequence[j:j+5] for j in range(max(0,len(sequence)-4))); counters.append(counter)
        for word in sorted(counter): index[word].append(i)
    upper=0
    for word in sorted(index):
        check_time(); n=len(index[word]); upper+=n*(n-1)//2
    counts['posting_pair_events_upper_bound']=upper
    checkpoint()
    if upper>limits.posting_pair_events: stop('posting_pair_events_preflight')
    edges=[]
    for i,(gene,a,_) in enumerate(entries):
        check_time(); common=Counter()
        for word in sorted(counters[i]):
            for j in index[word]:
                if j<=i: continue
                counts['posting_pair_events']+=1
                if counts['posting_pair_events']>limits.posting_pair_events: stop('posting_pair_events')
                check_time(); b=entries[j][1]; budget=max(len(a),len(b))//10
                if abs(len(a)-len(b))<=budget: common[j]+=min(counters[i][word],counters[j][word])
                if counts['posting_pair_events']%10000==0: checkpoint()
        for j in sorted(common):
            check_time(); other,b,_=entries[j]; budget=max(len(a),len(b))//10
            if budget==0 or max(max(0,len(a)-4),max(0,len(b)-4))-common[j]>5*budget: continue
            counts['candidate_pairs_seen']+=1
            if counts['candidate_pairs_seen']>limits.candidate_pairs: stop('candidate_pairs')
            cells=len(a)*len(b); proposed=counts['alignment_cells_completed']+cells
            if proposed>limits.alignment_cells:
                counts['rejected_next_alignment_cells']=cells; stop('alignment_cells')
            counts['pending_alignment_cells']=cells; checkpoint()
            distance=exact_distance(a,b)
            counts.pop('pending_alignment_cells')
            counts['exact_pairs_completed']+=1; counts['alignment_cells_completed']=proposed
            checkpoint()
            check_time()
            if 10*distance<=max(len(a),len(b)):
                union(gene,other); edges.append({'left_representative':gene,'right_representative':other,'distance':distance,'left_length':len(a),'right_length':len(b)})
    components=defaultdict(list)
    for gene in genes: check_time(); components[find(gene)].append(gene)
    result=receipt('complete')
    if result['elapsed_seconds']>limits.elapsed_seconds: stop('elapsed_seconds')
    result.update({'components':[{'group_ID':g,'GeneIDs':members} for g,members in sorted(components.items())],
                   'identity_groups':identity,'near_protein_edges':edges,'normalized_distance_denominator':'maximum full protein length',
                   'candidate_bound':'abs(lengthdiff)<=D and max(kmer_counts)-multiset_intersection<=5*D; D=floor(maxlength/10)',
                   'broader_family_independence_claimed':False})
    return result

def assign_partitions(components,*,digest=None):
    if digest is None: digest=lambda text:hashlib.sha256(text.encode('utf-8')).digest()
    records=[]; seen=set(); seen_digest={}
    for component in sorted(components,key=lambda c:c['group_ID']):
        group=component['group_ID']; members=component['GeneIDs']
        if not members or members!=sorted(set(members)) or group!=min(members): raise IntegrityError('Noncanonical component')
        if any(not re.fullmatch(r'ENSG\d+',g) or g in seen for g in members): raise IntegrityError('Component overlap/invalid GeneID')
        seen.update(members); values={}
        for kind in ('outer','inner'):
            if kind=='inner' and values['outer'][1]==0: continue
            text='cnot3-focus-v1|'+kind+'|'+group; raw=digest(text)
            if not isinstance(raw,bytes) or len(raw)!=32: raise IntegrityError('Expected complete32-byte SHA256 digest')
            if raw in seen_digest and seen_digest[raw]!=text: raise IntegrityError('Partition hash collision')
            seen_digest[raw]=text; values[kind]=(raw.hex(),int.from_bytes(raw,'big')%5)
        holdout=values['outer'][1]==0
        for gene in members:
            records.append({'GeneID':gene,'group_ID':group,'outer':'holdout' if holdout else 'training',
                            'inner_fold':None if holdout else values['inner'][1],
                            'outer_SHA256':values['outer'][0],
                            'inner_SHA256':None if holdout else values['inner'][0]})
    return sorted(records,key=lambda r:r['GeneID'])

def design_guards(eligibility_report,assignments):
    if len(assignments)!=eligibility_report['eligible_genes'] or len({r['GeneID'] for r in assignments})!=len(assignments): raise IntegrityError('Assignment/eligibility count mismatch')
    if {r['GeneID'] for r in assignments}!=set(eligibility_report['eligible_GeneIDs']): raise IntegrityError('Assignment/eligibility GeneID mismatch')
    train=[r for r in assignments if r['outer']=='training']; holdout=[r for r in assignments if r['outer']=='holdout']
    if any(r['outer'] not in ('training','holdout') or (r['outer']=='training' and r['inner_fold'] not in range(5)) or (r['outer']=='holdout' and r['inner_fold'] is not None) for r in assignments): raise IntegrityError('Invalid assignment state')
    groups=defaultdict(set)
    for r in assignments: groups[r['group_ID']].add((r['outer'],r['inner_fold']))
    if any(len(v)!=1 for v in groups.values()): raise IntegrityError('Component crosses a partition')
    checks={'reproduction_count_5481':eligibility_report['input_reproduction_genes']==5481,
            'reference_coverage_at_least80pct':eligibility_report['reference_coverage_pass'],
            'training_genes_at_least2000':len(train)>=2000,'holdout_genes_at_least300':len(holdout)>=300,
            'training_groups_at_least100':len({r['group_ID'] for r in train})>=100,
            'holdout_groups_at_least30':len({r['group_ID'] for r in holdout})>=30}
    folds=[]
    for fold in range(5):
        assessment=[r for r in train if r['inner_fold']==fold]; fitting=[r for r in train if r['inner_fold']!=fold]
        fitting_groups=len({r['group_ID'] for r in fitting})
        checks['inner'+str(fold)+'_assessment_genes_at_least100']=len(assessment)>=100
        checks['inner'+str(fold)+'_fitting_groups_at_least50']=fitting_groups>=50
        folds.append({'fold':fold,'assessment_genes':len(assessment),'assessment_groups':len({r['group_ID'] for r in assessment}),
                      'fitting_genes':len(fitting),'fitting_groups':fitting_groups})
    return {'status':'design_coverage_passed' if all(checks.values()) else 'insufficient_design_coverage',
            'checks':checks,'training_genes':len(train),'holdout_genes':len(holdout),
            'training_groups':len({r['group_ID'] for r in train}),'holdout_groups':len({r['group_ID'] for r in holdout}),
            'inner_folds':folds,'separate_minimum_inner_assessment_groups_chosen':False,'repartitioning_allowed':False}
