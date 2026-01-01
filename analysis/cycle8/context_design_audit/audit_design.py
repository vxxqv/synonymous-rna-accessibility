"""Exact metadata-only cPRIME descriptor and common-window design audit."""
from __future__ import annotations
import argparse, collections, hashlib, itertools, json, platform, re, time
from fractions import Fraction as Q
from pathlib import Path

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rref(rows):
    a = [[Q(x) for x in r] for r in rows]
    if not a: return [], []
    piv=[]; k=0
    for j in range(len(a[0])):
        z=next((i for i in range(k,len(a)) if a[i][j]),None)
        if z is None: continue
        a[k],a[z]=a[z],a[k]
        div=a[k][j]; a[k]=[x/div for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]:
                f=a[i][j]; a[i]=[u-f*v for u,v in zip(a[i],a[k])]
        piv.append(j); k+=1
        if k==len(a): break
    return a,piv

def rank(rows): return len(rref(rows)[1])
def transpose(rows): return list(map(list,zip(*rows)))
def nullspace(rows):
    a,piv=rref(rows); n=len(rows[0]); out=[]
    for j in range(n):
        if j in piv: continue
        v=[Q(0)]*n;v[j]=Q(1)
        for i,p in enumerate(piv):v[p]=-a[i][j]
        if any(sum(x*y for x,y in zip(r,v)) for r in rows): raise ValueError('Invalid exact null vector')
        out.append(v)
    return out

def frac(v): return [str(x) for x in v]
def parse_id(s):
    m=re.fullmatch(r'(chr\w+)-(\d+)-([ACGT])-([ACGT])',s)
    if not m or m[3]==m[4]:raise ValueError('Not a documented SNP: '+s)
    return {'id':s,'chrom':m[1],'pos':int(m[2]),'ref':m[3],'alt':m[4]}

def align(state,design):
    changes=sorted([parse_id(state['var_id']),parse_id(state['syn_var'])],key=lambda x:x['pos'])
    if len({x['pos'] for x in changes})!=2:raise ValueError('Coincident descriptors')
    raw=design['Target_sequence']; refs=[];alts=[]; brackets=[];last=0
    for m in re.finditer(r'\(([ACGT])/([ACGT])\)',raw):
        fixed=raw[last:m.start()]
        if re.search('[^ACGT]',fixed):raise ValueError('Unsupported context syntax')
        refs.extend(fixed);alts.extend(fixed)
        brackets.append((len(refs),m[1],m[2]));refs.append(m[1]);alts.append(m[2]);last=m.end()
    tail=raw[last:]
    if re.search('[^ACGT]',tail):raise ValueError('Unsupported tail syntax')
    refs.extend(tail);alts.extend(tail)
    if len(brackets)!=2:raise ValueError('Expected exactly two brackets')
    origins={v['pos']-b[0] for v,b in zip(changes,brackets)}
    if len(origins)!=1 or any((v['ref'],v['alt'])!=(b[1],b[2]) for v,b in zip(changes,brackets)):
        raise ValueError('Reference alignment or orientation not authenticated')
    origin=origins.pop(); ref=''.join(refs);alt=''.join(alts)
    word=state['motifs']['DNA'][0];hits=[m.start() for m in re.finditer('(?='+word+')',alt)]
    if len(hits)!=1:raise ValueError('Counted DNA word is not uniquely supported by source design')
    start=origin+hits[0];end=start+len(word)-1
    if not all(start<=v['pos']<=end for v in changes):raise ValueError('Word omits a specified change')
    if any(word[v['pos']-start]!=v['alt'] for v in changes):raise ValueError('Alternate mismatch')
    observed={str(start+i):b for i,b in enumerate(word)}
    return {'id':state['combined_var_id'],'exon':int(state['exon_number']),
      'intended_allele':state['var_id'],'marker_allele':state['syn_var'],
      'alleles':[v['id'] for v in changes],'changes':changes,'design_origin':origin,
      'design_end':origin+len(ref)-1,'design_reference_sha256':hashlib.sha256(ref.encode()).hexdigest(),
      'observed_DNA_word':word,'observed_start':start,'observed_end':end,'observed_bases':observed,
      'unknown_context':'Every coordinate outside the counted DNA word is unobserved for this local state. Source design is not an observed complete haplotype.',
      'reference_word':ref[hits[0]:hits[0]+len(word)],'design_reference':ref}

def graph_rank(alleles,pairs):
    adj={a:set() for a in alleles}
    for a,b in pairs:adj[a].add(b);adj[b].add(a)
    unseen=set(alleles); comps=[]
    while unseen:
        start=min(unseen);color={start:0};q=[start];bip=True
        for u in q:
            for v in sorted(adj[u]):
                if v not in color:color[v]=1-color[u];q.append(v)
                elif color[v]==color[u]:bip=False
        unseen-=set(q);e=sum(len(adj[v]) for v in q)//2
        comps.append({'alleles':sorted(q),'vertices':len(q),'edges':e,'bipartite':bip,
                      'rank':len(q)-int(bip),'color_classes':[[v for v in sorted(q) if color[v]==c] for c in [0,1]]})
    return sum(c['rank'] for c in comps),comps

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if not __debug__:raise ValueError('Optimization disabled assertions')
    if args.out.exists():raise FileExistsError(args.out)
    starttime=time.perf_counter()
    base=args.root/'work/research_cycle8'; gate=base/'cprime_gate'
    paths={'eligibility':base/'cprime_pilot/eligibility.json','design':gate/'multi_exon_selected_metadata.json',
           'assay_metadata':gate/'selected_assay_metadata.json','search_words':gate/'search_sequence_metadata.json',
           'contract':Path(__file__).with_name('design_contract.json')}
    hashes={k:sha(p) for k,p in paths.items()}
    elig=json.loads(paths['eligibility'].read_text()); designs=json.loads(paths['design'].read_text())['rows']
    for role,name in [('assay_metadata','selected_assay_metadata.json'),('search_words','search_sequence_metadata.json')]:
        if hashes[role]!=elig['input_sha256'][name]:raise ValueError('Eligibility source drift')
    byid=collections.defaultdict(list)
    for d in designs:byid[d['combined_var_id']].append(d)
    mapped=[]
    for key,s in sorted(elig['states'].items()):
        ds=byid[key]
        if not ds or len({d['Target_sequence'] for d in ds})!=1:raise ValueError('Missing or ambiguous source design')
        mapped.append(align(s,ds[0]))
    if len(mapped)!=66:raise ValueError('Fixed inventory changed')
    # Source reference consistency is tested only where designed windows overlap.
    refmap={}
    for m in mapped:
        for i,b in enumerate(m['design_reference']):
            key=(m['changes'][0]['chrom'],m['design_origin']+i)
            if key in refmap and refmap[key]!=b:raise ValueError('Conflicting source reference')
            refmap[key]=b
    exons=sorted({m['exon'] for m in mapped}); alleles=sorted({a for m in mapped for a in m['alleles']})
    labels=[f'exon:{e}' for e in exons]+alleles
    X=[[int(m['exon']==e) for e in exons]+[int(a in m['alleles']) for a in alleles] for m in mapped]
    rr,piv=rref(X); r=len(piv); left=nullspace(transpose(X));right=nullspace(X)
    gr,comps=graph_rank(alleles,[m['alleles'] for m in mapped])
    if r!=gr:raise ValueError('Independent incidence rank mismatch')
    # Every row has two within-exon descriptors; the intercept columns are redundant.
    for e in exons:
        ae={a for m in mapped if m['exon']==e for a in m['alleles']}
        if any(Q(sum(row[len(exons)+alleles.index(a)] for a in ae),2)!=row[exons.index(e)] for row in X):raise ValueError('Intercept identity failed')
    # Full recorded-pair terms, not post-outcome selected interactions.
    pairkeys=sorted({tuple(m['alleles']) for m in mapped})
    Z=[row+[int(set(p).issubset(m['alleles'])) for p in pairkeys] for row,m in zip(X,mapped)]
    zr=rank(Z)
    # Complement each descriptor within its exon; this is an invertible basis operation.
    allele_exon={a:next(m['exon'] for m in mapped if a in m['alleles']) for a in alleles}
    Xrec=[]
    for row in X:
        Xrec.append(row[:len(exons)]+[row[exons.index(allele_exon[a])]-row[len(exons)+j] for j,a in enumerate(alleles)])
    if rank(Xrec)!=r or rank([a+b for a,b in zip(X,Xrec)])!=r:raise ValueError('Reference recoding changed span')
    # Cell models, constructed directly rather than rank inferred only by a formula.
    zeros=[0]*len(labels)
    shared=X+X
    offsets=[row+[0]*len(exons) for row in X]+[row+row[:len(exons)] for row in X]
    context=[row+zeros for row in X]+[zeros+row for row in X]
    context_pair=[row+[0]*len(Z[0]) for row in Z]+[[0]*len(Z[0])+row for row in Z]
    ranks={'strict_shared':rank(shared),'shared_slopes_context_exon_offsets':rank(offsets),
           'context_specific_additive':rank(context),'context_specific_all_observed_pairs':rank(context_pair)}
    if not(ranks['strict_shared']==r and ranks['context_specific_additive']==2*r and ranks['context_specific_all_observed_pairs']==2*zr):raise ValueError('Block model identity')
    # Descriptor rectangles with two two-level choices, followed by exact common-window challenge.
    edge_to_row={frozenset(m['alleles']):i for i,m in enumerate(mapped)}
    if len(edge_to_row)!=len(mapped):raise ValueError('Duplicate physical descriptor pair')
    adjacency={a:set() for a in alleles}
    for a,b in pairkeys:adjacency[a].add(b);adjacency[b].add(a)
    rectangles={}
    for a,b in itertools.combinations(alleles,2):
        for c,d in itertools.combinations(sorted(adjacency[a]&adjacency[b]),2):
            inds=[edge_to_row[frozenset(p)] for p in [(a,c),(a,d),(b,c),(b,d)]]
            if len(set(inds))!=4:continue
            key=tuple(sorted(inds))
            if key in rectangles:continue
            rows=[mapped[i] for i in inds];lo=max(m['observed_start'] for m in rows);hi=min(m['observed_end'] for m in rows)
            variable=sorted({int(pos) for pos in range(lo,hi+1) if len({m['observed_bases'][str(pos)] for m in rows})>1}) if lo<=hi else []
            all_specified={v['pos'] for m in rows for v in m['changes']}
            complete=lo<=hi and all(lo<=p<=hi for p in all_specified)
            signs=[1,-1,-1,1]
            v=[0]*len(mapped)
            for i,s in zip(inds,signs):v[i]=s
            if any(sum(v[i]*X[i][j] for i in range(len(X))) for j in range(len(labels))):raise ValueError('Rectangle is not an additive contrast')
            # A proper two nucleotide-site factorial requires two variable coordinates,
            # exactly two levels each, and all four combinations within the common word.
            combos={tuple(m['observed_bases'][str(p)] for p in variable) for m in rows}
            levels=[sorted({m['observed_bases'][str(p)] for m in rows}) for p in variable]
            nucleotide_factorial=complete and len(variable)==2 and all(len(z)==2 for z in levels) and len(combos)==4
            rectangles[key]={'row_indices':inds,'ids':[m['id'] for m in rows],'descriptor_choices':[[a,b],[c,d]],
              'contrast_signs':signs,'common_observed_start':lo,'common_observed_end':hi,
              'common_observed_covers_all_specified_changes':complete,'variable_coordinates':variable,
              'variable_levels':levels,'two_site_nucleotide_factorial':nucleotide_factorial,
              'all_unshown_context_unknown':True}
    # Parameter gauge witness on the first bipartite incidence component.
    comp=next(c for c in comps if c['bipartite']);gauge=[Q(0)]*len(labels)
    for side,sign in zip(comp['color_classes'],[1,-1]):
        for a in side:gauge[len(exons)+alleles.index(a)]=Q(sign)
    if any(sum(a*b for a,b in zip(row,gauge)) for row in X):raise ValueError('Gauge witness failed')
    exsummary=[]
    for e in exons:
        inds=[i for i,m in enumerate(mapped) if m['exon']==e]
        local=[X[i] for i in inds]; localrank=rank(local)
        exsummary.append({'exon':e,'states':len(inds),'groups':len({mapped[i]['intended_allele'] for i in inds}),'rank':localrank,'residual_df':len(inds)-localrank})
    # Known-answer design checks, no observed phenotypes.
    rectangle=[[1,1,0,1,0],[1,1,0,0,1],[1,0,1,1,0],[1,0,1,0,1]]
    tree=rectangle[:3]; interaction=[0,0,0,1]
    if rank(tree)!=3 or rank(rectangle)!=3 or rank([x+[y] for x,y in zip(rectangle,interaction)])!=4:raise ValueError('Synthetic rank checks failed')
    missing={'00':0,'10':1,'01':1} # Both additive 11=2 and interacting 11=7 fit these rows.
    synthetics={'tree_rows':3,'tree_rank':3,'rectangle_rows':4,'rectangle_additive_rank':3,
       'rectangle_interaction_rank':4,'missing_corner_observations':missing,'compatible_missing_values':[2,7],
       'paired_context_rank_additive':rank([x+[0]*5 for x in rectangle]+[[0]*5+x for x in rectangle]),
       'all_exact_checks_passed':True}
    # Coefficient-vs-prediction estimability is explicit and not inferred from a fit.
    certificates={'row_ids':[m['id'] for m in mapped],'column_labels':labels,'matrix':X,'pivot_columns':piv,
       'rref':[[str(v) for v in row] for row in rr],'left_nullspace':[frac(v) for v in left],
       'right_nullspace':[frac(v) for v in right],'gauge_witness':frac(gauge),
       'interaction_pair_labels':[list(p) for p in pairkeys]}
    result={'status':'metadata_design_audit_complete','outcome_values_read':False,'states':len(mapped),
       'groups':len(elig['groups']),'exons':exons,'descriptor_alleles':len(alleles),'additive_columns':len(labels),
       'additive_rank':r,'additive_residual_df':len(mapped)-r,'all_pair_rank':zr,'all_pair_columns':len(Z[0]),
       'interaction_span_gain':zr-r,'context_ranks':ranks,'context_rows':2*len(mapped),'by_exon':exsummary,
       'incidence_components':comps,'independent_incidence_rank':gr,'rectangles':list(rectangles.values()),
       'descriptor_rectangles':len(rectangles),'verified_two_site_common_window_factorials':sum(v['two_site_nucleotide_factorial'] for v in rectangles.values()),
       'reference_recoding_span_equal':True,'source_inputs_sha256':hashes,'synthetic_checks':synthetics,
       'limits':['Recorded pair descriptors are not complete genotypes.','All context outside each counted DNA word is unknown.',
                 'Cell contexts and source rows do not imply independent preparations.','No source expression values or LP fits were read or run.']}
    if {k:sha(p) for k,p in paths.items()}!=hashes:raise ValueError('Source changed during audit')
    args.out.mkdir(parents=True)
    for m in mapped:m.pop('design_reference')
    for name,obj in [('mapping.json',mapped),('certificates.json',certificates),('result.json',result)]:
        (args.out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    receipt={'source_sha256':sha(Path(__file__)),'input_sha256':hashes,'outputs_sha256':{n:sha(args.out/n) for n in ['mapping.json','certificates.json','result.json']},
       'python':platform.python_version(),'runtime_seconds':time.perf_counter()-starttime,'all_exact_checks_passed':True,'outcomes_read':False}
    (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['states','groups','descriptor_alleles','additive_rank','additive_residual_df','interaction_span_gain','descriptor_rectangles','verified_two_site_common_window_factorials','context_ranks','by_exon']}))

if __name__=='__main__':main()
