"""Independent exact certificate and all-four-state local-window checks.

Does not import audit_design.py and does not read any phenotype data.
"""
import argparse, collections, hashlib, itertools, json
from fractions import Fraction
from pathlib import Path

def hashfile(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def modular_rank(rows,p=1000003):
    a=[list(r) for r in rows];n=0
    for j in range(len(a[0])):
        i=next((i for i in range(n,len(a)) if a[i][j]%p),None)
        if i is None:continue
        a[i],a[n]=a[n],a[i];inv=pow(a[n][j],-1,p)
        a[n]=[(v*inv)%p for v in a[n]]
        for k in range(n+1,len(a)):
            v=a[k][j];a[k]=[(x-v*y)%p for x,y in zip(a[k],a[n])]
        n+=1
        if n==len(a):break
    return n

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if not __debug__:raise RuntimeError('Assertions required')
    if a.out.exists():raise FileExistsError(a.out)
    paths=[a.inputs/n for n in ['mapping.json','certificates.json','result.json']]
    hashes={p.name:hashfile(p) for p in paths}
    mapped,cert,result=[json.loads(p.read_text()) for p in paths]
    X=cert['matrix'];left=[[Fraction(z) for z in v] for v in cert['left_nullspace']]
    right=[[Fraction(z) for z in v] for v in cert['right_nullspace']]
    checks=0
    for v in left:
        for j in range(len(X[0])):
            if sum(v[i]*X[i][j] for i in range(len(X))):raise AssertionError('Left certificate')
            checks+=1
    for v in right:
        for row in X:
            if sum(x*y for x,y in zip(row,v)):raise AssertionError('Right certificate')
            checks+=1
    # The supplied null vectors have unique free-variable unit coordinates,
    # proving independence over Q without a second rational elimination.
    def independent_pivots(vs):
        chosen=[]
        for i,v in enumerate(vs):
            j=next((j for j,z in enumerate(v) if z==1 and all(w[j]==0 for k,w in enumerate(vs) if k!=i)),None)
            if j is None:raise AssertionError('Missing null-basis unit coordinate')
            chosen.append(j)
        return chosen
    lp=independent_pivots(left);rp=independent_pivots(right)
    mr=modular_rank(X)
    # Rank mod p <= rank over Q <= nrows - independent left-null dimension.
    if mr!=len(X)-len(left) or mr!=len(X[0])-len(right):raise AssertionError('Exact rank sandwich')
    quartets=0;common_covers=0;two_site=[]
    for e in sorted({m['exon'] for m in mapped}):
        rows=[m for m in mapped if m['exon']==e]
        for four in itertools.combinations(rows,4):
            quartets+=1
            lo=max(m['observed_start'] for m in four);hi=min(m['observed_end'] for m in four)
            changes={v['pos'] for m in four for v in m['changes']}
            if lo>hi or not all(lo<=p<=hi for p in changes):continue
            common_covers+=1
            columns={p:{m['observed_bases'][str(p)] for m in four} for p in range(lo,hi+1)}
            var=[p for p,v in columns.items() if len(v)>1]
            if len(var)!=2 or any(len(columns[p])!=2 for p in var):continue
            comb={tuple(m['observed_bases'][str(p)] for p in var) for m in four}
            if len(comb)==4:two_site.append({'ids':[m['id'] for m in four],'positions':var})
    if len(two_site)!=result['verified_two_site_common_window_factorials']:raise AssertionError('Factorial enumeration disagrees')
    rect=result['rectangles'][0];rows=[mapped[i] for i in rect['row_indices']]
    lo=rect['common_observed_start'];hi=rect['common_observed_end']
    words=[''.join(m['observed_bases'][str(p)] for p in range(lo,hi+1)) for m in rows]
    hd=[[sum(x!=y for x,y in zip(u,v)) for v in words] for u in words]
    if hd!=[[0,2,2,4],[2,0,4,2],[2,4,0,2],[4,2,2,0]]:raise AssertionError('Composite contrast distance')
    # A specified intended effect differs along a verified null direction.
    gauge=[Fraction(z) for z in cert['gauge_witness']]
    nonident=[cert['column_labels'][i] for i,z in enumerate(gauge) if z]
    # Arbitrary strictly increasing output calibration can manufacture the sole
    # observed-scale rectangle contrast from an additive latent field.
    latent=[0,1,2,3];square=[0,1,4,9]
    signs=rect['contrast_signs']
    if sum(s*x for s,x in zip(signs,latent))!=0 or sum(s*x for s,x in zip(signs,square))!=4:raise AssertionError('Monotone-scale counterexample')
    report={'status':'PASS','input_sha256':hashes,'source_sha256':hashfile(Path(__file__)),
      'rank_mod_prime':mr,'prime':1000003,'exact_rank_lower_bound':mr,
      'exact_rank_upper_bound_left':len(X)-len(left),'exact_rank_upper_bound_right':len(X[0])-len(right),
      'left_null_independence_coordinates':lp,'right_null_independence_coordinates':rp,
      'exact_null_dot_products':checks,'all_within_exon_quartets':quartets,
      'quartets_with_common_word_covering_every_specified_change':common_covers,
      'two_site_nucleotide_factorials':two_site,'rectangle_common_words':words,'rectangle_hamming_distances':hd,
      'gauge_nonidentifiable_alleles':nonident,
      'scale_counterexample':{'latent_additive':latent,'strictly_order_preserving_square_output':square,'output_rectangle_contrast':4},
      'outcome_values_read':False,'limit':'Independent algorithm and representation, same author and inherited project context. No blinded or separately authored biological review claimed.'}
    if hashes!={p.name:hashfile(p) for p in paths}:raise AssertionError('Input changed')
    a.out.mkdir(parents=True);(a.out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','exact_rank_lower_bound','exact_rank_upper_bound_left','exact_null_dot_products','all_within_exon_quartets','quartets_with_common_word_covering_every_specified_change','rectangle_hamming_distances']}))

if __name__=='__main__':main()
