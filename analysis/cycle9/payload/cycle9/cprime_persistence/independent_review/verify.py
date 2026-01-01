"""Independent exact reconstruction from the already exposed B1 rows. Standard library only."""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python refused before input access')
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations,product
from collections import defaultdict,Counter
import json,csv,hashlib

ROOT=Path.cwd(); HERE=Path(__file__).resolve().parent
PILOT=HERE.parent
NAMES=('universal_three','universal_any','assignment_dependent','no_twofold','equality_excluded','continuous_twofold')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,msg):
    if not x:raise ValueError(msg)
def classifications(a,lo,hi):
    require(0<=lo<=hi<=3,'invalid assignment count')
    return {'universal_three':lo==3,'universal_any':lo>0,'assignment_dependent':lo!=hi,'no_twofold':hi==0,'equality_excluded':a>1,'continuous_twofold':a>=2}
def profile(vectors):
    require(len(vectors) in (2,3) and all(len(x)==3 and min(x)>0 for x in vectors),'invalid profile')
    lower=[min(x) for x in vectors];upper=[max(x) for x in vectors]
    a=max(lower)/min(upper)
    b=max(upper[i]/lower[j] for i in range(len(vectors)) for j in range(len(vectors)) if i!=j)
    results={}
    # Only common relabelling is quotiented out. A row's RNA and DNA stay paired.
    for tail in product(tuple(permutations(range(3))),repeat=len(vectors)-1):
        perm=((0,1,2),)+tail
        folds=tuple(max(vectors[i][perm[i][j]] for i in range(len(vectors)))/min(vectors[i][perm[i][j]] for i in range(len(vectors))) for j in range(3))
        results[perm]=(folds,sum(x>=2 for x in folds))
    counts=[x[1] for x in results.values()]
    return a,b,results,classifications(a,min(counts),max(counts))
def aggregate(rows):
    labels={r['intended'] for r in rows}
    counts={f:{'neither':0,'K562_only':0,'KBM7_only':0,'both':0} for f in NAMES}
    for label in labels:
        selected=[r for r in rows if r['intended']==label]
        require(len(selected)==2 and {r['cell'] for r in selected}=={'K562','KBM7'},'missing or duplicate cell')
        require(len({r['exon'] for r in selected})==1,'inconsistent exon')
        z={r['cell']:r for r in selected}
        for f in NAMES:
            key={(False,False):'neither',(True,False):'K562_only',(False,True):'KBM7_only',(True,True):'both'}[(z['K562'][f],z['KBM7'][f])]
            counts[f][key]+=1
    return counts
def synthetic():
    tests=[]
    def ok(v,msg):require(v,msg);tests.append(msg)
    for n in (2,3):
        vec=[tuple(F(2**(j+i)) for j in range(3)) for i in range(n)]
        a,b,f,flags=profile(vec)
        full=Counter()
        for perms in product(tuple(permutations(range(3))),repeat=n):
            k=sum(max(vec[i][perms[i][j]] for i in range(n))/min(vec[i][perms[i][j]] for i in range(n))>=2 for j in range(3))
            full[k]+=1
        reduced=Counter(x[1] for x in f.values())
        ok(full==Counter({k:6*v for k,v in reduced.items()}),'full permutation vs fixed first marker '+str(n))
        changed=[tuple(reversed(x)) for x in vec]
        ok(Counter(x[1] for x in profile(changed)[2].values())==reduced,'row-order invariance '+str(n))
    _,_,_,z=profile([(F(1),)*3,(F(2),)*3]);ok(z['universal_three'] and z['continuous_twofold'],'weak exact twofold tie accepted')
    _,_,_,z=profile([(F(1),)*3,(F(2)-F(1,10**20),)*3]);ok(not z['universal_any'] and not z['continuous_twofold'],'exact below two rejected')
    _,_,_,z=profile([(F(1),)*3,(F(2)+F(1,10**20),)*3]);ok(z['universal_three'],'exact above two accepted')
    a,_,_,z=profile([(F(1),F(4),F(16)),(F(2),F(8),F(32))]);ok(a<1 and z['universal_three'] and not z['continuous_twofold'],'discrete robustness does not imply continuous exclusion')
    ok(not classifications(F(1),0,0)['equality_excluded'],'equality tie retained')
    toy=[]
    for i,pair in enumerate(((False,False),(True,False),(False,True),(True,True))):
        for cell,v in zip(('K562','KBM7'),pair):toy.append({'intended':str(i),'cell':cell,'exon':17,**classifications(F(2 if v else 1),3 if v else 0,3 if v else 0)})
    ok(aggregate(toy)['universal_three']==dict(neither=1,K562_only=1,KBM7_only=1,both=1),'four cross-cell quadrants')
    bad=[dict(r) for r in toy];bad[0]['exon']=27
    for rows,name in ((toy[1:],'missing cell'),(toy+toy[:1],'duplicate cell'),(bad,'inconsistent exon')):
        try:aggregate(rows)
        except ValueError:tests.append(name+' rejected')
        else:raise ValueError(name+' incorrectly accepted')
    for v in ([],[(F(1),F(2),F(3))],[(F(1),F(2)),(F(1),F(2))],[(F(0),F(1),F(2)),(F(1),)*3]):
        try:profile(v)
        except ValueError:tests.append('invalid profile rejected')
        else:raise ValueError('invalid profile incorrectly accepted')
    paired=[F(x,y) for x,y in ((1,1),(4,2),(9,3))]
    unpaired=[F(x,y) for x,y in zip((1,4,9),(3,2,1))]
    ok(paired!=unpaired,'independent numerator denominator permutation changes estimand')
    return tests
def main():
    tests=synthetic();out=HERE/'run_v1';require(not out.exists(),'output exists')
    pins={'work/research_cycle8/cprime_results/primary_v1/result.json':'eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf','work/research_cycle8/cprime_followup/results_v1/result.json':'851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e','work/research_cycle9/cprime_persistence/protocol.md':'13274dd480f4c36b278370a1dd3907e86473beae0a59e5d0def1a01136fba613','work/research_cycle9/cprime_persistence/analyze.py':'f99e6d1cb7397c931d7b0b5ac44bfb6d74da824f971de6c1f5a9b4709f1d3aff','work/research_cycle9/cprime_persistence/results_v1/result.json':'46fdaffece0fdc6d6803d6e492ccec5a301954978fde39b96f04e689d76de6a8'}
    for p,h in pins.items():require(sha(ROOT/p)==h,'input hash changed: '+p)
    primary=json.loads((ROOT/'work/research_cycle8/cprime_results/primary_v1/result.json').read_text())
    buckets=defaultdict(lambda:defaultdict(list));ids=set();source_row_count=0
    for row in primary['source_rows']:
        require(row['source_row'] not in ids,'duplicate source row');ids.add(row['source_row'])
        require(row['usable'] and row['missing_reason'] is None,'unusable source row')
        ratio=F(row['RNA'])/F(row['DNA']);require(ratio==F(row['ratio']),'paired source ratio mismatch')
        buckets[(row['intended'],row['cell'])][row['state']].append((row['source_row'],ratio,int(row['exon'])))
        source_row_count+=1
    require(source_row_count==396 and len(buckets)==60,'incomplete B1 rows')
    own=[];assignment_records={};all_exons={}
    for (target,cell),states in sorted(buckets.items()):
        names=sorted(states);values=[sorted(states[n]) for n in names]
        require(all(len(x)==3 for x in values),'missing preparation row')
        exons={v[2] for x in values for v in x};require(len(exons)==1,'inconsistent exon map')
        a,b,assignments,flags=profile([tuple(v[1] for v in x) for x in values])
        lo=min(v[1] for v in assignments.values());hi=max(v[1] for v in assignments.values())
        own.append({'intended':target,'exon':next(iter(exons)),'cell':cell,'A_exact':str(a),'B_exact':str(b),'k_min':lo,'k_max':hi,'assignments':len(assignments),**flags})
        assignment_records[(target,cell)]=(names,[[x[0] for x in v] for v in values],assignments)
    # Independent reconstruction is complete before decoding source summaries or author answers.
    cross=aggregate(own);strata=[];omissions=[]
    for exon in sorted({r['exon'] for r in own}):
        local=[r for r in own if r['exon']==exon];other=[r for r in own if r['exon']!=exon]
        strata.append({'exon':exon,'profiles':len(local),'targets':len(local)//2,'profile_totals':{f:sum(r[f] for r in local) for f in NAMES},'paired_cell':aggregate(local)})
        omissions.append({'omitted_exon':exon,'remaining_targets':len(other)//2,'profile_totals':{f:sum(r[f] for r in other) for f in NAMES},'paired_cell':aggregate(other)})
    prior=json.loads((ROOT/'work/research_cycle8/cprime_followup/results_v1/result.json').read_text());seen=set();assignment_checks=0;fold_checks=0
    for group in prior['groups']:
        key=(group['intended'],group['cell']);require(key not in seen,'duplicate B2 group');seen.add(key)
        names,rowids,expected=assignment_records[key];require(names==group['markers'] and rowids==group['source_rows'],'B2 row identity mismatch')
        row=next(r for r in own if (r['intended'],r['cell'])==key)
        require(row['A_exact']==group['A_exact'] and row['B_exact']==group['B_exact'],'B2 interval ratio mismatch')
        received={}
        for item in group['assignments']:
            perm=tuple(tuple(x) for x in item['permutations']);require(perm not in received,'duplicate B2 assignment')
            received[perm]=(tuple(F(x) for x in item['folds']),item['preparations_at_least_twofold']);fold_checks+=len(item['folds'])
        require(received==expected,'B2 assignment mismatch');assignment_checks+=len(expected)
        require((row['k_min'],row['k_max'])==(group['minimum_preparations_at_least_twofold'],group['maximum_preparations_at_least_twofold']),'B2 count summary mismatch')
    require(seen==set(buckets) and assignment_checks==720 and fold_checks==2160,'incomplete assignment coverage')
    stored=json.loads((PILOT/'results_v1/result.json').read_text())
    require(stored['profiles']==60 and stored['targets']==30 and stored['assignments']==720 and stored['exons']==[17,25,27,39,42],'denominator mismatch')
    require(stored['paired_cell']==cross and stored['by_exon']==strata and stored['all_leave_one_exon_out']==omissions,'aggregate mismatch')
    with (PILOT/'results_v1/all_60_profiles.csv').open(newline='') as f:exported=list(csv.DictReader(f))
    require(len(exported)==60,'profile export incomplete')
    for row in own:
        found=[x for x in exported if x['intended']==row['intended'] and x['cell']==row['cell']];require(len(found)==1,'profile export duplicate/missing')
        require({k:str(v) for k,v in row.items()}==found[0],'profile export mismatch')
    for p,h in pins.items():require(sha(ROOT/p)==h,'input mutated during audit')
    both=[r['intended'] for r in own if r['cell']=='K562' and r['universal_three'] and next(x for x in own if x['intended']==r['intended'] and x['cell']=='KBM7')['universal_three']]
    out.mkdir();result={'status':'PASS_EXACT_INDEPENDENT_RECONSTRUCTION','profiles':60,'targets':30,'source_rows_checked':396,'assignments_checked':720,'folds_checked':2160,'paired_cell':cross,'by_exon':strata,'all_leave_one_exon_out':omissions,'both_cell_universal_three_targets':both,'synthetic_checks':tests,'mismatches':[],'input_hashes_before_after':pins,'reviewer_code_sha256':sha(Path(__file__)),'no_production_import':True,'scope':'Exposed cohort descriptive cross-tabulation; no new source outcomes, no independent biological replication or calibrated confidence.'}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'independent_profiles.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(own[0]));w.writeheader();w.writerows(own)
    print(json.dumps({'status':result['status'],'paired_cell':cross,'both_cell_targets':both,'synthetic_checks':len(tests),'verification_sha256':sha(out/'verification.json')},indent=2))
if __name__=='__main__':main()
