from pathlib import Path
from fractions import Fraction
from collections import Counter, defaultdict
import csv, hashlib, json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'work/research_cycle8/cprime_followup/results_v1/result.json'
PRIMARY=ROOT/'work/research_cycle8/cprime_results/primary_v1/result.json'
FLAGS=['universal_three','universal_any','assignment_dependent','no_twofold','equality_excluded','continuous_twofold']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def classify(A,kmin,kmax):
    assert 0<=kmin<=kmax<=3
    return dict(zip(FLAGS,[kmin==3,kmin>=1,kmin<kmax,kmax==0,A>1,A>=2]))
def tables(rows):
    by=defaultdict(dict)
    for r in rows:
        assert r['cell'] not in by[r['intended']]
        by[r['intended']][r['cell']]=r
    result={f:Counter({'neither':0,'K562_only':0,'KBM7_only':0,'both':0}) for f in FLAGS}
    for target,cells in by.items():
        assert set(cells)=={'K562','KBM7'},target
        assert cells['K562']['exon']==cells['KBM7']['exon']
        for f in FLAGS:
            a,b=cells['K562'][f],cells['KBM7'][f]
            result[f]['both' if a and b else 'K562_only' if a else 'KBM7_only' if b else 'neither']+=1
    return {f:dict(x) for f,x in result.items()}

def synthetic():
    rows=[]
    for i,(a,b) in enumerate([(False,False),(False,True),(True,False),(True,True)]):
        for c,z in [('K562',a),('KBM7',b)]:
            rows.append({'intended':str(i),'cell':c,'exon':17,**classify(Fraction(2 if z else 1),3 if z else 0,3 if z else 0)})
    t=tables(rows)
    for f in ['universal_three','universal_any','equality_excluded','continuous_twofold']:
        assert t[f]==dict(neither=1,K562_only=1,KBM7_only=1,both=1)
    assert classify(Fraction(2),0,0)['continuous_twofold']
    assert not classify(Fraction(1),0,0)['equality_excluded']
    for bad in [rows+rows[:1],rows[1:]]:
        try:tables(bad)
        except AssertionError:pass
        else:raise AssertionError('duplicate/missing cell accepted')
    bad=[dict(r) for r in rows];bad[0]['exon']=27
    try:tables(bad)
    except AssertionError:pass
    else:raise AssertionError('inconsistent exon accepted')
    return {'status':'PASS','paired_cell_quadrants':4,'invalid_join_guards':3,'exact_threshold_checks':2}

def main():
    test=synthetic()
    (HERE/'synthetic_receipt.json').write_text(json.dumps(test,indent=2)+'\n')
    assert sha(SOURCE)=='851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e'
    assert sha(PRIMARY)=='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'
    source=json.loads(SOURCE.read_text()); primary=json.loads(PRIMARY.read_text())
    mapping=defaultdict(set)
    for r in primary['source_rows']:mapping[r['intended']].add(int(r['exon']))
    rows=[];assignment_total=0
    for g in source['groups']:
        assert len(mapping[g['intended']])==1
        counts=[sum(Fraction(v)>=2 for v in a['folds']) for a in g['assignments']]
        assert all(k==a['preparations_at_least_twofold'] for k,a in zip(counts,g['assignments']))
        mn,mx=min(counts),max(counts)
        assert mn==g['minimum_preparations_at_least_twofold'] and mx==g['maximum_preparations_at_least_twofold']
        assignment_total+=len(counts)
        rows.append({'intended':g['intended'],'exon':next(iter(mapping[g['intended']])),'cell':g['cell'],
                     'A_exact':g['A_exact'],'B_exact':g['B_exact'],'k_min':mn,'k_max':mx,
                     'assignments':len(counts),**classify(Fraction(g['A_exact']),mn,mx)})
    assert len(rows)==60 and assignment_total==720
    cross=tables(rows); exons=sorted({r['exon'] for r in rows});assert exons==[17,25,27,39,42]
    strat=[];omissions=[]
    for e in exons:
        er=[r for r in rows if r['exon']==e]
        strat.append({'exon':e,'profiles':len(er),'targets':len(er)//2,'profile_totals':{f:sum(r[f] for r in er) for f in FLAGS},'paired_cell':tables(er)})
        other=[r for r in rows if r['exon']!=e]
        omissions.append({'omitted_exon':e,'remaining_targets':len(other)//2,'profile_totals':{f:sum(r[f] for r in other) for f in FLAGS},'paired_cell':tables(other)})
    out=HERE/'results_v1';out.mkdir(exist_ok=False)
    with (out/'all_60_profiles.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'status':'COMPLETE_EXPLORATORY_REQUIRES_INDEPENDENT_REVIEW','profiles':60,'targets':30,'cells':['K562','KBM7'],'exons':exons,'assignments':assignment_total,'paired_cell':cross,'by_exon':strat,'all_leave_one_exon_out':omissions,'limits':'Finite exposed cohort; no calibrated uncertainty, prevalence, independent validation, shared ordering or mechanism inference.'}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    manifest={'protocol_sha256':sha(HERE/'protocol.md'),'code_sha256':sha(Path(__file__)),'input_sha256':sha(SOURCE),'primary_metadata_sha256':sha(PRIMARY),'result_sha256':sha(out/'result.json'),'synthetic':test}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'paired_cell':cross,'by_exon':[{k:e[k] for k in ['exon','profiles','targets','profile_totals']} for e in strat]}))

if __name__=='__main__':main()
