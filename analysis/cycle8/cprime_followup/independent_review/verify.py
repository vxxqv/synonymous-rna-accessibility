"""Separately authored exact cPRIME follow-up audit, no production imports."""
from __future__ import annotations
import argparse, collections, hashlib, itertools, json, platform, random, time
from fractions import Fraction as F
from pathlib import Path

B1_SHA='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'
RESULT_SHA='851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e'
N=0
def require(v,why):
    global N
    if not v:raise AssertionError(why)
    N+=1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,obj):p.write_text(json.dumps(obj,indent=2,default=str)+'\n',encoding='utf-8')

def maximum_matching(adj,nright=3):
    """Augmenting paths, independent of full-bijection enumeration."""
    owner=[None]*nright
    def augment(u,seen):
        for v in sorted(adj[u]):
            if v in seen:continue
            seen.add(v)
            if owner[v] is None or augment(owner[v],seen):owner[v]=u;return True
        return False
    for u in range(len(adj)):augment(u,set())
    return sum(v is not None for v in owner),[(u,v) for v,u in enumerate(owner) if u is not None]

def assignment_oracle(values):
    m=len(values);ids=range(3);records={}
    # Build complete row matchings directly from selected tuples, not an author helper.
    for choice in itertools.product(tuple(itertools.permutations(ids)),repeat=m-1):
        perms=((0,1,2),)+choice
        rows=list(zip(*(tuple(vals[i] for i in p) for vals,p in zip(values,perms))))
        folds=tuple(max(row)/min(row) for row in rows);count=sum(f>=2 for f in folds)
        records[perms]=(folds,count)
    if m==2:
        good=[{j for j in ids if max(values[0][i],values[1][j])<2*min(values[0][i],values[1][j])} for i in ids]
        bad=[set(ids)-s for s in good]
        ng,gm=maximum_matching(good);nb,bm=maximum_matching(bad)
        bounds=(3-ng,nb);matching={'good_edges':[sorted(x) for x in good],'bad_edges':[sorted(x) for x in bad],
          'maximum_good_matching':gm,'maximum_bad_matching':bm,'minimum_count':bounds[0],'maximum_count':bounds[1]}
        require(bounds==(min(v[1] for v in records.values()),max(v[1] for v in records.values())),'Matching and enumeration disagree')
    else:matching=None
    return records,matching

def max_min_boundary(fields):
    """Max-min Floyd closure gives maximum bottleneck of any closed walk.

    The optimum closed-walk bottleneck equals that of a simple directed cycle.
    Zero initial diagonal represents absence of a zero-edge cycle.
    """
    ids=sorted(next(iter(fields.values())));n=len(ids)
    W=[[F(0) if i==j else max(f[ids[j]][0]/f[ids[i]][1] for f in fields.values()) for j in range(n)] for i in range(n)]
    D=[row[:] for row in W]
    for k in range(n):
        for i in range(n):
            for j in range(n):D[i][j]=max(D[i][j],min(D[i][k],D[k][j]))
    t=max([F(1)]+[D[i][i] for i in range(n)])
    return t,ids,W,[D[i][i] for i in range(n)]

def sqrt_bracket_check(q,b):
    lo,hi=map(F,b)
    require(0<=lo<hi and lo*lo<=q<hi*hi,'Invalid square-root bracket')
    require(hi-lo==F(1,10**12),'Unexpected display precision')

def synthetic():
    # Every 3x3 binary graph checks the matching implementation independently.
    for mask in range(512):
        adj=[{j for j in range(3) if mask & (1<<(3*i+j))} for i in range(3)]
        size,_=maximum_matching(adj)
        brute=max(sum(p[i] in adj[i] for i in range(3)) for p in itertools.permutations(range(3)))
        require(size==brute,'Synthetic bipartite matching')
    rng=random.Random(4268)
    for n in range(2,6):
        for rep in range(12):
            fields={c:{str(i):tuple(sorted([F(rng.randint(1,9)),F(rng.randint(1,9))])) for i in range(n)} for c in ['A','B']}
            t,ids,_,_=max_min_boundary(fields)
            brute=min(max([F(1)]+[field[p[i]][0]/field[p[j]][1] for field in fields.values() for i in range(n) for j in range(i,n)]) for p in itertools.permutations(ids))
            require(t==brute,'Synthetic maximum-cycle closure')
    cases=[([[F(1)]*3,[F(2)]*3],(3,3)),([[F(1),F(2),F(4)]]*2,(0,3)),([[F(1),F(4),F(16)],[F(2),F(8),F(32)]],(3,3))]
    for vals,expected in cases:
        rec,match=assignment_oracle(vals);require((min(x[1] for x in rec.values()),max(x[1] for x in rec.values()))==expected,'Known assignment case')
    # Touching intervals have no strict constraint at equality.
    fields={'A':{'x':(F(1),F(1)),'y':(F(2),F(2))},'B':{'x':(F(2),F(2)),'y':(F(1),F(1))}}
    t,_,W,_=max_min_boundary(fields);require(t==2 and all(v<=t for row in W for v in row),'Equality removes cycle edges')
    return {'matching_graphs':512,'boundary_total_order_comparisons':48,'assignment_hand_cases':3,'touching_case':1}

def review(b1,out):
    source_by={r['source_row']:r for r in b1['source_rows']};sc={(r['state'],r['cell']):r for r in b1['state_cells']}
    require(len(source_by)==396 and len(sc)==132,'Primary universe')
    require(all(r['complete'] and r['usable_rows']==3 for r in sc.values()),'Finite formulas outside complete domain')
    meta={}
    for r in source_by.values():
        key=r['state'];e=r['exon']
        require(key not in meta or meta[key]==e,'Conflicting state exon')
        meta[key]=e
        require(F(r['ratio'])==F(r['RNA'])/F(r['DNA']),'Preserved within-row ratio')
    original=b1['sensitivities'];require([F(s['gamma']) for s in original]==[F(1),F(5,4),F(3,2),F(2)],'Grid changed')
    groups={(g['intended'],g['cell']):g for g in original[0]['group_cell_results']}
    proposed={(g['intended'],g['cell']):g for g in out['groups']}
    require(len(groups)==60 and set(proposed)==set(groups),'Group coverage')
    profile_records=[];assignment_count=0;two=0;three=0;near=[];margin_records=[];grid_records=0
    for key in sorted(groups):
        base=groups[key];saved=proposed[key];markers=base['markers'];cell=key[1]
        require(markers==saved['markers'] and len(markers) in [2,3],'Marker identity/order')
        rs=[sc[m,cell] for m in markers];values=[]
        for r in rs:
            rows=r['source_rows'];require(len(set(rows))==3,'Distinct preparation row identities')
            require(all(source_by[i]['state']==r['state'] and source_by[i]['cell']==cell for i in rows),'Ratio row join')
            q=[F(source_by[i]['ratio']) for i in rows];values.append(q)
            require(min(q)==F(r['lower']) and max(q)==F(r['upper']),'Source envelope reconstruction')
        lower=[min(v) for v in values];upper=[max(v) for v in values]
        A=max(lower)/min(upper);B=max(upper[i]/lower[j] for i in range(len(markers)) for j in range(len(markers)) if i!=j)
        require(F(saved['A_exact'])==A and F(saved['B_exact'])==B,'A/B mismatch')
        require(F(saved['equality_feasible_from_t'])==max(F(1),A),'Equality threshold mismatch')
        require(saved['minimum_at_least_twofold_t_interval']==(None if A<2 else ['1',str(A/2)]),'Lower twofold transition')
        require(saved['maximum_at_most_twofold_t_interval']==(None if B>2 else ['1',str(2/B)]),'Upper twofold transition')
        sqrt_bracket_check(max(F(1),A),saved['equality_gamma_bracket'])
        require(saved['source_rows']==[r['source_rows'] for r in rs],'Preparation identities changed')
        records,matching=assignment_oracle(values);two+=len(markers)==2;three+=len(markers)==3
        keys=set();saved_stats=[]
        for item in saved['assignments']:
            perms=tuple(tuple(x) for x in item['permutations'])
            require(perms in records and perms not in keys,'Assignment absent or duplicate')
            keys.add(perms);f,count=records[perms]
            require(tuple(map(F,item['folds']))==f and item['preparations_at_least_twofold']==count,'Assignment exact fold/count mismatch')
            saved_stats.append((f,count))
        require(keys==set(records),'Missing assignments')
        require(saved['assignment_count']==len(records),'Assignment denominator')
        assignment_count+=len(records);allfolds=[f for fs,count in records.values() for f in fs]
        mincount=min(count for fs,count in records.values());maxcount=max(count for fs,count in records.values())
        expected={'minimum_preparations_at_least_twofold':mincount,'maximum_preparations_at_least_twofold':maxcount,
          'worst_minimum_fold':min(min(fs) for fs,count in records.values()),'best_minimum_fold':max(min(fs) for fs,count in records.values()),
          'minimum_maximum_fold':min(max(fs) for fs,count in records.values()),'maximum_maximum_fold':max(max(fs) for fs,count in records.values())}
        for name,val in expected.items():require(F(saved[name])==val,'Assignment summary mismatch '+name)
        distances=sorted({f-2 for f in allfolds if abs(f-2)<=F(1,10**12)})
        if distances:near.append({'intended':key[0],'cell':cell,'exact_fold_minus_two':list(map(str,distances)),
           'minimum_preparations':mincount,'maximum_preparations':maxcount,'has_exact_tie':F(0) in distances})
        if abs(A-2)<=F(1,10**12):margin_records.append({'intended':key[0],'cell':cell,'quantity':'A','exact_minus_two':str(A-2),'exact_at_least_two':A>=2})
        if mincount==3:margin_records.append({'intended':key[0],'cell':cell,'quantity':'smallest_recorded_assignment_fold_minus_two','exact_minus_two':str(min(allfolds)-2),'near_display_boundary':min(allfolds)-2<=F(1,10**12)})
        for s in original:
            g=next(g for g in s['group_cell_results'] if (g['intended'],g['cell'])==key);t=F(s['gamma'])**2
            lo=max(F(1),A/t);hi=B*t;point=max(sorted(v)[1] for v in values)/min(sorted(v)[1] for v in values)
            require(F(g['minimum_fold'])==lo and F(g['maximum_fold'])==hi,'Original grid profile disagreement')
            require(g['common_response_feasible']==(lo==1) and g['complete'] is True and g['maximum_attained'] is True,'Original grid logical fields')
            require(F(g['point_fold'])==point,'Original point baseline')
            require(g['markers']==markers and g['marker_count']==len(markers),'Original fixed ascertainment')
            for th,claims in g['thresholds'].items():
                require(claims=={'minimum_at_least':lo>=F(th),'maximum_at_most':hi<=F(th)},'Original threshold grid')
            grid_records+=1
        profile_records.append({'intended':key[0],'cell':cell,'markers':markers,'A':str(A),'B':str(B),
          'minimum_count':mincount,'maximum_count':maxcount,'matching_certificate':matching,'assignment_count':len(records)})
    require((two,three,assignment_count,grid_records)==(48,12,720,240),'Complete record accounting')
    exon_results=[];original_status_checks=0
    for saved in out['exons']:
        e=saved['exon'];fields={c:{v:(F(sc[v,c]['lower']),F(sc[v,c]['upper'])) for v in sorted(meta) if meta[v]==e} for c in ['K562','KBM7']}
        t,ids,W,diag=max_min_boundary(fields);n=len(ids);idx={v:i for i,v in enumerate(ids)}
        require(F(saved['t_exact'])==t,'Independent maximum-cycle boundary mismatch')
        sqrt_bracket_check(t,saved['gamma_bracket'])
        order=saved['feasible_order'];require(len(order)==n and set(order)==set(ids),'Boundary order coverage')
        position={v:i for i,v in enumerate(order)}
        required={(ids[i],ids[j]):W[i][j] for i in range(n) for j in range(n) if i!=j and W[i][j]>t}
        archived={(r['before'],r['after']):F(r['weight']) for r in saved['edges_at_boundary']}
        require(required==archived,'Edges at equality boundary')
        require(all(position[a]<position[b] for a,b in required),'Order violates strict boundary edge')
        require(set(saved['scaled_response_witnesses'])==set(fields),'Witness cell coverage')
        for cell,d in fields.items():
            witness=saved['scaled_response_witnesses'][cell];require(set(witness)==set(ids),'Witness state coverage')
            largest=F(0)
            for v in order:
                largest=max(largest,d[v][0]);y=F(witness[v])
                require(y==largest and d[v][0]<=y<=t*d[v][1] and y>0,'Boundary rational monotone completion')
        cert=saved['lower_certificate']
        if t==1:require(cert is None,'Spurious below-domain certificate')
        else:
            require(bool(cert),'Missing below-boundary cycle')
            require(all(cert[i]['after']==cert[(i+1)%len(cert)]['before'] for i in range(len(cert))),'Cycle not closed')
            require(len({r['before'] for r in cert})==len(cert),'Repeated cycle vertices')
            weights=[]
            for edge in cert:
                a,b=edge['before'],edge['after'];w=W[idx[a]][idx[b]];weights.append(w)
                origins=sorted(c for c,d in fields.items() if d[b][0]/d[a][1]==w)
                require(F(edge['weight'])==w and edge['contexts']==origins,'Cycle inequality/origin')
            require(min(weights)==t,'Cycle does not certify every t below boundary')
        for s in original:
            previous=next(z for z in s['exon_shared_order_results'] if z['exon']==e)
            require(previous['status']==('feasible' if F(s['gamma'])**2>=t else 'infeasible'),'Original common-order grid disagreement')
            original_status_checks+=1
        exon_results.append({'exon':e,'state_count':n,'t_exact':str(t),'gamma_bracket':saved['gamma_bracket'],'closure_diagonal':list(map(str,diag)),'cycle_edges':0 if cert is None else len(cert),'boundary_certificate_passed':True})
    require(sorted(x['exon'] for x in exon_results)==[17,25,27,39,42] and original_status_checks==20,'Exon/status coverage')
    summary={'group_cells':len(profile_records),'enumerated_group_cell_assignments':assignment_count,
      'all_three_preparations_twofold_under_every_assignment':sum(r['minimum_count']==3 for r in profile_records),
      'at_least_one_preparation_twofold_under_every_assignment':sum(r['minimum_count']>=1 for r in profile_records),
      'assignment_changes_twofold_preparation_count':sum(r['minimum_count']!=r['maximum_count'] for r in profile_records),
      'no_assignment_has_twofold_preparation':sum(r['maximum_count']==0 for r in profile_records)}
    for k,v in summary.items():require(out['summary'][k]==v,'Headline mismatch '+k)
    require(out['summary']['exon_boundaries']==[{k:z[k] for k in ['exon','t_exact','gamma_bracket']} for z in exon_results],'Boundary summary mismatch')
    return {'status':'PASS','summaries':summary,'two_marker_cases':two,'three_marker_cases':three,
      'assignment_records_checked':assignment_count,'exact_assignment_fold_values_checked':3*assignment_count,
      'original_group_grid_records_checked':grid_records,'original_order_statuses_checked':original_status_checks,
      'profile_records':profile_records,'exon_results':exon_results,'decimal_boundary_diagnostics':{
          'diagnostic_only_absolute_nearness_to_two':'1/1000000000000','not_a_new_decision_rule':True,
          'cases_with_assignment_folds_near_two':len(near),'case_records':near,'transition_and_all_three_margins':margin_records},
      'maximum_exact_numerical_discrepancy':'0','limitations':['Finite stored decimal ratios only.','Anonymous assignments are not authenticated replicate labels.','No calibrated uncertainty or independent-assay replication.','Point-set assignment robustness differs from continuous-envelope inference.']}

def main():
    if not __debug__:raise RuntimeError('Optimized execution prohibited')
    ap=argparse.ArgumentParser();ap.add_argument('--b1',type=Path,required=True);ap.add_argument('--result',type=Path,required=True)
    ap.add_argument('--protocol',type=Path,required=True);ap.add_argument('--derivation',type=Path,required=True);ap.add_argument('--analysis',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    inputs={k:getattr(a,k) for k in ['b1','result','protocol','derivation','analysis']};pins={k:sha(p) for k,p in inputs.items()}
    require(pins['b1']==B1_SHA and pins['result']==RESULT_SHA,'Held source mismatch')
    a.out.mkdir(parents=True);started=time.perf_counter()
    try:
        syn=synthetic();b1=json.loads(a.b1.read_text());result=json.loads(a.result.read_text())
        review_result=review(b1,result)
        require(pins=={k:sha(p) for k,p in inputs.items()},'Reviewed source changed')
        put(a.out/'reconstruction.json',review_result)
        receipt={'status':'PASS','checks':N,'input_sha256':pins,'reviewer_source_sha256':sha(Path(__file__)),
          'reconstruction_sha256':sha(a.out/'reconstruction.json'),'synthetic_checks':syn,'elapsed_seconds':time.perf_counter()-started,
          'python':platform.python_version(),'maximum_exact_discrepancy':'0','new_source_outcomes_read':False,
          'review_identity':'Separately authored computational review with inherited project context; not blinded or external.'}
        put(a.out/'receipt.json',receipt);print(json.dumps({'receipt':receipt,'summary':review_result['summaries'],
          'near_two_cases':review_result['decimal_boundary_diagnostics']['cases_with_assignment_folds_near_two'],'boundary_margins':review_result['decimal_boundary_diagnostics']['transition_and_all_three_margins']}))
    except Exception as exc:
        put(a.out/'failure.json',{'status':'BLOCKED','checks_before_failure':N,'error':str(exc),'type':type(exc).__name__,'input_sha256':pins});raise

if __name__=='__main__':main()
