"""Portable exact branch-D reconstruction. Standard library, no production import.

Default inputs are immutable repository-relative gzip JSON results. Explicit
--primary/--result/--freeze paths may instead be supplied. No Git/network needed.
"""
import sys
if not __debug__:
    raise RuntimeError('Optimized execution refused')
from collections import deque
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse,csv,gzip,hashlib,json,math

PRIMARY='3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'
RESULT='0c2ad1a6a602f80e44ff05eaf3321522d9eaf0940382c599309a04d33f4b383e'
FREEZE='340624e795e600abbbaeb3daefa89b054ceca7e356855014d320df20b043a10f'
LABELS=('fixed_route','adaptive_routes_only','mixed_scenario_connectivity',
        'joint_obstruction_box_connected','upper_box_obstruction')


def need(ok,message):
    if not ok:raise AssertionError(message)


def sha(data):return hashlib.sha256(data).hexdigest()


def load(path,pin):
    original=path.read_bytes()
    raw=gzip.decompress(original) if original[:2]==b'\x1f\x8b' else original
    need(sha(raw)==pin,'Input decompressed hash mismatch: '+path.name)
    return json.loads(raw),{'stored_sha256':sha(original),'decompressed_sha256':sha(raw)}


def closure(scores,adj):
    """Independent BFS at every distinct level, not a Floyd max-min closure."""
    n=len(scores);out=[[None]*n for _ in range(n)]
    for threshold in sorted(set(scores),reverse=True):
        unseen={v for v in range(n) if scores[v]>=threshold}
        while unseen:
            start=min(unseen);component={start};queue=deque([start])
            while queue:
                for v in adj[queue.popleft()]:
                    if scores[v]>=threshold and v not in component:
                        component.add(v);queue.append(v)
            unseen-=component
            for s in component:
                for t in component:
                    if out[s][t] is None:out[s][t]=threshold
    return [[-math.inf if v is None else v for v in row] for row in out]


def cell(label,cut,ceiling,E):
    lower=max(Q(0),cut);upper=min(E,ceiling)
    lower_closed=cut<Q(0)
    nonempty=upper>lower or (upper==lower and lower_closed)
    return {'label':label,'lower':str(lower),'upper':str(upper),
            'lower_closed':lower_closed,'upper_closed':True,'nonempty':nonempty,
            'width':str(max(Q(0),upper-lower))}


def member(c,tau):
    return c['nonempty'] and (tau>=Q(c['lower']) if c['lower_closed'] else tau>Q(c['lower'])) and tau<=Q(c['upper'])


def rebuild(primary,result):
    rows=primary['vertices_data'];n=len(rows)
    need(n==72 and [r['vertex'] for r in rows]==list(range(n)),'Primary vertex inventory')
    sequences=[r['sequence'] for r in rows]
    need(len(set(sequences))==72 and all(len(s)==63 for s in sequences),'Sequence dimensions')
    need(all(len(r['scenario_ratios'])==27 for r in rows),'Scenario coverage')
    need(all(v is not None for r in rows for v in r['scenario_ratios']),'Expected complete source, no silent unknown exclusion')
    fields=[[Q(r['scenario_ratios'][k]) for r in rows] for k in range(27)]
    need(all(v>0 for f in fields for v in f),'Expected positive complete fields')
    edges=[(s,t) for s,t in combinations(range(n),2) if sum(a!=b for a,b in zip(sequences[s],sequences[t]))==1]
    need(len(edges)==204,'Graph edge count')
    adj=[set() for _ in range(n)]
    for s,t in edges:adj[s].add(t);adj[t].add(s)
    lower=[min(f[v] for f in fields) for v in range(n)]
    upper=[max(f[v] for f in fields) for v in range(n)]
    need([Q(r['lower']) for r in rows]==lower and [Q(r['upper']) for r in rows]==upper,'Unchanged envelope arithmetic')
    matrices=[closure(f,adj) for f in fields]
    low=closure(lower,adj);high=closure(upper,adj)
    pairs=list(combinations(range(n),2))
    need(len(result['pairs'])==len(primary['bounded_pairs'])==2556,'All endpoint pair records required')
    for key,expected in [('vertices',72),('scenarios',27),('edges',204),('pair_count',2556),('unknown_cells',[])]:
        need(result[key]==expected,'Result metadata '+key)
    need(result['primary_result_sha256']==PRIMARY,'Result source identity')
    counts={label:0 for label in LABELS};adaptive=set();mixed=set();nonfull=set();full=set()
    endpoint_cases=0;critical_memberships=0;records=[]
    for index,(s,t) in enumerate(pairs):
        caps=[matrix[s][t] for matrix in matrices]
        E=min(lower[s],lower[t]);fix=low[s][t];ad=min(caps);mx=max(caps);up=high[s][t]
        need(fix<=ad<=mx<=up and ad<=E,'Capacity ordering')
        cuts=[-math.inf,fix,ad,mx,up,math.inf]
        cells=[cell(label,cuts[k],cuts[k+1],E) for k,label in enumerate(LABELS)]
        expected={'s':s,'t':t,'scenario_capacities':list(map(str,caps)),
                  'E':str(E),'C_fix':str(fix),'C_ad':str(ad),'C_max':str(mx),'C_upper':str(up),
                  'status':'exact','intervals':cells,
                  'all_scenario_obstruction':cell('all_scenario_obstruction',mx,math.inf,E),
                  'upper_box_obstruction':cell('upper_box_obstruction',up,math.inf,E)}
        need(result['pairs'][index]==expected,'Complete pair record mismatch at '+str((s,t)))
        old=primary['bounded_pairs'][index]
        need((old['a'],old['b'])==(s,t) and Q(old['endpoint_lower'])==E and Q(old['bottleneck_lower'])==fix and Q(old['bottleneck_upper'])==up,'Primary pair bound mismatch')
        for k,cap in enumerate(caps):
            label=f'BR{k//3+1}_TR{k%3+1}'
            prior=primary['all_scenario_pairs'][label][index]
            need(Q(prior['bottleneck_lower'])==cap==Q(prior['bottleneck_upper']),'Original independently verified scenario capacity drift')
        for c in cells:
            counts[c['label']]+=c['nonempty']
        if cells[1]['nonempty']:adaptive.add((s,t))
        if cells[2]['nonempty']:mixed.add((s,t))
        (full if fix==E else nonfull).add((s,t))
        # Exact tie and between-cut membership verification over the entire
        # finite critical decomposition, not a user-selected numerical grid.
        points=sorted({Q(0),E,*[c for c in (fix,ad,mx,up) if 0<=c<=E]})
        tests=points+[(a+b)/2 for a,b in zip(points,points[1:])]+[E+1]
        for tau in tests:
            wanted='ineligible_endpoints' if tau>E else (
                LABELS[0] if tau<=fix else LABELS[1] if tau<=ad else
                LABELS[2] if tau<=mx else LABELS[3] if tau<=up else LABELS[4])
            actual=[c['label'] for c in cells if member(c,tau)]
            need(actual==([] if wanted=='ineligible_endpoints' else [wanted]),'Threshold tie/partition mismatch')
            critical_memberships+=1
        need(fix>0 and member(cells[0],Q(0)),'tau0 or positive fixed-route range')
        endpoint_cases+=1
        records.append(expected)
    need(counts==result['summary']['category_nonempty_pair_counts'],'Summary category counts')
    all_block=sum(r['all_scenario_obstruction']['nonempty'] for r in records)
    box_block=sum(r['upper_box_obstruction']['nonempty'] for r in records)
    need(all_block==result['summary']['all_scenario_obstruction_pairs'] and box_block==result['summary']['upper_box_obstruction_pairs'],'Obstruction summary counts')
    need(adaptive|mixed==nonfull,'Full-range fixed route equivalence')
    return {'status':'PASS_ALL_EXACT_RECORDS','pair_count':2556,'scenario_capacity_count':69012,
            'complete_record_comparisons':2556,'critical_membership_checks':critical_memberships,
            'category_nonempty_pair_counts':counts,'all_scenario_obstruction_pairs':all_block,
            'upper_box_obstruction_pairs':box_block,'adaptive_mixed_overlap':len(adaptive&mixed),
            'adaptive_mixed_union':len(adaptive|mixed),'C_fix_equals_E_pairs':len(full),
            'C_fix_less_than_E_pairs':len(nonfull),'tau0_fixed_route_pairs':endpoint_cases,
            'positive_width_fixed_route_pairs':sum(Q(r['intervals'][0]['width'])>0 for r in records),
            'all27_fields_retained':True,'unknown_vertices_dropped':False,
            'qualification':'Exact descriptive quantifiers on27 dependent preparation fields. No confidence calibration, causal mechanism, measured evolutionary route or general full-synonymous-space conclusion.'},records


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--primary',type=Path);p.add_argument('--result',type=Path);p.add_argument('--freeze',type=Path)
    args=p.parse_args();repo=args.repo.resolve();out=args.out.resolve()
    need(repo.is_dir(),'Repository directory missing')
    need(not out.exists(),'Refusing existing output')
    need(not out.is_relative_to(repo),'Output must be outside repository')
    freeze_path=args.freeze or repo/'analysis/cycle7/freeze.json'
    primary_path=args.primary or repo/'results/cycle6/full_result_v1.json.gz'
    result_path=args.result or repo/'results/cycle7/scenario_geometry.json.gz'
    freeze,freeze_hash=load(freeze_path,FREEZE)
    frozen={name:sha((freeze_path.parent/name).read_bytes())==expected for name,expected in freeze['files'].items()}
    need(all(frozen.values()),'Frozen source changed')
    need(freeze['primary_result_sha256']==PRIMARY and freeze['primary_rules_changed'] is False,'Freeze scope changed')
    primary,primary_hash=load(primary_path,PRIMARY);result,result_hash=load(result_path,RESULT)
    audit,records=rebuild(primary,result)
    audit['inputs']={'primary':primary_hash,'result':result_hash,'freeze':freeze_hash}
    audit['frozen_source_checks']=frozen;audit['verifier_sha256']=sha(Path(__file__).read_bytes())
    audit['historical_git_verified']=False
    out.mkdir(parents=True)
    with (out/'all_pairs.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['s','t','E','C_fix','C_ad','C_max','C_upper',*[f'BR{k//3+1}_TR{k%3+1}' for k in range(27)],*LABELS])
        for r in records:
            writer.writerow([r[k] for k in ('s','t','E','C_fix','C_ad','C_max','C_upper')]+r['scenario_capacities']+[c['nonempty'] for c in r['intervals']])
    audit['all_pairs_csv_sha256']=sha((out/'all_pairs.csv').read_bytes())
    (out/'verification.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in ('inputs','frozen_source_checks')}))


if __name__=='__main__':main()
