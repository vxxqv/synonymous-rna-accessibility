"""Separate simple-path oracle; synthetic inputs only, no biological file reader."""
import sys
if not __debug__:raise RuntimeError('Optimized Python is not allowed')
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
import datetime
import hashlib
import json
import math
import platform
import time
import argparse
import scenario_geometry as impl

def need(ok,msg):
    if not ok:raise AssertionError(msg)

def paths(n,edges,s,t):
    adj=[set() for _ in range(n)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    found=[]
    def visit(path):
        if path[-1]==t:found.append(tuple(path));return
        for v in sorted(adj[path[-1]]):
            if v not in path:visit(path+[v])
    visit([s]);return found

def oracle(fields,allpaths,s,t):
    # Direct path/scenario incidence. No Floyd closure or author helper is used.
    bottlenecks=[[min(row[v] for v in path) for path in allpaths] for row in fields]
    caps=tuple(max(row,default=-math.inf) for row in bottlenecks)
    fixed=max((min(row[k] for row in bottlenecks) for k in range(len(allpaths))),default=-math.inf)
    upper=max((min(max(row[v] for row in fields) for v in path) for path in allpaths),default=-math.inf)
    return {'E':min(row[v] for row in fields for v in (s,t)),
            'C_fix':fixed,'C_ad':min(caps),'C_max':max(caps),'C_upper':upper,'scenario_capacities':caps}

def oracle_label(fields,allpaths,s,t,tau):
    if any(row[s]<tau or row[t]<tau for row in fields):return 'ineligible_endpoints'
    accepted=[[all(row[v]>=tau for v in path) for path in allpaths] for row in fields]
    if any(all(row[k] for row in accepted) for k in range(len(allpaths))):return 'fixed_route'
    if all(any(row) for row in accepted):return 'adaptive_routes_only'
    if any(any(row) for row in accepted):return 'mixed_scenario_connectivity'
    if any(all(any(row[v]>=tau for row in fields) for v in path) for path in allpaths):return 'joint_obstruction_box_connected'
    return 'upper_box_obstruction'

def run():
    counts={'complete_landscapes':0,'complete_pair_metric_checks':0,'threshold_classifications':0,
            'unknown_landscapes':0,'completion_metric_bound_checks':0,'completion_certificates':0}
    known=[]
    examples=[('parallel',[(0,1),(1,3),(0,2),(2,3)],[[2,1,0,2],[2,0,1,2]],(0,1,1,1,2)),
              ('series',[(0,1),(1,2),(2,3)],[[2,1,0,2],[2,0,1,2]],(0,0,0,1,2)),
              ('mixed',[(0,1),(1,3),(0,2),(2,3)],[[2,1,0,2],[2,Q(1,2),0,2]],(Q(1,2),Q(1,2),1,1,2))]
    for name,edges,fields,expected in examples:
        a=impl.analyze(fields,edges);rec=next(x for x in a['pairs'] if (x['s'],x['t'])==(0,3))
        need(tuple(rec[k] for k in ('C_fix','C_ad','C_max','C_upper','E'))==expected,'Known '+name)
        known.append({'name':name,'record':impl.json_ready(rec)})
    pairgraph=((0,1),(0,2),(1,2))
    allfields=list(product((0,1,2),repeat=3))
    thresholds=[Q(i,2) for i in range(6)]
    for mask in range(8):
        edges=[e for i,e in enumerate(pairgraph) if mask>>i&1]
        pathmap={(s,t):paths(3,edges,s,t) for s,t in pairgraph}
        for fields in product(allfields,repeat=2):
            a=impl.analyze(fields,iter(edges));counts['complete_landscapes']+=1
            need(a['pair_count']==3,'Three-vertex pair coverage')
            for rec in a['pairs']:
                s,t=rec['s'],rec['t'];ps=pathmap[s,t];ref=oracle(fields,ps,s,t)
                need(all(rec[key]==val for key,val in ref.items()),'Exact path oracle')
                counts['complete_pair_metric_checks']+=1
                for tau in thresholds:
                    expected=oracle_label(fields,ps,s,t,tau)
                    need(impl.classify_exact(rec,tau)==expected,'Threshold quantifier classification')
                    intervals=[entry for entry in rec['intervals'] if entry['nonempty'] and
                               (tau>=entry['lower'] if entry['lower_closed'] else tau>entry['lower']) and tau<=entry['upper']]
                    need(len(intervals)==(0 if expected=='ineligible_endpoints' else 1),'Interval partition')
                    if intervals:need(intervals[0]['label']==expected,'Interval label')
                    counts['threshold_classifications']+=1
    edges=[(0,1),(1,2)];pathmap={(s,t):paths(3,edges,s,t) for s,t in pairgraph}
    for flat in product((0,1,None),repeat=6):
        missing=[i for i,x in enumerate(flat) if x is None]
        if not missing:continue
        fields=[flat[:3],flat[3:]];a=impl.analyze(fields,edges);counts['unknown_landscapes']+=1
        need(a['unknown_cells']==[(i//3,i%3) for i in missing],'Unknown indexing')
        for completion in product((0,1),repeat=len(missing)):
            filled=list(flat)
            for i,x in zip(missing,completion):filled[i]=x
            fs=[filled[:3],filled[3:]]
            for rec in a['pairs']:
                s,t=rec['s'],rec['t'];ps=pathmap[s,t];ref=oracle(fs,ps,s,t)
                need(all(lo<=ref[key]<=hi for key,(lo,hi) in rec['bounds'].items()),'Unknown metric soundness')
                counts['completion_metric_bound_checks']+=1
                for tau in (Q(0),Q(1,2),Q(1),Q(3,2)):
                    label=impl.classify_bounded(rec,tau)
                    if not label.startswith('unresolved'):
                        need(label==oracle_label(fs,ps,s,t,tau),'False unknown certificate')
                        counts['completion_certificates']+=1
    # Permute every vertex of one informative diamond and reverse scenarios.
    import itertools
    fields=examples[0][2];edges=examples[0][1]
    base=impl.analyze(fields,edges)
    for perm in itertools.permutations(range(4)):
        fs=[[row[old] for old in perm] for row in reversed(fields)]
        inverse={old:new for new,old in enumerate(perm)}
        es=[(inverse[a],inverse[b]) for a,b in edges]
        a=impl.analyze(fs,es)
        for rec in a['pairs']:
            oldpair=tuple(sorted((perm[rec['s']],perm[rec['t']])))
            ref=next(x for x in base['pairs'] if (x['s'],x['t'])==oldpair)
            need(all(rec[k]==ref[k] for k in impl.METRICS),'Permutation invariance')
    # Full application shape, using only an abstract product and constant scores.
    states=list(product(range(3),range(3),range(2),range(2),range(2)))
    es=[(i,j) for i,j in combinations(range(72),2) if sum(abs(a-b) for a,b in zip(states[i],states[j]))==1]
    need(len(es)==204,'Abstract product graph size')
    fields=[[r%3]*72 for r in range(27)]
    a=impl.analyze(fields,es)
    need(a['pair_count']==2556 and a['scenarios']==27 and not a['unknown_cells'],'Full shape coverage')
    for rec in a['pairs']:
        need(rec['scenario_capacities']==tuple(Q(r%3) for r in range(27)),'All27 scenario order')
        need(tuple(rec[k] for k in impl.METRICS)==(0,0,0,2,2),'Constant known metrics')
    invalid=[([],[]),([[]],[]),([[0],[0,1]],[]),([[float('nan')]],[]),([[-1]],[]),([[True]],[]),
             ([[0,1]],[(0,0)]),([[0,1]],[(0,2)]),([[0,1]],[(0,1),(1,0)])]
    for fields,edges in invalid:
        try:impl.analyze(fields,edges)
        except (ValueError,TypeError):pass
        else:raise AssertionError('Invalid input accepted')
    need(impl.analyze([[0]],[])['pair_count']==0,'Singleton empty pair family')
    return {'status':'PASS','assertion_tested':'Exact path quantifiers, complete interval partition, and conservative unknown-completion bounds on the declared synthetic domain.',
            'counts':counts,'known_examples':known,'vertex_permutations':24,
            'full_shape':{'vertices':72,'edges':204,'scenarios':27,'pairs':2556,'scenario_pair_entries':69012},
            'invalid_inputs_rejected':len(invalid),'observed_data_opened':False,
            'independence':'Distinct simple-path oracle by the same author; not a separately authored or blinded review.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    need(not args.out.exists(),'Refusing existing output')
    here=Path(__file__).resolve().parent
    sources=[here/'protocol.md',here/'derivation.md',here/'scenario_geometry.py',Path(__file__).resolve()]
    hashes={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    started=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.perf_counter()
    result=run();result['runtime_seconds']=time.perf_counter()-t
    result['started_at']=started;result['python']=sys.version;result['platform']=platform.platform()
    result['source_hashes']=hashes
    need(all(hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected for path,expected in hashes.items()),'Inputs changed during run')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','counts','runtime_seconds','observed_data_opened')}))

if __name__=='__main__':main()
