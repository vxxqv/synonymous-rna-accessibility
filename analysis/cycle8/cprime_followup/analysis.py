from fractions import Fraction as Q
from itertools import permutations,product
from collections import defaultdict
from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,csv,random
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'cprime_results/primary_v1/result.json'
SOURCE_SHA='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'

def topology(nodes,edges):
    adj={x:[] for x in nodes};deg={x:0 for x in nodes}
    for a,b in edges:adj[a].append(b);deg[b]+=1
    ready=sorted(x for x in nodes if not deg[x]);order=[]
    while ready:
        x=ready.pop(0);order.append(x)
        for y in sorted(adj[x]):
            deg[y]-=1
            if not deg[y]:ready.append(y);ready.sort()
    return order if len(order)==len(nodes) else None

def cycle(nodes,edges):
    adj={x:[] for x in nodes}
    for a,b in edges:adj[a].append(b)
    color={x:0 for x in nodes};stack=[]
    def dfs(x):
        color[x]=1;stack.append(x)
        for y in sorted(adj[x]):
            if color[y]==1:return stack[stack.index(y):]+[y]
            if not color[y]:
                found=dfs(y)
                if found:return found
        stack.pop();color[x]=2
        return None
    for x in nodes:
        if not color[x]:
            found=dfs(x)
            if found:return found
    return None

def weights(intervals):
    nodes=sorted(next(iter(intervals.values())));w={};origins={}
    assert all(set(d)==set(nodes) for d in intervals.values())
    for a,b in permutations(nodes,2):
        vals={c:d[b][0]/d[a][1] for c,d in intervals.items()}
        w[a,b]=max(vals.values());origins[a,b]=sorted(c for c,v in vals.items() if v==w[a,b])
    return nodes,w,origins

def completion(order,intervals,t):
    result={}
    for c,d in intervals.items():
        v=Q(0);res={}
        for x in order:
            v=max(v,d[x][0]);assert 0<v<=t*d[x][1];res[x]=str(v)
        result[c]=res
    return result

def boundary(intervals):
    nodes,w,origins=weights(intervals);values=sorted({Q(1)}|{v for v in w.values() if v>1})
    last=Q(1);lastcycle=None
    for t in values:
        edges=[e for e,v in w.items() if v>t];order=topology(nodes,edges)
        if order is not None:
            if t>1:
                before=[e for e,v in w.items() if v>=t]
                cyc=cycle(nodes,before);assert cyc and min(w[a,b] for a,b in zip(cyc,cyc[1:]))==t
            else:cyc=None
            return dict(t_exact=str(t),gamma_bracket=root_bracket(t),
                lower_certificate=None if cyc is None else [dict(before=a,after=b,weight=str(w[a,b]),contexts=origins[a,b]) for a,b in zip(cyc,cyc[1:])],
                feasible_order=order,scaled_response_witnesses=completion(order,intervals,t),
                edges_at_boundary=[dict(before=a,after=b,weight=str(w[a,b])) for a,b in edges])
    raise AssertionError('Finite graph has no boundary')

def root_bracket(q,places=12):
    scale=10**places
    # Integer comparisons certify the bracket; no floating-point decision.
    from math import isqrt
    k=isqrt((q.numerator*scale*scale)//q.denominator)
    lo=Q(k,scale);hi=Q(k+1,scale)
    assert lo*lo<=q<hi*hi
    return [format(Decimal(k)/Decimal(scale),'f'),format(Decimal(k+1)/Decimal(scale),'f')]

def assignment_profile(values):
    assert len(values) in (2,3) and all(len(x)==3 for x in values)
    records=[]
    for rest in product(list(permutations(range(3))),repeat=len(values)-1):
        perm=[(0,1,2),*rest]
        folds=[max(values[m][perm[m][j]] for m in range(len(values)))/min(values[m][perm[m][j]] for m in range(len(values))) for j in range(3)]
        records.append(dict(permutations=[list(x) for x in perm],folds=[str(x) for x in folds],preparations_at_least_twofold=sum(x>=2 for x in folds)))
    counts=[r['preparations_at_least_twofold'] for r in records]
    return dict(assignments=records,assignment_count=len(records),minimum_preparations_at_least_twofold=min(counts),maximum_preparations_at_least_twofold=max(counts),
                worst_minimum_fold=str(min(min(map(Q,r['folds'])) for r in records)),
                best_minimum_fold=str(max(min(map(Q,r['folds'])) for r in records)),
                minimum_maximum_fold=str(min(max(map(Q,r['folds'])) for r in records)),
                maximum_maximum_fold=str(max(max(map(Q,r['folds'])) for r in records)))

def tests():
    rng=random.Random(20261025);checks=0
    for n in range(2,6):
        for rep in range(20):
            iv={str(c):{str(i):tuple(sorted([Q(rng.randint(1,8)),Q(rng.randint(1,8))])) for i in range(n)} for c in range(2)}
            ans=boundary(iv);t=Q(ans['t_exact']);nodes,w,_=weights(iv)
            # Independent exhaustive order criterion: max earlier L <= t*later U.
            need=[]
            for order in permutations(nodes):
                need.append(max([Q(1)]+[iv[c][order[i]][0]/iv[c][order[j]][1] for c in iv for i in range(n) for j in range(i,n)]))
            assert t==min(need);checks+=1
    assert assignment_profile([[Q(1)]*3,[Q(2)]*3])['minimum_preparations_at_least_twofold']==3
    v=assignment_profile([[Q(1),Q(2),Q(4)],[Q(1),Q(2),Q(4)]])
    assert v['minimum_preparations_at_least_twofold']==0 and v['maximum_preparations_at_least_twofold']==3
    for values in ([[Q(1),Q(3),Q(9)],[Q(2),Q(4),Q(8)]],[[Q(1),Q(3),Q(9)],[Q(2),Q(4),Q(8)],[Q(1),Q(2),Q(6)]]):
        a=assignment_profile(values);b=assignment_profile([list(reversed(x)) for x in values]);assert sorted((r['preparations_at_least_twofold'],tuple(sorted(map(Q,r['folds'])))) for r in a['assignments'])==sorted((r['preparations_at_least_twofold'],tuple(sorted(map(Q,r['folds'])))) for r in b['assignments']);checks+=1
    return {'total_order_instances':80,'assignment_hand_and_relabel_checks':4,'status':'PASS'}

def main():
    out=HERE/'results_v1'
    if out.exists():raise FileExistsError(out)
    test=tests()
    raw=SOURCE.read_bytes();assert hashlib.sha256(raw).hexdigest()==SOURCE_SHA
    d=json.loads(raw);sc={(r['state'],r['cell']):r for r in d['state_cells']};rawrows={r['source_row']:r for r in d['source_rows']}
    groups=[]
    for g in d['sensitivities'][0]['group_cell_results']:
        rr=[sc[s,g['cell']] for s in g['markers']];assert all(r['complete'] for r in rr)
        A=max(Q(r['lower']) for r in rr)/min(Q(r['upper']) for r in rr)
        B=max(Q(a['upper'])/Q(b['lower']) for a,b in permutations(rr,2))
        vals=[[Q(rawrows[i]['ratio']) for i in r['source_rows']] for r in rr]
        p=assignment_profile(vals)
        groups.append(dict(intended=g['intended'],cell=g['cell'],markers=g['markers'],source_rows=[r['source_rows'] for r in rr],
            A_exact=str(A),B_exact=str(B),equality_feasible_from_t=str(max(Q(1),A)),equality_gamma_bracket=root_bracket(max(Q(1),A)),
            minimum_at_least_twofold_t_interval=None if A<2 else ['1',str(A/2)],
            maximum_at_most_twofold_t_interval=None if B>2 else ['1',str(Q(2)/B)],**p))
    exons=[]
    state_exon={r['state']:r['exon'] for r in d['source_rows']}
    for exon in sorted(set(state_exon.values())):
        iv=defaultdict(dict)
        for r in d['state_cells']:
            if state_exon[r['state']]==exon:iv[r['cell']][r['state']]=(Q(r['lower']),Q(r['upper']))
        exons.append({'exon':exon,**boundary(dict(iv))})
    summary={'group_cells':len(groups),'enumerated_group_cell_assignments':sum(g['assignment_count'] for g in groups),
        'all_three_preparations_twofold_under_every_assignment':sum(g['minimum_preparations_at_least_twofold']==3 for g in groups),
        'at_least_one_preparation_twofold_under_every_assignment':sum(g['minimum_preparations_at_least_twofold']>=1 for g in groups),
        'assignment_changes_twofold_preparation_count':sum(g['minimum_preparations_at_least_twofold']!=g['maximum_preparations_at_least_twofold'] for g in groups),
        'no_assignment_has_twofold_preparation':sum(g['maximum_preparations_at_least_twofold']==0 for g in groups),
        'exon_boundaries':[{'exon':r['exon'],'t_exact':r['t_exact'],'gamma_bracket':r['gamma_bracket']} for r in exons]}
    result={'status':'EXPLORATORY_UNREVIEWED','input_sha256':SOURCE_SHA,'test_receipt':test,'groups':groups,'exons':exons,'summary':summary}
    out.mkdir();(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    flat=[{k:v for k,v in g.items() if k not in ['assignments','markers','source_rows']} for g in groups]
    with (out/'all_group_profiles.csv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    manifest={'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())},'source_sha256':SOURCE_SHA,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
