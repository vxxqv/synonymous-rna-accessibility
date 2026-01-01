"""Pre-outcome falsification fixtures and independent tiny reference calculations."""
if not __debug__:
    raise RuntimeError('Synthetic verification requires Python assertions enabled')
from fractions import Fraction as F
from collections import deque
import itertools
import json
from pathlib import Path
import math
from ordinal_geometry import features, order_certificate, widest_paths, bounded_pairs, graph_edges

HERE=Path(__file__).resolve().parent
checks=[]
q2=list(itertools.product(range(2),repeat=2))
x2=features(q2)
e2=graph_edges(q2,[{frozenset((0,1))}]*2)

def fixture(name,x,l,u,want):
    r=order_certificate(x,l,u)
    assert r['status']==want and r['exact_checked'], (name,r)
    checks.append({'name':name,'status':r['status'],'forced_orders':r['forced_orders']})
    return r

fixture('additive',x2,[0,2,1,3],[0,2,1,3],'feasible')
fixture('nonlinear_monotone',x2,[0,8,1,27],[0,8,1,27],'feasible')
fixture('constant_tied',x2,[1]*4,[1]*4,'feasible')
positive=fixture('reciprocal_sign',x2,[2,0,0,2],[2,0,0,2],'infeasible')
fixture('intervals_hide_point_failure',x2,[0]*4,[2]*4,'feasible')
fixture('unknown_states',x2,[0]*4,[None]*4,'feasible')
fixture('duplicate_features_incompatible',[[1,0],[1,0]],[0,2],[1,3],'infeasible')
fixture('partial_infinite_upper',x2,[2,0,0,2],[None,0,0,None],'infeasible')
fixture('touching_intervals_admit_constant',[[1,0],[1,0]],[0,1],[1,2],'feasible')
for perm in itertools.permutations(range(4)):
    vals=[2,0,0,2]
    assert order_certificate([x2[i] for i in perm],[vals[i] for i in perm],[vals[i] for i in perm])['status']=='infeasible'
checks.append({'name':'all_Q2_row_permutations','cases':24})

# All 81 three-valued Q2 fields versus independent bounded-integer latent enumeration.
for values in itertools.product(range(3),repeat=4):
    feasible=False
    for b,c in itertools.product(range(-4,5),repeat=2):
        latent=[b*s[0]+c*s[1] for s in q2]
        if all(latent[i]<latent[j] for i in range(4) for j in range(4) if values[i]<values[j]):
            feasible=True; break
    got=order_certificate(x2,values,values)
    assert (got['status']=='feasible')==feasible, values
    # A positive monotone transform changes magnitudes but not admissibility.
    transformed=[v*v*v+7 for v in values]
    assert order_certificate(x2,transformed,transformed)['status']==got['status']
checks.append({'name':'all_ternary_Q2_fields_and_monotone_transforms','cases':81,'reference':'integer latent enumeration [-4,4]^2; finite verification only'})

def threshold_reachable(scores,edges,i,j,t):
    if scores[i]<t or scores[j]<t:return False
    seen={i};todo=deque([i])
    while todo:
        v=todo.popleft()
        for a,b in edges:
            w=b if a==v else a if b==v else None
            if w is not None and w not in seen and scores[w]>=t:
                seen.add(w);todo.append(w)
    return j in seen

for values in itertools.product(range(3),repeat=4):
    result=widest_paths(values,e2)
    for i,j in itertools.combinations(range(4),2):
        ref=max(t for t in set(values) if threshold_reachable(values,e2,i,j,t))
        assert result[i][j]==ref
checks.append({'name':'widest_paths_vs_threshold_BFS','fields':81,'pairs_per_field':6})

states=list(itertools.product(range(3),range(3),range(2),range(2),range(2)))
xe=features(states)
edges=graph_edges(states,[{frozenset((0,1)),frozenset((1,2))}]*2+[{frozenset((0,1))}]*3)
assert len(states)==72 and len(edges)==204 and len(xe[0])==8
dip=[int(s[0]!=1) for s in states]
fixture('codon_gate_dip_full72',xe,dip,dip,'feasible')
bp=bounded_pairs(list(map(F,dip)),list(map(F,dip)),edges)
assert any(F(z['certified_interval_width'])>0 for z in bp)
checks.append({'name':'codon_gate_can_disconnect_without_interaction','vertices':72,'edges':204,'parameters':8})

# An unmeasured low bridge is allowed high by its upper bound, so prevents exclusion.
lower=list(map(F,[2,0,0,2]));upper=list(map(F,[2,0,0,2]));upper[1]=None
pairs=bounded_pairs(lower,upper,e2)
assert next(z for z in pairs if (z['a'],z['b'])==(0,3))['certified_interval_width']=='0'
checks.append({'name':'unknown_bridge_prevents_false_exclusion'})
iter_pairs=bounded_pairs([2,0,2],[2,3,2],iter([(0,1),(1,2)]))
assert next(z for z in iter_pairs if (z['a'],z['b'])==(0,2))['certified_nonempty'] is False
checks.append({'name':'edge_iterator_is_materialized_before_reuse'})

# Median normalized inputs preserve a ratio inside the separate-input ratio range.
for a in itertools.product((1,3,7),repeat=3):
    for b in itertools.product((2,5),repeat=3):
        totals=(13,19,23)
        af=[F(v,t) for v,t in zip(a,totals)]
        bf=[F(v,t) for v,t in zip(b,totals)]
        ratio=sorted(af)[1]/sorted(bf)[1]
        lo=min(v/w for v,w in zip(af,bf));hi=max(v/w for v,w in zip(af,bf))
        assert lo<=ratio<=hi
checks.append({'name':'median_input_ratio_within_scenario_envelope','cases':216})

report={'status':'PASS','scope':'Synthetic tests only; no biological outcome file opened',
        'checks':checks,'positive_fixture_certificate':positive,
        'proof_limit':'Finite enumeration validates implementations and examples, not general statements.'}
(HERE/'synthetic_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':report['status'],'checks':len(checks),'output':'synthetic_verification.json'}))
