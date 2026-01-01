"""Independent exact min/max Floyd oracle. No author implementation imports."""
if not __debug__:
    raise RuntimeError('Assertions required')
from fractions import Fraction as F
from itertools import product
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'joint_math'
PINS = {'claim_and_proof.md': 'd5ef26923060fd587bb35b0c2d03d6ae39c34fda1e9a967f87e723e00d5cadb1',
        'prior_map.md': '252e2ba00b776e5feabe24f91a074dbe09100bf99430a3d548b593d8df9e1b30'}
for name, pin in PINS.items():
    if hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() != pin:
        raise RuntimeError('Reviewed revision changed: '+name)

def violation(traits, scales, intervals):
    values = [F(0)]
    for x, a, (lower, upper) in zip(traits, scales, intervals):
        assert a > 0
        if lower is not None:
            values.append((lower-x)/a)
        if upper is not None:
            values.append((x-upper)/a)
    return max(values)

def connected(vertices, edges, s, t, allowed):
    if s not in allowed or t not in allowed:
        return False
    seen = {s}
    stack = [s]
    while stack:
        a = stack.pop()
        for u, v in edges:
            b = v if a == u else u if a == v else None
            if b is not None and b in allowed and b not in seen:
                seen.add(b)
                stack.append(b)
    return t in seen

def evaluate(example):
    ids = [v['id'] for v in example['vertices']]
    assert len(set(ids)) == len(ids)
    pos = {v:i for i,v in enumerate(ids)}
    scales = [F(x) for x in example['scales']]
    intervals = [(None if a is None else F(a), None if b is None else F(b))
                 for a,b in example['acceptable_intervals']]
    costs = {v['id']:violation([F(x) for x in v['traits']], scales, intervals)
             for v in example['vertices']}
    N = len(ids)
    dist = [[None]*N for _ in ids]
    for a in ids:
        dist[pos[a]][pos[a]] = costs[a]
    edges = [tuple(e) for e in example['edges']]
    for u,v in edges:
        assert u in pos and v in pos
        dist[pos[u]][pos[v]] = dist[pos[v]][pos[u]] = max(costs[u],costs[v])
    for k in range(N):
        for i in range(N):
            if dist[i][k] is None:
                continue
            for j in range(N):
                if dist[k][j] is None:
                    continue
                new = max(dist[i][k],dist[k][j])
                if dist[i][j] is None or new < dist[i][j]:
                    dist[i][j] = new
    answer = dist[pos[example['s']]][pos[example['t']]]
    critical = sorted(set(costs.values()))
    # A separate threshold traversal verifies attainability and the lower certificate.
    threshold_answer = next((c for c in critical if connected(ids,edges,example['s'],example['t'],
                             {a for a in ids if costs[a] <= c})),None)
    assert threshold_answer == answer
    labels = {v['id']:v['label'] for v in example['vertices']}
    hamming_edges = {tuple(sorted((a,b))) for a in ids for b in ids if a < b
                     and sum(x!=y for x,y in zip(labels[a],labels[b])) == 1}
    assert {tuple(sorted(e)) for e in edges} == hamming_edges
    computed = 'infinity' if answer is None else str(answer)
    assert computed == example['expected_R'], (example['name'],computed)
    return {'name':example['name'],'vertex_count':N,'edge_count':len(edges),
            'computed_R':computed,'expected_R':example['expected_R'],
            'minmax_Floyd_and_threshold_agree':True,'complete_induced_hamming_edges_checked':True,
            'critical_levels':[str(c) for c in critical],
            'endpoint_violations':[str(costs[example[k]]) for k in ('s','t')]}

raw_path = SOURCE/'oracle_examples.json'
raw = json.loads(raw_path.read_text())
results = [evaluate(x) for x in raw['examples']]

# Independent margin fixtures, all at the exact boundary or outside it.
extras = []
vertices = list(product((0,1),repeat=2))
for margin in [F(0),F(1,4),F(1),F(2)]:
    e = {'name':'margin '+str(margin),
         'vertices':[{'id':i,'label':list(x),'traits':[str(x[0]-x[1])]} for i,x in enumerate(vertices)],
         'edges':[[i,j] for i,x in enumerate(vertices) for j,y in enumerate(vertices)
                  if i<j and sum(a!=b for a,b in zip(x,y))==1],
         's':0,'t':3,'scales':['1'],'acceptable_intervals':[[-margin,margin]],
         'expected_R':str(max(F(0),1-margin))}
    extras.append(evaluate(e))

# Local path inside its endpoint box, but outside the endpoint affine line.
a,b,z=(F(0),F(0)),(F(1),F(1)),(F(0),F(1))
assert all(min(x,y)<=v<=max(x,y) for x,y,v in zip(a,b,z))
extras.append({'name':'endpoint box is larger than endpoint segment','eta_box':'0',
               'distance_infinity_to_segment':'1/2','intermediate_outside_endpoint_effect_span':True})

# Only controlling model error at endpoints is insufficient.
epsilon = F(1,4)
base = raw['examples'][0]
assert base['name'].startswith('Leu')
additive = {v['id']:-epsilon for v in base['vertices']}
true = {v['id']:(F(0) if v['id'] in (base['s'],base['t']) else -2*epsilon) for v in base['vertices']}
assert max(abs(true[a]-additive[a]) for a in true)==epsilon
assert {v['id']:F(v['traits'][0]) for v in raw['examples'][-1]['vertices']} == true
assert raw['examples'][-1]['edges'] == base['edges']
extras.append({'name':'sharp global approximation coefficient','epsilon':str(epsilon),
               'constant_model_rank':0,'constant_model_M':'0','constant_model_eta':'0',
               'true_R':results[-1]['computed_R'],'equals_2epsilon':results[-1]['computed_R']==str(2*epsilon)})

record = {'status':'passed','source_pins':PINS,
          'oracle_examples_sha256':hashlib.sha256(raw_path.read_bytes()).hexdigest(),
          'implementation':'Exact Fraction min/max Floyd, cross-checked by critical-level traversal; author verify_finite.py not imported or executed',
          'examples':results,'additional_boundary_checks':extras,
          'limitations':['Finite checks corroborate examples only; no finite verification is offered as proof of C5.1-C5.3.',
                         'Nucleotide synonymy application uses the explicitly declared standard-code domains; Hamming checks alone do not authenticate a biological clone.']}
(HERE/'checks.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':record['status'],'examples_checked':len(results),'extra_checks':len(extras)},indent=2))
