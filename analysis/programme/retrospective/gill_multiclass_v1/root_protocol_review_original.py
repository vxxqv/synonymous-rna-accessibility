"""Independent root review of sequence-only projection and membership.

No workbook or measured outcome is opened. The full indicator design and
Moore-Penrose projection provide a separate numerical check of the proposed
factorial-mean residual operator on every eligible class architecture.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
inventory = json.loads((HERE/'inventory.json').read_text(encoding='utf-8'))
proposal = json.loads((HERE/'freeze_proposal.json').read_text(encoding='utf-8'))
rows = [x for x in inventory['classes'] if x['candidate_reproducibility_eligible'] and x['candidate_partition'] in ('development','validation')]
assert len(rows) == 2176
assert sum(x['observed_vertices'] for x in rows) == 38338
assert {x['class_id'] for x in rows} == set(proposal['class_ids_by_partition']['development']) | set(proposal['class_ids_by_partition']['validation'])
assert set(proposal['class_ids_by_partition']['development']).isdisjoint(proposal['class_ids_by_partition']['validation'])
seen = set()
checked = set()
max_error = 0.
for row in rows:
    keys = row['observed_keys']
    assert not seen.intersection(keys)
    seen.update(keys)
    factors = [(k[:1],k[1:4],k[4:7],k[7:]) for k in keys]
    levels = [sorted({f[j] for f in factors}) for j in range(4)]
    sizes = tuple(len(v) for v in levels)
    assert list(sizes) == row['codon_factor_sizes']
    assert np.prod(sizes) == len(keys)
    assert 1+sum(s-1 for s in sizes) == row['codon_main_effect_rank']
    if sizes in checked:
        continue
    checked.add(sizes)
    n = len(keys)
    design = np.column_stack([np.ones(n)]+[np.array([f[j] == level for f in factors],float) for j in range(4) for level in levels[j][1:]])
    reference = design @ np.linalg.pinv(design)
    proposed = sum(np.array([[float(a[j] == b[j])/(n/sizes[j]) for b in factors] for a in factors]) for j in range(4))-3*np.ones((n,n))/n
    error = float(np.max(np.abs(reference-proposed)))
    max_error = max(max_error,error)
    assert error < 1e-12
    assert np.max(np.abs(proposed@proposed-proposed)) < 1e-12
    assert np.max(np.abs(proposed-proposed.T)) < 1e-12
    assert np.linalg.matrix_rank(design) == row['codon_main_effect_rank']
    # Normalization decomposition is algebraic, checked on deterministic synthetic fields.
    dummy = np.arange(n*6,dtype=float).reshape(n,6)**2 % 31
    residual = (np.eye(n)-proposed)@dummy
    u=(residual[:,0]+residual[:,1]-residual[:,3]-residual[:,4])/2
    v=residual[:,2]-residual[:,5]
    r=(residual[:,0]+residual[:,1])/2
    d=(residual[:,3]+residual[:,4])/2
    terms = [r@residual[:,2],d@residual[:,5],-(r@residual[:,5]),-(d@residual[:,2])]
    assert np.isclose(u@v,sum(terms),atol=1e-9,rtol=1e-12)
result = {'status':'passed sequence-only protocol checks','classes':len(rows),'distinct_keys':len(seen),
          'factor_architectures':len(checked),'maximum_projection_error':max_error,
          'workbook_opened':False,'outcome_values_opened':False,
          'scope':'separately authored root check; numerical matrix identities support algebra, no biological independence inferred',
          'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'inventory.json',HERE/'freeze_proposal.json',Path(__file__)]}}
(HERE/'root_protocol_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))
