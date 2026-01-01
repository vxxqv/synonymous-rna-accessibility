"""Subsequent receipt comparison and explicitly labelled additional challenges."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
from check_exact import check_case, cap_bfs

HERE = Path(__file__).resolve().parent
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
own_path = HERE/'independent_results.json'
author_path = HERE.parent/'recorded_check/results.json'
own = json.loads(own_path.read_text(encoding='utf-8'))
author = json.loads(author_path.read_text(encoding='utf-8'))
lookup = {c['id']: c['exact'] for c in author['cases']}
assert own['status'] == author['status'] == 'completed'
assert len(lookup) == len(own['cases']) == 81
keys = ['fixed', 'guaranteed', 'possible', 'barrier_minimum', 'barrier_maximum']
comparisons = 0
for case in own['cases']:
    reference = lookup[case['id']]
    assert case['status'] == reference['status'], case['id']
    if case['status'] == 'connected':
        for key in keys:
            assert Q(case[key]) == Q(reference[key]), (case['id'], key, case[key], reference[key])
            comparisons += 1

suggested = {
    'id': 'root_suggested_genuine_two_dimensional_corner_failure',
    'n': 6, 's': 0, 't': 5,
    'edges': [[0, j] for j in range(1, 5)]+[[j, 5] for j in range(1, 5)],
    'box': [[-1, 1], [-1, 1]],
    'coefficients': [[1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, 0, 0]],
}
extra = check_case(suggested)
assert [extra[k] for k in keys] == [Q(-1), Q(0), Q(1), Q(0), Q(1)]
corners = []
for x in [-1, 1]:
    for y in [-1, 1]:
        scores = [Q(a+b*x+c*y) for a, b, c in suggested['coefficients']]
        c = cap_bfs(6, suggested['edges'], 0, 5, scores)
        corners.append({'theta': [x, y], 'capacity': c, 'barrier': 1-c})
        assert c == 1
origin = cap_bfs(6, suggested['edges'], 0, 5, [Q(a) for a, _, _ in suggested['coefficients']])
assert origin == 0

# This original fixture's name is not evidence that a corner shortcut failed.
original = next(c for c in json.loads((HERE.parent/'cases.json').read_text())
                if c['id'] == 'two_dimensional_corner_failure')
original_corners = []
for x in [-1, 1]:
    for y in [-1, 1]:
        scores = [Q(a+b*x+c*y) for a, b, c in original['coefficients']]
        original_corners.append(cap_bfs(4, original['edges'], 0, 3, scores))
assert min(original_corners) == 0 and max(original_corners) == 2

result = {
    'status': 'completed', 'author_receipt_sha256': digest(author_path),
    'independent_receipt_sha256': digest(own_path),
    'case_status_comparisons': 81, 'exact_value_comparisons': comparisons,
    'mismatches': 0,
    'original_fixture_naming_defect': {
        'id': original['id'], 'capacity_at_corners': original_corners,
        'conclusion': 'Does not demonstrate failure of parameter corner extrema; all capacity and barrier extrema already occur at corners.'},
    'additional_case_provenance': 'Suggested by root during review, after independent first receipt. Not a blinded reviewer challenge and not part of the frozen 81 author cases.',
    'additional_case': suggested, 'exact_additional_result': extra,
    'corner_results': corners, 'origin_capacity': origin,
}
output = HERE/'comparison.json'
if output.exists():
    raise RuntimeError('Refusing to replace comparison receipt')
output.write_text(json.dumps(result, indent=2, default=str)+'\n', encoding='utf-8')
print(json.dumps({'cases':81, 'exact_value_comparisons':comparisons, 'mismatches':0,
                  'genuine_2d_corner_failure_verified':True}))
