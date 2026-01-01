"""Independent exact short-support enumeration. No author analysis imports."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys

SOURCE_SHA = '3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'
TSTAR = Q(9128860950384, 8609574801263)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def weights(rows):
    """Positive null vector for a minimal support of size <=3, else None.

    Rank-one three-row positive dependencies cannot be minimal: zero rows or
    a pair of opposite rays already give a proper positive dependency.
    """
    n, d = len(rows), len(rows[0])
    if n == 1:
        w = [1] if not any(rows[0]) else None
    elif n == 2:
        k = next((k for k in range(d) if rows[0][k] or rows[1][k]), None)
        if k is None:
            return None
        w = [-rows[1][k], rows[0][k]]
    elif n == 3:
        w = None
        for p, q in combinations(range(d), 2):
            a, b = [r[p] for r in rows], [r[q] for r in rows]
            trial = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
            if any(trial):
                w = trial
                break
    else:
        raise ValueError('Only short supports are defined')
    if w is None or not (all(v > 0 for v in w) or all(v < 0 for v in w)):
        return None
    w = [Q(v, sum(w)) for v in w]
    if any(sum(w[k]*rows[k][j] for k in range(n)) for j in range(d)):
        return None
    return w

def enumerate_minimal(rows):
    found, counts = [], {}
    for size in (1, 2, 3):
        total = skipped = 0
        entries = []
        for support in combinations(range(len(rows)), size):
            total += 1
            if any(set(old['support']).issubset(support) for old in found):
                skipped += 1
                continue
            w = weights([rows[k] for k in support])
            if w is not None:
                entries.append({'support': list(support), 'weights': w})
        found.extend(entries)
        counts[str(size)] = {'combinations': total, 'nonminimal_skipped': skipped, 'certificates': len(entries)}
    return found, counts

def synthetic():
    checks = 0
    cases = [([[0, 0]], [Q(1)]), ([[1, 0]], None),
             ([[2, 0], [-3, 0]], [Q(3, 5), Q(2, 5)]),
             ([[1, 0], [2, 0]], None),
             ([[1, 0], [0, 1], [-1, -1]], [Q(1, 3)]*3),
             ([[1, 0], [0, 1], [1, 1]], None),
             ([[1, 0, 0], [0, 1, 0], [-1, -1, 1]], None),
             ([[1, 0], [-1, 0], [2, 0]], None)]
    for rows, expected in cases:
        require(weights(rows) == expected, 'Synthetic null-vector failure')
        checks += 1
    found, counts = enumerate_minimal([[0, 0], [1, 0], [-1, 0]])
    require([f['support'] for f in found] == [[0], [1, 2]], 'Minimality fixture failed')
    checks += 1
    # Exhaustive rational products and a constructive split of each multiplier.
    for gamma, eta in product([Q(1), Q(3, 2), Q(2)], repeat=2):
        for r in [1/(gamma*eta), Q(1), gamma*eta, Q(1, 2)*(1/(gamma*eta)+gamma*eta)]:
            lo, hi = max(1/gamma, r/eta), min(gamma, r*eta)
            require(lo <= hi and 1/eta <= r/lo <= eta, 'Product factorization failed')
            checks += 1
    return checks

def run(source, author_result, robustness):
    require(digest(source) == SOURCE_SHA, 'Frozen source hash mismatch')
    j = json.loads(Path(source).read_bytes())
    vertices = j['vertices_data']
    require(len(vertices) == 72 and [v['vertex'] for v in vertices] == list(range(72)), 'Vertex identity mismatch')
    states = [tuple(v['codon_states']) for v in vertices]
    starts = (23, 38, 29, 32, 41)
    require(all(tuple(v['sequence'][k:k+3] for k in starts) == s for v, s in zip(vertices, states)), 'Sequence/codon-state mismatch')
    levels = [sorted({s[k] for s in states}) for k in range(5)]
    columns = [(k, level) for k, group in enumerate(levels) for level in group[1:]]
    x = [[1] + [int(s[k] == level) for k, level in columns] for s in states]
    require(all(len(row) == 8 for row in x), 'Design width mismatch')
    lower, upper = ([Q(v[field]) for v in vertices] for field in ('lower', 'upper'))
    require(all(0 < lo <= hi for lo, hi in zip(lower, upper)), 'Positive interval assumption fails')
    require(lower[46] == upper[46] == 1, 'Anchor is not the fixed WT identity')
    edges = []
    for i, k in product(range(72), repeat=2):
        if upper[i] < lower[k]:
            exponent = int(i != 46) + int(k != 46)
            require(exponent in (1, 2), 'Invalid widening exponent')
            threshold = (lower[k]/upper[i]) ** (2//exponent)
            edges.append({'i': i, 'j': k, 'row': [b-a for a, b in zip(x[i], x[k])], 'loss_t': threshold})
    require(len(edges) == 74, 'Forced-order count differs')
    certs, counts = enumerate_minimal([e['row'] for e in edges])
    for cert in certs:
        cert['edges'] = [edges[k] for k in cert['support']]
        cert['cancellation'] = [sum(w*e['row'][d] for w, e in zip(cert['weights'], cert['edges'])) for d in range(8)]
        cert['vertices'] = sorted({v for e in cert['edges'] for v in (e['i'], e['j'])})
        cert['persistence_t_exclusive'] = min(e['loss_t'] for e in cert['edges'])
    author = json.loads(Path(author_result).read_bytes())
    encoded_edges = json.loads(json.dumps(edges, default=str))
    require(encoded_edges == author['all_forced_orders'], 'Author forced-order vectors/thresholds mismatch')
    require(counts == author['enumeration'], 'Author enumeration counts mismatch')
    require(len(certs) == len(author['certificates']), 'Author certificate count mismatch')
    for c, a in zip(certs, author['certificates']):
        require(c['support'] == a['edge_indices'], 'Support mismatch')
        require(list(map(str, c['weights'])) == a['weights'], 'Weight mismatch')
        require(c['vertices'] == a['vertices'], 'Vertex support mismatch')
        require(str(c['persistence_t_exclusive']) == a['persistence_t_exclusive'], 'Persistence mismatch')
    require(len(certs) == 1 and len(certs[0]['support']) == 3, 'Unique global minimum claim fails')
    require(certs[0]['persistence_t_exclusive'] == TSTAR, 'Existing sharp threshold mismatch')
    # Check the supplied boundary primal directly; no LP package or author code.
    previous = json.loads(Path(robustness).read_bytes())['full']['multiplicative_squared']
    require(Q(previous['parameter_exact']) == TSTAR, 'Previous threshold mismatch')
    beta = list(map(Q, previous['feasible_at_boundary']['beta']))
    require(len(beta) == 8, 'Boundary beta dimension mismatch')
    latent = [sum(a*b for a, b in zip(row, beta)) for row in x]
    boundary_orders = []
    completion_checks = 0
    for i, k in product(range(72), repeat=2):
        power = int(i != 46) + int(k != 46)
        if upper[i]**2 * TSTAR**power < lower[k]**2:
            require(latent[i] < latent[k], 'Boundary primal violates a forced order')
            boundary_orders.append((i, k))
        if latent[i] <= latent[k]:
            require(lower[i]**2 <= upper[k]**2 * TSTAR**power, 'No isotonic interval completion at boundary')
            completion_checks += 1
    # All original certificate comparisons hold below t*, one is tied at t*.
    boundary_ties = []
    for e in certs[0]['edges']:
        i, k = e['i'], e['j']
        p = int(i != 46) + int(k != 46)
        require(upper[i]**2 < lower[k]**2, 'Certificate comparison not initially strict')
        require(upper[i]**2 * ((1+TSTAR)/2)**p < lower[k]**2, 'Certificate lost below threshold')
        if upper[i]**2 * TSTAR**p == lower[k]**2:
            boundary_ties.append((i, k))
    require(boundary_ties == [(42, 13)], 'Unexpected threshold tie')
    return {'status': 'passed', 'source_sha256': digest(source), 'author_result_sha256': digest(author_result),
            'previous_robustness_sha256': digest(robustness), 'states': 72, 'design_columns': 8,
            'levels': levels, 'columns': [['intercept']] + [list(c) for c in columns],
            'forced_orders': len(edges), 'enumeration': counts, 'certificates': certs,
            'minimum_comparisons': 3, 'boundary_tied_comparisons': boundary_ties,
            'boundary_forced_orders': len(boundary_orders), 'boundary_beta': beta,
            'isotonic_completion_pair_checks': completion_checks, 'synthetic_checks': synthetic(),
            'author_analysis_imported_or_executed': False,
            'limitations': ['Finite support audit covers every support of sizes 1,2,3 only.',
                            'Global feasibility at the boundary is independently checked from supplied exact beta.',
                            'No new source outcomes, calibrated uncertainty, molecular mechanism or independent biological validation.']}

def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is refused')
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', required=True)
    ap.add_argument('--author-result', required=True)
    ap.add_argument('--robustness', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    require(not out.exists(), 'Refusing output overwrite')
    result = run(args.source, args.author_result, args.robustness)
    result['reviewer_code_sha256'] = digest(__file__)
    result['python'] = sys.version
    out.write_text(json.dumps(result, indent=2, default=str)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'enumeration': result['enumeration'], 'sha256': digest(out)}))

if __name__ == '__main__':
    main()
