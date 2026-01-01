"""Separately authored exact LP and Boolean audit. No author solver imports."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
from time import perf_counter
from datetime import datetime, timezone
import hashlib
import json
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
START = perf_counter()
COUNTS = {'lp_calls': 0, 'largest_intermediate': 0, 'projection_steps': 0}


def guard():
    if perf_counter() - START > 180:
        raise RuntimeError('Declared 180 second cap reached')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduced(rows):
    """Canonicalize A.x <= b using positive scaling, never sign reversal."""
    out = {}
    for row, b in rows:
        scale = next((abs(v) for v in row if v), None)
        if scale is None:
            if b < 0:
                return None
            continue
        key = tuple(v / scale for v in row)
        value = b / scale
        if key not in out or value < out[key]:
            out[key] = value
    COUNTS['largest_intermediate'] = max(COUNTS['largest_intermediate'], len(out))
    if len(out) > 20000:
        raise RuntimeError('Declared 20000 inequality cap reached')
    return list(out.items())


def eliminate(rows, j):
    positive, negative, zero = [], [], []
    for row, b in rows:
        c = row[j]
        tail = row[:j] + row[j+1:]
        if c > 0:
            positive.append((tail, b, c))
        elif c < 0:
            negative.append((tail, b, c))
        else:
            zero.append((tail, b))
    for r, b, c in positive:
        for rr, bb, cc in negative:
            zero.append((tuple(x/c-y/cc for x, y in zip(r, rr)), b/c-bb/cc))
    COUNTS['projection_steps'] += 1
    return reduced(zero)


def projected_objective(rows, d):
    """Project theta out, retaining z as the final coordinate. Exact, no LP library."""
    guard()
    COUNTS['lp_calls'] += 1
    rows = reduced(rows)
    if rows is None:
        return None
    for remaining in range(d, 0, -1):
        # Elimination order controls intermediate size but cannot affect the projection.
        def size(j):
            pos = sum(r[j] > 0 for r, _ in rows)
            neg = sum(r[j] < 0 for r, _ in rows)
            return pos*neg, j
        j = min(range(remaining), key=size)
        rows = eliminate(rows, j)
        if rows is None:
            return None
    lower, upper = None, None
    for row, b in rows:
        c = row[0]
        if c > 0:
            upper = b/c if upper is None else min(upper, b/c)
        elif c < 0:
            lower = b/c if lower is None else max(lower, b/c)
        elif b < 0:
            return None
    if lower is not None and upper is not None and lower > upper:
        return None
    return lower, upper


def constraints(case):
    d = len(case['box'])
    out = []
    for j, (lo, hi) in enumerate(case['box']):
        unit = tuple(Q(int(i == j)) for i in range(d))
        out.extend([(unit, Q(hi)), (tuple(-x for x in unit), -Q(lo))])
    out.extend((tuple(map(Q, a)), Q(b)) for a, b in case.get('inequalities', []))
    return out


def connected(n, edges, s, t, mask):
    if not mask >> s & 1 or not mask >> t & 1:
        return False
    seen = {s}
    todo = [s]
    while todo:
        u = todo.pop()
        for a, b in edges:
            v = b if a == u else a if b == u else None
            if v is not None and mask >> v & 1 and v not in seen:
                seen.add(v)
                todo.append(v)
    return t in seen


def minimal_sets(n, edges, s, t):
    full = (1 << n)-1
    paths, cuts = [], []
    for mask in sorted(range(1 << n), key=lambda x: (x.bit_count(), x)):
        if connected(n, edges, s, t, mask) and not any(p & mask == p for p in paths):
            paths.append(mask)
        if not connected(n, edges, s, t, full ^ mask) and not any(k & mask == k for k in cuts):
            cuts.append(mask)
    convert = lambda values: [tuple(i for i in range(n) if mask >> i & 1) for mask in values]
    return convert(paths), convert(cuts)


def lp(base, affine, mode, d):
    # mode min: minimize z >= every affine expression; mode max is the dual hypograph.
    rows = [(a+(Q(0),), b) for a, b in base]
    for a0, *a in affine:
        if mode == 'min':
            rows.append((tuple(a)+(Q(-1),), -a0))
        else:
            rows.append((tuple(-v for v in a)+(Q(1),), a0))
    bounds = projected_objective(rows, d)
    if bounds is None:
        return None
    result = bounds[0 if mode == 'min' else 1]
    if result is None:
        raise AssertionError('Unexpected unbounded requested objective')
    return result


def check_case(case):
    n, s, t, edges = (case[k] for k in ('n', 's', 't', 'edges'))
    d = len(case['box'])
    base = constraints(case)
    coeffs = [tuple(map(Q, row)) for row in case['coefficients']]
    feasibility = lp(base, [tuple([Q(0)]*(d+1))], 'min', d)
    if feasibility is None:
        return {'id': case['id'], 'status': 'empty_parameter_set'}
    paths, cuts = minimal_sets(n, edges, s, t)
    if not paths:
        return {'id': case['id'], 'status': 'disconnected'}
    lower = [lp(base, [row], 'min', d) for row in coeffs]
    difference = lambda i, j: tuple(a-b for a, b in zip(coeffs[i], coeffs[j]))
    return {
        'id': case['id'], 'status': 'connected',
        'fixed': max(min(lower[v] for v in p) for p in paths),
        'guaranteed': min(lp(base, [coeffs[v] for v in k], 'min', d) for k in cuts),
        'possible': max(lp(base, [coeffs[v] for v in p], 'max', d) for p in paths),
        'barrier_minimum': min(lp(base, [difference(j, v) for v in p], 'min', d)
                               for p in paths for j in {s, t}),
        'barrier_maximum': max(lp(base, [difference(j, v) for v in k for j in {s, t}], 'max', d)
                               for k in cuts),
        'minimal_path_sets': paths, 'minimal_separators': cuts,
    }


def cap_bfs(n, edges, s, t, values):
    for h in sorted(set(values), reverse=True):
        mask = sum((1 << i) for i, x in enumerate(values) if x >= h)
        if connected(n, edges, s, t, mask):
            return h
    return None


def graph_suite():
    records = {'graphs': 0, 'connected_endpoint_graphs': 0, 'unique_minimal': 0,
               'multiple_minimal': 0, 'retained_sets_checked': 0, 'separation_witnesses': 0}
    for n in range(2, 6):
        pairs = list(combinations(range(n), 2))
        for graphmask in range(1 << len(pairs)):
            guard()
            edges = [edge for i, edge in enumerate(pairs) if graphmask >> i & 1]
            records['graphs'] += 1
            paths, _ = minimal_sets(n, edges, 0, n-1)
            if not paths:
                continue
            records['connected_endpoint_graphs'] += 1
            for mask in range(1 << n):
                contains = any(all(mask >> v & 1 for v in p) for p in paths)
                assert connected(n, edges, 0, n-1, mask) == contains
                records['retained_sets_checked'] += 1
            if len(paths) == 1:
                records['unique_minimal'] += 1
                continue
            records['multiple_minimal'] += 1
            p, q = map(set, paths[:2])
            inter = p & q
            assert not connected(n, edges, 0, n-1, sum(1 << v for v in inter))
            caps = []
            for theta in [Q(-1), Q(-1, 2), Q(0), Q(1, 2), Q(1)]:
                field = [Q(1) if i in inter else theta if i in p-q else -theta if i in q-p else Q(-1)
                         for i in range(n)]
                caps.append(cap_bfs(n, edges, 0, n-1, field))
            assert caps == [Q(1), Q(1, 2), Q(0), Q(1, 2), Q(1)]
            lower = [Q(1) if i in inter else Q(-1) for i in range(n)]
            assert cap_bfs(n, edges, 0, n-1, lower) == -1
            records['separation_witnesses'] += 1
    assert minimal_sets(1, [], 0, 0) == ([(0,)], [(0,)])
    records['coincident_endpoint_cases'] = 1
    return records


def literal_true(lit, assignment):
    i, positive = lit
    return bool(assignment[i]) == positive


def construction(formula, mode):
    endpoints = [0, 1]
    nextv = 2
    groups, labels, edges = [], {}, []
    for clause in formula:
        group = []
        for lit in clause:
            labels[nextv] = lit
            group.append(nextv)
            nextv += 1
        groups.append(group)
    if mode == 'possible':
        sequence = [[0]] + groups + [[1]]
        for left, right in zip(sequence, sequence[1:]):
            edges += list(product(left, right))
    else:
        for group in groups:
            path = [0]+group+[1]
            edges += list(zip(path, path[1:]))
    return nextv, edges, labels


def formula_suite():
    allclauses = [tuple(enumerate(signs)) for signs in product([False, True], repeat=3)]
    formulas = [[c for j, c in enumerate(allclauses) if mask >> j & 1] for mask in range(1, 256)]
    formulas += [
        [((0, True),), ((0, False),)],
        [((0, True), (0, True), (0, True)), ((0, False), (0, False), (0, False))],
        [((0, True), (0, False))],
        [((0, True), (1, True)), ((0, False), (1, True)), ((1, False),)],
    ]
    records = {'formula_cases': 0, 'unsatisfiable': 0, 'boolean_scenarios': 0, 'ternary_scenarios': 0}
    for formula in formulas:
        guard()
        d = 1+max(i for c in formula for i, _ in c)
        assignments = list(product([0, 1], repeat=d))
        sat = [a for a in assignments if all(any(literal_true(lit, a) for lit in c) for c in formula)]
        possible_graph = construction(formula, 'possible')
        guaranteed_graph = construction(formula, 'guaranteed')
        def connects(graph, theta, mode):
            n, edges, labels = graph
            values = [Q(1)]*n
            for v, (i, positive) in labels.items():
                val = theta[i] if positive else 1-theta[i]
                values[v] = val if mode == 'possible' else 1-val
            threshold = Q(1) if mode == 'possible' else Q(1, 2)
            mask = sum(1 << i for i, value in enumerate(values) if value >= threshold)
            return connected(n, edges, 0, 1, mask)
        observed = []
        for a in assignments:
            p = connects(possible_graph, a, 'possible')
            g = connects(guaranteed_graph, a, 'guaranteed')
            truth = all(any(literal_true(lit, a) for lit in c) for c in formula)
            assert p == truth and g == (not truth)
            observed.append(p)
            records['boolean_scenarios'] += 1
        assert any(observed) == bool(sat)
        # Exact tie and fractional probes. Universal implication itself is proved in the review.
        for theta in product([Q(0), Q(1, 2), Q(1)], repeat=d):
            g = connects(guaranteed_graph, theta, 'guaranteed')
            if not sat:
                assert g
            rounded = tuple(v >= Q(1, 2) for v in theta)
            false_clause = any(not any(literal_true(lit, rounded) for lit in c) for c in formula)
            if false_clause:
                assert g
            records['ternary_scenarios'] += 1
        records['formula_cases'] += 1
        records['unsatisfiable'] += not bool(sat)
    return records


def main():
    started = datetime.now(timezone.utc).isoformat()
    inputs = [HERE/'check_exact.py', HERE/'contract.json', HERE/'proof_reconstruction.md']
    inputs += [HERE.parent/name for name in ('claims.md', 'affine_polytope.py', 'protocol.json', 'cases.json')]
    before = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    freeze = {'source_hashes': before, 'created_before_run': started,
              'author_results_read': False, 'command': [sys.executable, str(Path(__file__).resolve())]}
    (HERE/'freeze.json').write_text(json.dumps(freeze, indent=2)+'\n', encoding='utf-8')
    out = HERE/'independent_results.json'
    if out.exists():
        raise RuntimeError('Refusing to replace existing independent receipt')
    results = {'provenance': 'Separately authored with inherited context, not blinded', 'status': 'running', 'cases': []}
    failure = None
    try:
        # Known LP controls and infeasible constant/strict bound cases.
        assert projected_objective([((Q(1),), Q(2)), ((Q(-1),), Q(3))], 0) == (Q(-3), Q(2))
        assert projected_objective([((Q(1),), Q(0)), ((Q(-1),), Q(-1))], 0) is None
        cases = json.loads((HERE.parent/'cases.json').read_text(encoding='utf-8'))
        assert len(cases) == 81
        for case in cases:
            results['cases'].append(check_case(case))
        results['graphs'] = graph_suite()
        results['formulas'] = formula_suite()
        results['status'] = 'completed'
    except Exception as exc:
        failure = f'{type(exc).__name__}: {exc}'
        results['status'] = 'failed'
        results['failure'] = failure
    elapsed = perf_counter()-START
    after = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    if before != after:
        results['status'] = 'inputs-changed'
        failure = 'Input changed during audit'
    results['counters'] = COUNTS
    results['runtime_seconds'] = elapsed
    results['author_results_read'] = False
    out.write_text(json.dumps(results, indent=2, default=str)+'\n', encoding='utf-8')
    gitroot = ROOT/'outputs/synonymous-rna-accessibility'
    commit = subprocess.check_output(['git', '-C', str(gitroot), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = bool(subprocess.check_output(['git', '-C', str(gitroot), 'status', '--porcelain'], text=True).strip())
    manifest = {
        'schema_version': 2, 'claim_id': 'U1_U3_SEPARATELY_AUTHORED_REVIEW_V1',
        'repository': {'commit': commit, 'dirty': dirty},
        'command': [sys.executable, str(Path(__file__).resolve())],
        'environment': {'software': [{'name': 'Python', 'version': platform.python_version()}],
                        'hardware': platform.platform()+'; '+platform.machine()},
        'mathematics': {'assertion_tested': 'Exact independent LP extrema for all 81 supplied cases and Boolean graph/reduction checks',
                        'coefficient_domain': 'Rational, fractions.Fraction',
                        'conventions': 'Closed vertex threshold; endpoints included; minimal vertex path sets and separators',
                        'inputs': list(before), 'bounds': {'affine_cases':81, 'graph_vertices_max':5, 'latent_dimension_max':3},
                        'non_claims': ['Universal proof from finite checks', 'Blinded review', 'Biological calibration', 'Generic machinery novelty']},
        'randomness': {'used':False, 'generator':'', 'seed':None},
        'run': {'started_at':started, 'runtime_seconds':elapsed, 'exit_status':int(failure is not None), 'status':results['status'],
                'resource_limits': {'memory_bytes':None, 'cpu_seconds':180, 'max_cores':1, 'max_threads':1}},
        'outputs': [{'path':str(out.relative_to(ROOT)), 'sha256':sha(out), 'kind':'result'}],
        'checks': ['Exact Fourier-Motzkin elimination', 'No author solver imports', 'No author results read before receipt', COUNTS],
        'result': 'Independent exact calculations saved; comparison with author receipt is a subsequent action' if not failure else failure,
        'residual_risks': ['Same model family and inherited context; not a blinded review', 'Finite cases are not universal proof'],
        'input_artifacts': [{'path':p, 'sha256':h, 'sha256_after':after[p]} for p,h in before.items()],
        'declared_results': [str(out.relative_to(ROOT))],
    }
    (HERE/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status':results['status'], 'cases':len(results['cases']), 'seconds':elapsed, 'counters':COUNTS, 'failure':failure}))
    return int(failure is not None)


if __name__ == '__main__':
    raise SystemExit(main())
