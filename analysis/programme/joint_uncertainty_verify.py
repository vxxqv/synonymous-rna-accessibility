"""Independent simple-path/piecewise-envelope audit for affine joint uncertainty.

The reference enumerates simple paths and merges their line envelopes. It does
not call the production widest-path routine or its global score arrangement.
All arithmetic and comparisons use fractions; no datasets are opened.
"""
import argparse
from collections import Counter
import datetime
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

import joint_uncertainty_affine as implementation

ROOT = Path(__file__).resolve().parents[2]
COUNTS = Counter()


def all_paths(n, edges, source, target):
    neighbors = [[] for _ in range(n)]
    for i, j in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)
    paths = []
    def search(path):
        if path[-1] == target:
            paths.append(tuple(path))
        else:
            for v in neighbors[path[-1]]:
                if v not in path:
                    search(path+[v])
    search([source])
    return paths


def line_on(pieces, theta):
    return next((a, b) for left, right, a, b in pieces if left <= theta <= right)


def merge(first, second, choose_max):
    cuts = {x for piece in first+second for x in piece[:2]}
    for l, r, a, b in first:
        for p, q, c, d in second:
            left, right = max(l, p), min(r, q)
            if left < right and b != d:
                crossing = (c-a)/(b-d)
                if left < crossing < right:
                    cuts.add(crossing)
    result = []
    cuts = sorted(cuts)
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left+right)/2
        a, b = line_on(first, midpoint)
        c, d = line_on(second, midpoint)
        take_first = (a+b*midpoint >= c+d*midpoint) if choose_max else (a+b*midpoint <= c+d*midpoint)
        aa, bb = (a, b) if take_first else (c, d)
        if result and result[-1][1] == left and result[-1][2:] == (aa, bb):
            old = result.pop()
            result.append((old[0], right, aa, bb))
        else:
            result.append((left, right, aa, bb))
    return result


def independent_reference(n, edges, coefficients, source, target, lower, upper):
    paths = all_paths(n, edges, source, target)
    COUNTS['simple_paths_enumerated'] += len(paths)
    if not paths:
        return None, paths
    fixed = max(min(min(coefficients[v][0]+coefficients[v][1]*theta for v in path)
                    for theta in (lower, upper)) for path in paths)
    if lower == upper:
        value = max(min(coefficients[v][0]+coefficients[v][1]*lower for v in path) for path in paths)
        endpoint = min(coefficients[v][0]+coefficients[v][1]*lower for v in (source, target))
        return {'fixed': fixed, 'capacity_range': (value, value), 'barrier_range': (endpoint-value, endpoint-value)}, paths
    capacity = None
    for path in paths:
        path_function = [(lower, upper, *coefficients[path[0]])]
        for v in path[1:]:
            path_function = merge(path_function, [(lower, upper, *coefficients[v])], False)
        capacity = path_function if capacity is None else merge(capacity, path_function, True)
    endpoints = merge([(lower, upper, *coefficients[source])], [(lower, upper, *coefficients[target])], False)
    cuts = sorted({x for piece in capacity+endpoints for x in piece[:2]})
    capacities, barriers = [], []
    for theta in cuts:
        a, b = line_on(capacity, theta)
        c, d = line_on(endpoints, theta)
        capacities.append(a+b*theta)
        barriers.append(c+d*theta-a-b*theta)
    return {'fixed': fixed, 'capacity_range': (min(capacities), max(capacities)),
            'barrier_range': (min(barriers), max(barriers)), 'capacity_pieces': capacity}, paths


def check(n, edges, coefficients, source=0, target=None, lower=F(-1), upper=F(1), corner_check=True):
    target = n-1 if target is None else target
    coefficients = tuple((F(a), F(b)) for a, b in coefficients)
    result = implementation.analyze(n, edges, coefficients, source, target, lower, upper)
    reference, paths = independent_reference(n, edges, coefficients, source, target, lower, upper)
    if reference is None:
        assert result['status'] == 'disconnected' and result['capacity_range'] is None and result['barrier_range'] is None
        COUNTS['disconnected_instances'] += 1
    else:
        assert result['status'] == 'connected'
        assert result['fixed'] == reference['fixed']
        assert result['capacity_range'] == reference['capacity_range']
        assert result['barrier_range'] == reference['barrier_range']
        assert result['fixed'] <= result['guaranteed'] <= result['possible']
        assert result['box_guaranteed'] <= result['guaranteed'] and result['possible'] <= result['box_possible']
        assert result['box_barrier_range'][0] <= result['barrier_range'][0] <= result['barrier_range'][1] <= result['box_barrier_range'][1]
        fixed_path = result['fixed_path']
        assert fixed_path in paths
        assert min(min(coefficients[v][0]+coefficients[v][1]*theta for v in fixed_path) for theta in (lower, upper)) == result['fixed']
        assert len(result['breakpoints']) <= n*(n-1)//2+2
        for cell in result['cells']:
            assert cell['path'] in paths
            for k in (1, 2, 3):
                theta = cell['lower']+(cell['upper']-cell['lower'])*F(k, 4)
                scores = [a+b*theta for a, b in coefficients]
                exact = max(min(scores[v] for v in path) for path in paths)
                a, b = cell['capacity_line']
                c, d = cell['barrier_line']
                assert a+b*theta == exact == min(scores[v] for v in cell['path'])
                assert c+d*theta == min(scores[source], scores[target])-exact
                COUNTS['interior_cell_certificate_checks'] += 1
        if corner_check:
            lower_scores = [min(a+b*lower, a+b*upper) for a, b in coefficients]
            upper_scores = [max(a+b*lower, a+b*upper) for a, b in coefficients]
            cap, barrier = [], []
            for choices in itertools.product((0, 1), repeat=n):
                scores = [upper_scores[v] if choices[v] else lower_scores[v] for v in range(n)]
                value = max(min(scores[v] for v in path) for path in paths)
                cap.append(value)
                barrier.append(min(scores[source], scores[target])-value)
            assert (min(cap), max(cap)) == (result['box_guaranteed'], result['box_possible'])
            assert (min(barrier), max(barrier)) == result['box_barrier_range']
            COUNTS['box_corners_enumerated'] += 2**n
        COUNTS['connected_instances'] += 1
        COUNTS['strict_adaptivity_gaps'] += result['fixed'] < result['guaranteed']
    COUNTS['total_instances'] += 1
    return result


def serialized(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serialized(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialized(v) for v in value]
    return value


def main():
    if not __debug__:
        raise RuntimeError('exact verification requires assertions')
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    out = (ROOT/args.outdir).resolve()
    if ROOT not in out.parents or out.exists():
        raise ValueError('fresh repository output directory required')
    out.mkdir(parents=True)
    files = ['analysis/programme/joint_uncertainty_affine.py', 'analysis/programme/joint_uncertainty_verify.py']
    before = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}
    started, clock = datetime.datetime.now(datetime.timezone.utc).isoformat(), time.perf_counter()
    diamond_edges = ((0, 1), (1, 3), (0, 2), (2, 3))
    diamond = check(4, diamond_edges, ((1, 0), (0, 1), (0, -1), (1, 0)))
    assert (diamond['fixed'], diamond['guaranteed'], diamond['possible']) == (F(-1), F(0), F(1))
    assert diamond['barrier_range'] == (F(0), F(1)) and diamond['box_barrier_range'] == (F(0), F(2))
    assert diamond['barrier_maximizer'] == 0 and diamond['capacity_minimizer'] == 0
    assert all(row['capacity'] >= F(1, 2) for row in diamond['evaluations'] if row['theta'] in (-1, 1))
    assert diamond['guaranteed'] < F(1, 2)
    shifted = check(4, diamond_edges, ((4, 2), (3, 3), (3, 1), (4, 2)))
    assert shifted['barrier_range'] == diamond['barrier_range']
    scaled = check(4, diamond_edges, ((5, 0), (2, 3), (2, -3), (5, 0)))
    assert scaled['barrier_range'] == (F(0), F(3))
    reference_error = check(3, ((0, 1), (1, 2)), ((1, 1), (0, 0), (1, 0)))
    assert reference_error['barrier_range'] == (F(0), F(1))
    normalized_reference = check(3, ((0, 1), (1, 2)), ((0, 0), (-1, -1), (0, -1)))
    assert normalized_reference['barrier_range'] == reference_error['barrier_range']
    # Complete three-vertex universe: all graph topologies and all nine affine
    # coefficient choices per vertex from {-1,0,1} x {-1,0,1}.
    possible_edges = list(itertools.combinations(range(3), 2))
    coefficient_choices = tuple(itertools.product((-1, 0, 1), repeat=2))
    for mask in range(8):
        edges = tuple(e for i, e in enumerate(possible_edges) if mask & (1 << i))
        for coefficients in itertools.product(coefficient_choices, repeat=3):
            value = check(3, edges, coefficients)
            if value['status'] == 'connected':
                assert value['fixed'] == value['guaranteed']  # no 3-vertex adaptivity gap
    print('all5832 three-vertex affine graph instances passed', flush=True)
    possible_edges = list(itertools.combinations(range(4), 2))
    for mask in range(64):
        edges = tuple(e for i, e in enumerate(possible_edges) if mask & (1 << i))
        for slopes in itertools.product((-1, 0, 1), repeat=4):
            coefficients = tuple(zip((1, 0, 0, 1), slopes))
            check(4, edges, coefficients)
    print('all5184 four-vertex graph/slope instances passed', flush=True)
    seed = 2026102502
    rng = random.Random(seed)
    for trial in range(300):
        n = 2+trial % 6
        edges = tuple(e for e in itertools.combinations(range(n), 2) if rng.random() < .4)
        coefficients = tuple((F(rng.randint(-5, 5), rng.randint(1, 5)), F(rng.randint(-5, 5), rng.randint(1, 5))) for _ in range(n))
        lower = F(rng.randint(-5, 0), rng.randint(1, 5))
        upper = F(rng.randint(0, 5), rng.randint(1, 5))
        check(n, edges, coefficients, lower=lower, upper=upper)
    # Identical functions, equal slopes, singleton interval, same endpoint,
    # negative values, disconnected graph and duplicate edges.
    boundaries = [
        check(1, (), ((-2, 3),), 0, 0),
        check(3, ((0, 1), (1, 2)), ((1, 2), (1, 2), (1, 2))),
        check(3, ((0, 1), (1, 2)), ((1, 2), (-1, 2), (2, 2))),
        check(3, ((0, 1), (1, 2)), ((0, 1), (1, -1), (2, 3)), lower=F(2, 3), upper=F(2, 3)),
        check(3, ((0, 1),), ((1, 1), (0, -1), (1, 0))),
        check(3, ((0, 1), (1, 2)), ((0, 1), (1, -1), (2, 3)), 1, 1),
        implementation.analyze(2, ((0, 1), (1, 0)), ((1, 0), (0, 0)), 0, 1, -1, 1)]
    for parameters in [(0, (), (), 0, 0, 0, 1), (1, (), ((0, 0),), 0, 0, 2, 1),
                       (1, ((0, 0),), ((0, 0),), 0, 0, 0, 1), (1, (), ((0.5, 0),), 0, 0, 0, 1)]:
        try:
            implementation.analyze(*parameters)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid input accepted')
    COUNTS['invalid_input_cases'] = 4
    payload = {'claim_id': 'JOINT_AFFINE_UNCERTAINTY_EXACT_EXTREMA', 'counts': dict(COUNTS),
               'diamond': diamond, 'shared_affine_offset': shifted, 'positive_scale': scaled,
               'reference_coordinate_error': reference_error, 'after_shared_reference_normalization': normalized_reference,
               'boundary_cases': boundaries, 'seed': seed,
               'nonclaims': ['Deterministic uncertainty bounds do not provide probabilistic calibration.',
                             'Generic minimax inequality, robust bottleneck envelopes and adjustable robustness are prior principles.',
                             'No biological measurements or interventions are used.']}
    output = out/'joint_uncertainty_results.json'
    output.write_text(json.dumps(serialized(payload), indent=2)+'\n')
    after = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}
    if before != after:
        raise RuntimeError('source changed during verification')
    manifest = {'schema_version': 2, 'claim_id': payload['claim_id'],
                'repository': {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(), 'dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))},
                'command': [sys.executable, 'analysis/programme/joint_uncertainty_verify.py', '--outdir', str(out.relative_to(ROOT))],
                'environment': {'software': [{'name': 'CPython', 'version': platform.python_version()}], 'hardware': platform.platform()},
                'mathematics': {'assertion_tested': 'Exact affine capacity/barrier extrema and fixed/guaranteed/possible quantifiers agree with independent simple-path envelope algebra',
                                'coefficient_domain': 'Q, exact fractions and graph combinatorics',
                                'conventions': 'Finite static undirected graphs; all path vertices include endpoints; closed rational one-parameter interval',
                                'inputs': files, 'bounds': dict(COUNTS), 'non_claims': payload['nonclaims']},
                'randomness': {'used': True, 'generator': 'Python random.Random MT19937', 'seed': seed},
                'run': {'started_at': started, 'runtime_seconds': time.perf_counter()-clock, 'status': 'completed', 'exit_status': 0,
                        'resource_limits': {'memory_bytes': None, 'cpu_seconds': None, 'max_cores': None, 'max_threads': None}},
                'input_artifacts': [{'path': p, 'sha256': before[p], 'sha256_after': after[p]} for p in files],
                'outputs': [{'path': str(output.relative_to(ROOT)).replace('\\', '/'), 'sha256': hashlib.sha256(output.read_bytes()).hexdigest(), 'kind': 'result'}],
                'checks': [{'name': 'source-hashes-unchanged', 'passed': True}], 'result': 'completed',
                'declared_results': [str(output.relative_to(ROOT)).replace('\\', '/')], 'residual_risks': ['No hard resource limits; finite tests supplement separate proofs.']}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'counts': dict(COUNTS), 'seconds': manifest['run']['runtime_seconds']}), flush=True)


if __name__ == '__main__':
    main()
