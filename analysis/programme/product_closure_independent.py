"""Fresh bounded adversarial check of held product-closure source.

Independent proof obligations are exercised using explicit product BFS and a
separate whole-factor threshold closure. Source/proof remain read-only.
"""
import argparse
from collections import Counter, deque, defaultdict
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

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO
sys.path.insert(0, str(REPO))
from src import product_closure as tested

COUNT = Counter()


def explicit(graphs, scores):
    states = tuple(itertools.product(*(range(len(g)) for g in graphs)))
    number = {state: i for i, state in enumerate(states)}
    edges = [[] for _ in states]
    values = [sum(scores[i][state[i]] for i in range(len(graphs))) for state in states]
    for u, state in enumerate(states):
        for i, graph in enumerate(graphs):
            for v in graph[state[i]]:
                changed = list(state)
                changed[i] = v
                edges[u].append(number[tuple(changed)])
    return states, number, edges, values


def components(graph, values, threshold):
    label = {}
    groups = []
    for start in range(len(graph)):
        if start in label or values[start] < threshold:
            continue
        identifier = len(groups)
        label[start] = identifier
        group, queue = {start}, deque([start])
        while queue:
            for v in graph[queue.popleft()]:
                if values[v] >= threshold and v not in label:
                    label[v] = identifier
                    group.add(v)
                    queue.append(v)
        groups.append(group)
    return label, groups


def tree_route(graph, allowed, source, target):
    parent, queue = {source: None}, deque([source])
    while target not in parent:
        u = queue.popleft()
        for v in graph[u]:
            if v in allowed and v not in parent:
                parent[v] = u
                queue.append(v)
    path = [target]
    while path[-1] != source:
        path.append(parent[path[-1]])
    return list(reversed(path))


def independent_anchor_closure(graphs, scores, start, threshold, true_states):
    sets = [{v} for v in start]
    anchor = list(start)
    peaks = [scores[i][v] for i, v in enumerate(start)]
    while True:
        change = False
        for i, graph in enumerate(graphs):
            cutoff = threshold-sum(value for j, value in enumerate(peaks) if j != i)
            reached, queue = {start[i]}, deque([start[i]])
            assert scores[i][start[i]] >= cutoff
            while queue:
                for v in graph[queue.popleft()]:
                    if v not in reached and scores[i][v] >= cutoff:
                        reached.add(v)
                        queue.append(v)
            assert sets[i] <= reached
            change |= sets[i] != reached
            sets[i] = reached
            new_anchor = min(reached, key=lambda v: (-scores[i][v], v))
            for v in tree_route(graph, reached, anchor[i], new_anchor)[1:]:
                anchor[i] = v
                assert tuple(anchor) in true_states
                assert sum(scores[j][anchor[j]] for j in range(len(graphs))) >= threshold
                COUNT['independent_anchor_moves'] += 1
            peaks[i] = scores[i][new_anchor]
            assert tuple(anchor) in true_states
            assert sum(scores[j][anchor[j]] for j in range(len(graphs))) == sum(peaks)
            for j, group in enumerate(sets):
                assert all(scores[j][v]+sum(p for k, p in enumerate(peaks) if k != j) >= threshold for v in group)
                assert group <= {state[j] for state in true_states}
            COUNT['intermediate_anchor_invariants'] += 1
        if not change:
            return tuple(frozenset(group) for group in sets), tuple(anchor)


def replay(graphs, scores, source, moves, threshold, target):
    state = list(source)
    assert sum(scores[i][v] for i, v in enumerate(state)) >= threshold
    for i, v in moves:
        assert v in graphs[i][state[i]]
        state[i] = v
        assert sum(scores[j][w] for j, w in enumerate(state)) >= threshold
        COUNT['production_witness_moves'] += 1
    assert tuple(state) == tuple(target)


def check_model(graphs, scores, capacity=False):
    scores = tuple(tuple(F(v) for v in values) for values in scores)
    states, number, adjacency, values = explicit(graphs, scores)
    levels = sorted(set(values))
    thresholds = sorted(set(levels+[levels[0]-1, levels[-1]+1]+[(a+b)/2 for a, b in zip(levels, levels[1:])]))
    single_bound = sum((len(g)-1)**2 for g in graphs)
    maximum_capacity = [[None]*len(states) for _ in states]
    for threshold in thresholds:
        labels, groups = components(adjacency, values, threshold)
        true_signatures = [tuple(frozenset(states[u][i] for u in group) for i in range(len(graphs))) for group in groups]
        predicted, actual = defaultdict(set), defaultdict(set)
        for u, state in enumerate(states):
            if u not in labels:
                assert tested.product_connected(graphs, scores, state, state, threshold) is False
                assert tested.connection_witness(graphs, scores, state, state, threshold) is None
                COUNT['absent_endpoint_checks'] += 1
                continue
            truth = {states[v] for v in groups[labels[u]]}
            closure = tested.component_closure(graphs, scores, state, threshold, witness=True)
            assert closure.coordinates == true_signatures[labels[u]]
            predicted[closure.coordinates].add(u)
            actual[labels[u]].add(u)
            expected_anchor = tuple(min(group, key=lambda v: (-scores[i][v], v)) for i, group in enumerate(closure.coordinates))
            assert closure.anchor == expected_anchor
            assert tuple(closure.anchor) in truth
            assert sum(closure.peaks) == max(values[v] for v in groups[labels[u]])
            replay(graphs, scores, state, closure.moves, threshold, closure.anchor)
            assert len(closure.moves) <= single_bound
            stats = closure.statistics
            size = sum(map(len, graphs))
            assert stats['admitted_vertices'] <= size-len(graphs)
            assert stats['factor_heap_pushes'] <= size-len(graphs)
            assert stats['global_heap_pops'] == stats['admitted_vertices']
            assert stats['adjacency_entries_scanned'] <= sum(len(nb) for g in graphs for nb in g)
            independent, anchor = independent_anchor_closure(graphs, scores, state, threshold, truth)
            assert independent == closure.coordinates and anchor == closure.anchor
            COUNT['retained_start_closure_queries'] += 1
        assert {frozenset(group) for group in predicted.values()} == {frozenset(group) for group in actual.values()}
        COUNT['complete_partition_checks'] += 1
        # All pairs are classified by the full partition; replay one far pair
        # per retained component and one pair from distinct components.
        for group in groups:
            x, y = states[min(group)], states[max(group)]
            moves = tested.connection_witness(graphs, scores, x, y, threshold)
            assert moves is not None and len(moves) <= 2*single_bound
            replay(graphs, scores, x, moves, threshold, y)
            assert tested.product_connected(graphs, scores, x, y, threshold)
            if threshold in levels:
                for u, v in itertools.product(group, repeat=2):
                    old = maximum_capacity[u][v]
                    maximum_capacity[u][v] = threshold if old is None else max(old, threshold)
            COUNT['endpoint_witness_queries'] += 1
        if len(groups) > 1:
            x, y = states[min(groups[0])], states[min(groups[1])]
            assert not tested.product_connected(graphs, scores, x, y, threshold)
            assert tested.connection_witness(graphs, scores, x, y, threshold) is None
            COUNT['disconnected_endpoint_queries'] += 1
    if capacity:
        for u, x in enumerate(states):
            for v, y in enumerate(states):
                assert tested.exact_product_bottleneck(graphs, scores, x, y) == maximum_capacity[u][v]
                COUNT['exact_rational_capacity_queries'] += 1
    COUNT['models'] += 1
    COUNT['largest_product'] = max(COUNT['largest_product'], len(states))


def main():
    if not __debug__:
        raise RuntimeError('independent assertions must remain enabled')
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    out = (WORKSPACE/args.outdir).resolve()
    if WORKSPACE not in out.parents or out.exists():
        raise ValueError('fresh workspace output directory required')
    out.mkdir(parents=True)
    inputs = ['src/product_closure.py',
              'analysis/programme/product_closure_proof.md',
              'analysis/programme/product_closure_independent.py']
    before = {p: hashlib.sha256((WORKSPACE/p).read_bytes()).hexdigest() for p in inputs}
    started, clock = datetime.datetime.now(datetime.timezone.utc).isoformat(), time.perf_counter()
    p1, p2, p3 = ((),), ((1,), (0,)), ((1,), (0, 2), (1,))
    cycle5 = ((1, 4), (0, 2), (1, 3), (2, 4), (0, 3))
    path5 = ((1,), (0, 2), (1, 3), (2, 4), (3,))
    disconnected = ((1,), (0,), (3,), (2,))
    cases = [
        ((p1,), ((-3,),), True),
        ((cycle5,), ((F(-2, 3), F(7, 11), -2, 0, F(1, 2)),), True),
        ((p3, p3), ((3, 0, 1), (3, 0, 1)), True),
        ((path5, path5), ((4, -1, 2, -2, 6), (4, -2, 3, 0, 7)), False),
        ((p3, p3, p3), ((2, -1, 3), (3, 0, 1), (-1, -2, 2)), False),
        ((disconnected, p3), ((F(1, 2), -2, F(2, 3), 3), (0, -1, F(3, 5))), True),
        ((p2,)*5, ((0, 1), (1, -1), (0, 0), (-2, 3), (F(1, 2), F(-1, 3))), False),
        ((p2,), ((F(-1, 1000000007), F(1, 1000000009)),), True),
        ((p3,), ((F(1, 1000000007), F(-1, 1000000009), F(1, 1000000033)),), True)]
    for graphs, scores, capacity in cases:
        check_model(graphs, scores, capacity)
    seed = 2026102503
    rng = random.Random(seed)
    for k in range(24):
        d = 1+k % 4
        graphs, scores = [], []
        for _ in range(d):
            n = rng.randint(1, 4)
            graph = [[] for _ in range(n)]
            for i, j in itertools.combinations(range(n), 2):
                if rng.random() < .4:
                    graph[i].append(j)
                    graph[j].append(i)
            graphs.append(tuple(tuple(nb) for nb in graph))
            scores.append(tuple(F(rng.randint(-3, 3), rng.choice((1, 2, 3, 5))) for _ in range(n)))
        check_model(tuple(graphs), tuple(scores), capacity=d <= 2)
    # Explicitly falsify the stronger one-sided projection-membership claim.
    first = tested.component_closure((p3, p3), ((3, 0, 1), (3, 0, 1)), (0, 0), 2)
    second = tested.component_closure((p3, p3), ((3, 0, 1), (3, 0, 1)), (2, 2), 2)
    assert all(2 in coordinate for coordinate in first.coordinates)
    assert first.coordinates != second.coordinates
    assert not tested.product_connected((p3, p3), ((3, 0, 1), (3, 0, 1)), (0, 0), (2, 2), 2)
    payload = {'status': 'all independent checks passed', 'counts': dict(COUNT), 'seed': seed,
               'projection_product_counterexample': {'graphs': 'P3 x P3', 'scores': [[3, 0, 1], [3, 0, 1]], 'threshold': 2,
                 'x': [0, 0], 'y': [2, 2], 'x_projections': [sorted(a) for a in first.coordinates], 'y_projections': [sorted(a) for a in second.coordinates]},
               'scope': 'Abstract additive Cartesian products only; no biological datasets; finite checks supplement independent proof.'}
    result_path = out/'independent_results.json'
    result_path.write_text(json.dumps(payload, indent=2)+'\n')
    after = {p: hashlib.sha256((WORKSPACE/p).read_bytes()).hexdigest() for p in inputs}
    if before != after:
        raise RuntimeError('held source/proof changed during audit')
    manifest = {'schema_version': 2, 'claim_id': 'PRODUCT_CLOSURE_INDEPENDENT_REVIEW',
                'repository': {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(), 'dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True))},
                'command': [sys.executable, 'analysis/programme/product_closure_independent.py', '--outdir', str(out.relative_to(WORKSPACE))],
                'environment': {'software': [{'name': 'CPython', 'version': platform.python_version()}], 'hardware': platform.platform()},
                'mathematics': {'assertion_tested': 'Held closure source matches independent explicit product components, intermediate anchor invariants and rational widest capacities',
                                'coefficient_domain': 'Q and exact finite graph combinatorics', 'conventions': 'Nonempty undirected explicit factors; additive scores; inclusive superlevels; absent endpoints are not connected',
                                'inputs': inputs, 'bounds': dict(COUNT), 'non_claims': ['Finite tests do not replace the independent proof.', 'No algorithmic novelty claim.', 'No claim that projection containment implies component membership.']},
                'randomness': {'used': True, 'generator': 'Python MT19937', 'seed': seed},
                'run': {'started_at': started, 'runtime_seconds': time.perf_counter()-clock, 'status': 'completed', 'exit_status': 0, 'resource_limits': {'memory_bytes': None, 'cpu_seconds': None, 'max_cores': None, 'max_threads': None}},
                'input_artifacts': [{'path': p, 'sha256': before[p], 'sha256_after': after[p]} for p in inputs],
                'outputs': [{'path': str(result_path.relative_to(WORKSPACE)).replace('\\', '/'), 'sha256': hashlib.sha256(result_path.read_bytes()).hexdigest(), 'kind': 'result'}],
                'checks': [{'name': 'held-input-hashes-unchanged', 'passed': True}], 'result': 'completed',
                'residual_risks': ['No hard resource limits. Source priority comparison remains incomplete.'],
                'declared_results': [str(result_path.relative_to(WORKSPACE)).replace('\\', '/')]}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'counts': dict(COUNT), 'seconds': manifest['run']['runtime_seconds']}), flush=True)


if __name__ == '__main__':
    main()
