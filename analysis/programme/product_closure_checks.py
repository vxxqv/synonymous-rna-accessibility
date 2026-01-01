"""Bounded exhaustive falsification of additive Cartesian product closure.

The oracle constructs each small product explicitly and uses ordinary BFS.
Production code never constructs that product. Finite checks do not replace
the proof or establish novelty over existing additive merge-tree composition.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from heapq import heappop, heappush
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.product_closure import (component_closure, connection_witness,
    exact_product_bottleneck, product_connected)

P2 = ((1,), (0,))
P3 = ((1,), (0, 2), (1,))


def explicit_product(graphs, scores):
    vertices = list(product(*(range(len(g)) for g in graphs)))
    indices = {state: i for i, state in enumerate(vertices)}
    adjacency = []
    for state in vertices:
        neighbors = []
        for i, graph in enumerate(graphs):
            for v in graph[state[i]]:
                changed = list(state)
                changed[i] = v
                neighbors.append(indices[tuple(changed)])
        adjacency.append(tuple(sorted(neighbors)))
    values = [sum(scores[i][state[i]] for i in range(len(graphs))) for state in vertices]
    return vertices, adjacency, values


def bfs_components(adjacency, values, threshold):
    labels = [-1]*len(adjacency)
    components = []
    for start in range(len(adjacency)):
        if values[start] < threshold or labels[start] != -1:
            continue
        label = len(components)
        reached = {start}
        labels[start] = label
        queue = deque([start])
        while queue:
            for neighbor in adjacency[queue.popleft()]:
                if values[neighbor] >= threshold and labels[neighbor] == -1:
                    labels[neighbor] = label
                    reached.add(neighbor)
                    queue.append(neighbor)
        components.append(reached)
    return labels, components


def cyclic_reference(graphs, scores, start, threshold):
    reached = [{v} for v in start]
    peaks = [scores[i][v] for i, v in enumerate(start)]
    while True:
        changed = False
        for i, graph in enumerate(graphs):
            cutoff = threshold-sum(peaks[j] for j in range(len(graphs)) if j != i)
            component = {start[i]}
            queue = deque([start[i]])
            while queue:
                for v in graph[queue.popleft()]:
                    if scores[i][v] >= cutoff and v not in component:
                        component.add(v)
                        queue.append(v)
            if not reached[i] <= component:
                raise AssertionError('Cyclic coordinate closure was not monotone')
            changed |= reached[i] != component
            reached[i] = component
            peaks[i] = max(scores[i][v] for v in component)
        if not changed:
            return tuple(frozenset(a) for a in reached)


def check_moves(graphs, scores, start, threshold, moves, finish):
    state = list(start)
    if sum(scores[i][state[i]] for i in range(len(graphs))) < threshold:
        raise AssertionError('Witness begins below threshold')
    for i, v in moves:
        if v not in graphs[i][state[i]]:
            raise AssertionError('Witness contains a nonedge')
        state[i] = v
        if sum(scores[j][state[j]] for j in range(len(graphs))) < threshold:
            raise AssertionError('Witness crosses below threshold')
    if tuple(state) != tuple(finish):
        raise AssertionError('Witness ends at the wrong state')


def check_dataset(graphs, scores, thresholds, count):
    vertices, adjacency, values = explicit_product(graphs, scores)
    length_bound = sum((len(g)-1)**2 for g in graphs)
    for threshold in thresholds:
        labels, components = bfs_components(adjacency, values, threshold)
        projection = [tuple(frozenset(vertices[u][i] for u in component)
                            for i in range(len(graphs))) for component in components]
        signatures = {}
        for index, state in enumerate(vertices):
            if labels[index] < 0:
                continue
            closure = component_closure(graphs, scores, state, threshold, witness=True)
            expected = projection[labels[index]]
            if closure.coordinates != expected:
                raise AssertionError(('Wrong coordinate projection', graphs, scores, state, threshold))
            if cyclic_reference(graphs, scores, state, threshold) != expected:
                raise AssertionError('Heap and cyclic algorithms disagree')
            if closure.coordinates in signatures and signatures[closure.coordinates] != labels[index]:
                raise AssertionError('Different components share the same projection signature')
            signatures[closure.coordinates] = labels[index]
            expected_anchor = tuple(min(a, key=lambda v: (-scores[i][v], v))
                                    for i, a in enumerate(expected))
            if closure.anchor != expected_anchor:
                raise AssertionError('Noncanonical peak anchor')
            if sum(closure.peaks) != max(values[u] for u in components[labels[index]]):
                raise AssertionError('Joint peak fails to maximize score over its component')
            check_moves(graphs, scores, state, threshold, closure.moves, closure.anchor)
            if len(closure.moves) > length_bound:
                raise AssertionError('Constructive witness exceeds the polynomial bound')
            n = sum(map(len, graphs))
            m2 = sum(len(ns) for graph in graphs for ns in graph)
            if closure.statistics['admitted_vertices'] > n-len(graphs):
                raise AssertionError('A factor vertex was admitted more than once')
            if closure.statistics['adjacency_entries_scanned'] > m2:
                raise AssertionError('A factor adjacency list was repeatedly scanned')
            count['closure_queries'] += 1
            count['witness_moves_checked'] += len(closure.moves)
            count['largest_anchor_witness'] = max(count['largest_anchor_witness'], len(closure.moves))
        count['threshold_worlds'] += 1
        count['explicit_components'] += len(components)
        retained = sum(label >= 0 for label in labels)
        count['ordered_pairs_covered_by_partition_equality'] += retained**2
    count['factor_score_datasets'] += 1


def every_binary_factor():
    result = []
    for n in range(1, 4):
        edges = list(combinations(range(n), 2))
        for presence in product((0, 1), repeat=len(edges)):
            graph = [[] for _ in range(n)]
            for bit, (u, v) in zip(presence, edges):
                if bit:
                    graph[u].append(v)
                    graph[v].append(u)
            graph = tuple(tuple(ns) for ns in graph)
            for scores in product((0, 1), repeat=n):
                result.append((graph, scores))
    return result


def oracle_widest(adjacency, values, a, b):
    best = {a: values[a]}
    queue = [(-values[a], a)]
    while queue:
        neg, u = heappop(queue)
        if -neg != best[u]:
            continue
        if u == b:
            return -neg
        for v in adjacency[u]:
            value = min(-neg, values[v])
            if v not in best or value > best[v]:
                best[v] = value
                heappush(queue, (-value, v))
    return None


def rational_checks():
    count = Counter()
    alphabet = (Fraction(-2, 3), Fraction(7, 11))
    for left, right in product(product(alphabet, repeat=3), repeat=2):
        graphs = (P3, P3)
        scores = (left, right)
        vertices, adjacency, values = explicit_product(graphs, scores)
        for a, x in enumerate(vertices):
            for b, y in enumerate(vertices):
                actual = exact_product_bottleneck(graphs, scores, x, y)
                expected = oracle_widest(adjacency, values, a, b)
                if actual != expected:
                    raise AssertionError(('Rational widest capacity mismatch', scores, x, y))
                count['exact_bottleneck_queries'] += 1
                if not product_connected(graphs, scores, x, y, actual):
                    raise AssertionError('Widest threshold is not attainable')
                if product_connected(graphs, scores, x, y, actual+Fraction(1, 33)):
                    raise AssertionError('Above-widest threshold remains connected')
                moves = connection_witness(graphs, scores, x, y, actual)
                if moves is None:
                    raise AssertionError('Connected endpoints have no witness')
                check_moves(graphs, scores, x, actual, moves, y)
                if len(moves) > 2*sum((len(g)-1)**2 for g in graphs):
                    raise AssertionError('Endpoint witness exceeds bound')
                count['endpoint_witnesses_checked'] += 1
    # Disconnected factors and thresholds with a separate denominator.
    disconnected = (((), ()), P2)
    if exact_product_bottleneck(disconnected, ((1, 2), (3, 4)), (0, 0), (1, 1)) is not None:
        raise AssertionError('Baseline disconnection was lost')
    if product_connected((P2,), ((Fraction(1, 3), Fraction(2, 3)),), (0,), (1,), Fraction(1, 2)):
        raise AssertionError('A threshold with a new denominator was rounded')
    huge = Fraction(1, 2**101+3)
    if exact_product_bottleneck((P3,), ((huge, 0, 2*huge),), (0,), (2,)) != 0:
        raise AssertionError('Large rational denominator failed')
    count['boundary_checks'] = 3
    return dict(count)


def counterexample():
    graphs = (P3, P3)
    scores = ((2, 0, 1), (2, 0, 1))
    x, y, threshold = (0, 0), (2, 2), 2
    left = component_closure(graphs, scores, x, threshold)
    right = component_closure(graphs, scores, y, threshold)
    if not all(y[i] in left.coordinates[i] for i in range(2)):
        raise AssertionError('Projection counterexample is not inside the box')
    if product_connected(graphs, scores, x, y, threshold):
        raise AssertionError('Projection containment was mistaken for connectivity')
    return {'factor_scores': scores, 'threshold': threshold, 'x': x, 'y': y,
            'x_projection_sets': [sorted(a) for a in left.coordinates],
            'y_projection_sets': [sorted(a) for a in right.coordinates],
            'retained_y_inside_x_projection_product': True, 'connected': False}


def run(output):
    if output.exists():
        raise FileExistsError('A new output directory is required')
    started = perf_counter()
    results = {}
    factors = every_binary_factor()
    count = Counter()
    for left, right in combinations_with_replacement(factors, 2):
        check_dataset((left[0], right[0]), (left[1], right[1]), range(-1, 4), count)
    results['all_two_factors_up_to_three_vertices_binary_scores'] = dict(count)
    print(json.dumps({'completed': list(results)[-1], 'counts': dict(count)}), flush=True)
    count = Counter()
    for left, right in combinations_with_replacement(list(product((0, 1, 2), repeat=3)), 2):
        check_dataset((P3, P3), (left, right), range(-1, 6), count)
    results['two_P3_factors_ternary_scores'] = dict(count)
    print(json.dumps({'completed': list(results)[-1], 'counts': dict(count)}), flush=True)
    count = Counter()
    reversal_representatives = [s for s in product((0, 1, 2), repeat=3) if s <= s[::-1]]
    for scores in combinations_with_replacement(reversal_representatives, 3):
        check_dataset((P3,)*3, scores, range(-1, 8), count)
    results['three_P3_factors_ternary_scores'] = dict(count)
    print(json.dumps({'completed': list(results)[-1], 'counts': dict(count)}), flush=True)
    count = Counter()
    for scores in combinations_with_replacement(list(product((0, 1), repeat=2)), 4):
        check_dataset((P2,)*4, scores, range(-1, 6), count)
    results['four_P2_factors_binary_scores'] = dict(count)
    results['rational_capacity_and_endpoint_witnesses'] = rational_checks()
    results['projection_product_counterexample'] = counterexample()
    results['runtime_seconds'] = perf_counter()-started
    results['software'] = {'python': platform.python_version()}
    results['input_sha256'] = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in (
        'src/product_closure.py', 'analysis/programme/product_closure_checks.py',
        'analysis/programme/product_closure_contract_v1.json')}
    results['status'] = 'All stated finite assertions passed; proof and novelty require separate review'
    output.mkdir(parents=True, exist_ok=False)
    (output/'checks.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(results), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
