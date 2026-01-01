"""Exact threshold reachability for additive Cartesian graph products.

Input factors are explicit finite simple undirected graphs. Scores and thresholds
are integers or fractions.Fraction. The product is never materialized here.
Composition of additive merge trees is established prior work; this module
implements a direct component query with a constructive certificate.
"""
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from heapq import heappop, heappush
from math import lcm


@dataclass(frozen=True)
class Closure:
    coordinates: tuple
    peaks: tuple
    anchor: tuple
    moves: tuple
    statistics: dict


def _fraction(value):
    if not isinstance(value, (int, Fraction)):
        raise TypeError('Exact integers or Fraction values are required')
    return Fraction(value)


def _prepare(graphs, scores, endpoints, threshold=None):
    if not graphs or len(graphs) != len(scores):
        raise ValueError('At least one factor and one score list per factor are required')
    normalized = []
    rational = []
    for graph, values in zip(graphs, scores):
        n = len(graph)
        if not n or len(values) != n:
            raise ValueError('Every factor must be nonempty with one score per vertex')
        sets = [set(neighbors) for neighbors in graph]
        for u, neighbors in enumerate(graph):
            if len(sets[u]) != len(neighbors) or u in sets[u]:
                raise ValueError('Expected simple factor graphs')
            if any(not isinstance(v, int) or not 0 <= v < n for v in neighbors):
                raise ValueError('Invalid factor vertex index')
        if any(u not in sets[v] for u, neighbors in enumerate(sets) for v in neighbors):
            raise ValueError('Expected undirected factor graphs')
        normalized.append(tuple(tuple(sorted(ns)) for ns in sets))
        rational.append(tuple(_fraction(value) for value in values))
    for endpoint in endpoints:
        if len(endpoint) != len(graphs) or any(
                not isinstance(v, int) or not 0 <= v < len(graphs[i])
                for i, v in enumerate(endpoint)):
            raise ValueError('Invalid product endpoint')
    rational_threshold = None if threshold is None else _fraction(threshold)
    denominator = 1
    for values in rational:
        for value in values:
            denominator = lcm(denominator, value.denominator)
    if rational_threshold is not None:
        denominator = lcm(denominator, rational_threshold.denominator)
    integer_scores = tuple(tuple(value.numerator*(denominator//value.denominator)
                                 for value in values) for values in rational)
    integer_threshold = (None if rational_threshold is None else
                         rational_threshold.numerator*(denominator//rational_threshold.denominator))
    return tuple(normalized), integer_scores, integer_threshold, denominator


def _tree_path(parents, start, finish):
    left = []
    u = start
    while u is not None:
        left.append(u)
        u = parents[u]
    left_position = {v: i for i, v in enumerate(left)}
    right = []
    v = finish
    while v not in left_position:
        right.append(v)
        v = parents[v]
    return left[:left_position[v]+1] + right[::-1]


def _closure(graphs, scores, start, threshold, witness=False):
    dimension = len(graphs)
    peaks = [scores[i][start[i]] for i in range(dimension)]
    surplus = sum(peaks)-threshold
    if surplus < 0:
        raise ValueError('The starting endpoint is below threshold')
    reached = [{vertex} for vertex in start]
    discovered = [{vertex} for vertex in start]
    parents = [{vertex: None} for vertex in start]
    frontier = [[] for _ in graphs]
    global_heap = []
    anchor = list(start)
    moves = []
    statistics = {'admitted_vertices': 0, 'factor_heap_pushes': 0,
                  'global_heap_pops': 0, 'adjacency_entries_scanned': 0,
                  'peak_anchor_changes': 0}

    def expose(i, u):
        for v in graphs[i][u]:
            statistics['adjacency_entries_scanned'] += 1
            if v not in discovered[i]:
                discovered[i].add(v)
                parents[i][v] = u
                heappush(frontier[i], (-scores[i][v], v))
                statistics['factor_heap_pushes'] += 1

    for i, vertex in enumerate(start):
        expose(i, vertex)
        if frontier[i]:
            heappush(global_heap, (peaks[i]+frontier[i][0][0], i))

    while global_heap and global_heap[0][0] <= surplus:
        deficit, i = heappop(global_heap)
        statistics['global_heap_pops'] += 1
        neg_score, vertex = heappop(frontier[i])
        if deficit != peaks[i]+neg_score:
            raise AssertionError('A factor frontier representative became stale')
        reached[i].add(vertex)
        statistics['admitted_vertices'] += 1
        value = -neg_score
        if value > peaks[i] or (value == peaks[i] and vertex < anchor[i]):
            if witness:
                moves.extend((i, v) for v in _tree_path(parents[i], anchor[i], vertex)[1:])
            surplus += value-peaks[i]
            peaks[i] = value
            anchor[i] = vertex
            statistics['peak_anchor_changes'] += 1
        expose(i, vertex)
        if frontier[i]:
            heappush(global_heap, (peaks[i]+frontier[i][0][0], i))

    return Closure(tuple(frozenset(vertices) for vertices in reached), tuple(peaks),
                   tuple(anchor), tuple(moves), statistics)


def component_closure(graphs, scores, start, threshold, witness=False):
    """Return exact coordinate projections and a reachable canonical peak anchor.

    The reachable component need not equal the retained part of the Cartesian
    product of these projections. Query another endpoint with its own closure.
    """
    graph_data, integer_scores, integer_threshold, denominator = _prepare(
        graphs, scores, (start,), threshold)
    result = _closure(graph_data, integer_scores, tuple(start), integer_threshold, witness)
    return Closure(result.coordinates, tuple(Fraction(p, denominator) for p in result.peaks),
                   result.anchor, result.moves, result.statistics)


def product_connected(graphs, scores, x, y, threshold):
    """Exact connectivity of two endpoints in the induced superlevel product."""
    graph_data, values, level, _ = _prepare(graphs, scores, (x, y), threshold)
    if min(sum(values[i][x[i]] for i in range(len(graphs))),
           sum(values[i][y[i]] for i in range(len(graphs)))) < level:
        return False
    first = _closure(graph_data, values, tuple(x), level)
    second = _closure(graph_data, values, tuple(y), level)
    return first.coordinates == second.coordinates


def connection_witness(graphs, scores, x, y, threshold):
    """Return coordinate moves witnessing connectivity, or None if disconnected."""
    graph_data, values, level, _ = _prepare(graphs, scores, (x, y), threshold)
    if min(sum(values[i][x[i]] for i in range(len(graphs))),
           sum(values[i][y[i]] for i in range(len(graphs)))) < level:
        return None
    first = _closure(graph_data, values, tuple(x), level, True)
    second = _closure(graph_data, values, tuple(y), level, True)
    if first.coordinates != second.coordinates:
        return None
    if first.anchor != second.anchor:
        raise AssertionError('Equal projections did not yield a common canonical anchor')
    old_moves = []
    state = list(y)
    for i, vertex in second.moves:
        old_moves.append((i, state[i]))
        state[i] = vertex
    return first.moves+tuple(reversed(old_moves))


def exact_product_bottleneck(graphs, scores, x, y):
    """Exact widest capacity by integer threshold search without product expansion.

    Return None for different baseline Cartesian components. The bit complexity
    includes common denominator scaling and logarithmic numeric range search.
    """
    graph_data, values, _, denominator = _prepare(graphs, scores, (x, y))
    minima = []
    for i, graph in enumerate(graph_data):
        reached = {x[i]}
        queue = deque([x[i]])
        while queue:
            for neighbor in graph[queue.popleft()]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    queue.append(neighbor)
        if y[i] not in reached:
            return None
        minima.append(min(values[i][v] for v in reached))
    low = sum(minima)
    high = min(sum(values[i][x[i]] for i in range(len(graphs))),
               sum(values[i][y[i]] for i in range(len(graphs))))
    while low < high:
        middle = (low+high+1)//2
        left = _closure(graph_data, values, tuple(x), middle)
        right = _closure(graph_data, values, tuple(y), middle)
        if left.coordinates == right.coordinates:
            low = middle
        else:
            high = middle-1
    return Fraction(low, denominator)
