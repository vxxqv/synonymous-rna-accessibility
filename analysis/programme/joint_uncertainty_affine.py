"""Exact one-parameter affine uncertainty on finite, static undirected graphs.

Scores and the closed parameter interval are rational. No probabilities,
measurement calibration, biological data or optimization interventions enter.
Disconnected endpoint capacities use None, representing extended value -infinity;
their barriers are nonfinite and are reported by explicit status, not as zero.
"""
from fractions import Fraction
import heapq
import itertools


def rational(value):
    if not isinstance(value, (int, str, Fraction)) or isinstance(value, bool):
        raise ValueError('exact integer, rational string or Fraction required')
    return Fraction(value)


def graph_from_edges(n, edges):
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('positive integer vertex count required')
    graph = [set() for _ in range(n)]
    for edge in edges:
        if len(edge) != 2:
            raise ValueError('two distinct valid endpoints per edge required')
        i, j = edge
        if not isinstance(i, int) or not isinstance(j, int) or i == j or not (0 <= i < n and 0 <= j < n):
            raise ValueError('invalid graph edge')
        graph[i].add(j)
        graph[j].add(i)
    return [tuple(sorted(neighbors)) for neighbors in graph]


def widest_path(graph, scores, source, target):
    """Exact maximum vertex bottleneck, including both endpoints, with witness."""
    n = len(graph)
    if len(scores) != n or not (0 <= source < n and 0 <= target < n):
        raise ValueError('valid scores and distinguished endpoints required')
    best, previous, settled = [None]*n, [None]*n, set()
    best[source] = scores[source]
    queue = [(-scores[source], source)]
    while queue:
        negative, u = heapq.heappop(queue)
        if u in settled or -negative != best[u]:
            continue
        settled.add(u)
        if u == target:
            path = [u]
            while path[-1] != source:
                path.append(previous[path[-1]])
            return best[u], tuple(reversed(path))
        for v in graph[u]:
            if v in settled:
                continue
            candidate = min(best[u], scores[v])
            if best[v] is None or candidate > best[v]:
                best[v], previous[v] = candidate, u
                heapq.heappush(queue, (-candidate, v))
    return None, None


def score_at(coefficients, theta):
    return tuple(a + b*theta for a, b in coefficients)


def affine_breakpoints(coefficients, lower, upper):
    points = {lower, upper}
    for (a, b), (c, d) in itertools.combinations(coefficients, 2):
        if b != d:
            theta = (c-a)/(b-d)
            if lower <= theta <= upper:
                points.add(theta)
    return tuple(sorted(points))


def analyze(n, edges, coefficients, source, target, lower, upper):
    """Exact capacity/barrier ranges and fixed/guaranteed/possible values.

    Each open-cell certificate gives one path optimal throughout that cell,
    its affine capacity and the affine endpoint-relative barrier. Extrema
    include interval boundaries and all pairwise affine-score intersections.
    """
    graph = graph_from_edges(n, edges)
    if len(coefficients) != n or not (0 <= source < n and 0 <= target < n):
        raise ValueError('one affine score per vertex and valid endpoints required')
    coefficients = tuple((rational(a), rational(b)) for a, b in coefficients)
    lower, upper = rational(lower), rational(upper)
    if lower > upper:
        raise ValueError('closed nonempty interval required')
    low_scores = tuple(min(a+b*lower, a+b*upper) for a, b in coefficients)
    upper_scores = tuple(max(a+b*lower, a+b*upper) for a, b in coefficients)
    fixed, fixed_path = widest_path(graph, low_scores, source, target)
    box_possible, _ = widest_path(graph, upper_scores, source, target)
    if fixed is None:
        return {'status': 'disconnected', 'capacity_extended_value': '-infinity',
                'barrier_extended_value': '+infinity', 'fixed': None, 'guaranteed': None,
                'possible': None, 'capacity_range': None, 'barrier_range': None,
                'breakpoints': (), 'cells': (), 'fixed_path': None}
    breakpoints = affine_breakpoints(coefficients, lower, upper)
    evaluations = []
    for theta in breakpoints:
        scores = score_at(coefficients, theta)
        capacity, path = widest_path(graph, scores, source, target)
        endpoint_minimum = min(scores[source], scores[target])
        evaluations.append({'theta': theta, 'capacity': capacity,
                            'barrier': endpoint_minimum-capacity, 'path': path})
    cells = []
    for left, right in zip(breakpoints, breakpoints[1:]):
        midpoint = (left+right)/2
        scores = score_at(coefficients, midpoint)
        _, path = widest_path(graph, scores, source, target)
        bottleneck = min(path, key=lambda v: scores[v])
        endpoint = min((source, target), key=lambda v: scores[v])
        capacity_line = coefficients[bottleneck]
        barrier_line = tuple(x-y for x, y in zip(coefficients[endpoint], capacity_line))
        cells.append({'lower': left, 'upper': right, 'capacity_line': capacity_line,
                      'barrier_line': barrier_line, 'path': path,
                      'capacity_vertex': bottleneck, 'endpoint_vertex': endpoint})
    min_capacity = min(evaluations, key=lambda row: row['capacity'])
    max_capacity = max(evaluations, key=lambda row: row['capacity'])
    min_barrier = min(evaluations, key=lambda row: row['barrier'])
    max_barrier = max(evaluations, key=lambda row: row['barrier'])
    # Exact full-box barrier extrema use mixed endpoint/interior corners.
    barrier_high_scores = list(low_scores)
    barrier_low_scores = list(upper_scores)
    for v in (source, target):
        barrier_high_scores[v] = upper_scores[v]
        barrier_low_scores[v] = low_scores[v]
    high_capacity, _ = widest_path(graph, barrier_high_scores, source, target)
    low_capacity, _ = widest_path(graph, barrier_low_scores, source, target)
    box_barrier_range = (min(barrier_low_scores[source], barrier_low_scores[target])-low_capacity,
                         min(barrier_high_scores[source], barrier_high_scores[target])-high_capacity)
    return {'status': 'connected', 'fixed': fixed, 'fixed_path': fixed_path,
            'guaranteed': min_capacity['capacity'], 'possible': max_capacity['capacity'],
            'capacity_range': (min_capacity['capacity'], max_capacity['capacity']),
            'barrier_range': (min_barrier['barrier'], max_barrier['barrier']),
            'capacity_minimizer': min_capacity['theta'], 'capacity_maximizer': max_capacity['theta'],
            'barrier_minimizer': min_barrier['theta'], 'barrier_maximizer': max_barrier['theta'],
            'box_guaranteed': fixed, 'box_possible': box_possible, 'box_barrier_range': box_barrier_range,
            'lower_field': low_scores, 'upper_field': upper_scores,
            'breakpoints': breakpoints, 'evaluations': evaluations, 'cells': cells}
