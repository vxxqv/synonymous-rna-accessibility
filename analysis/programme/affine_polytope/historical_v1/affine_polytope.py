"""Exact identification for rational affine scores on a bounded parameter set.

The set is an explicit box intersected with optional closed half spaces.
Candidate intersections are exhaustive within a declared resource cap. No
probabilities, empirical outcomes, fitted covariance or statistical coverage.
"""
from fractions import Fraction as F
from itertools import combinations
import heapq
from math import comb


def frac(x):
    if isinstance(x, bool) or not isinstance(x, (int, str, F)):
        raise ValueError('Use integer, rational string or Fraction, never float')
    return F(x)


def solve_square(rows, rhs):
    """Return the unique exact intersection, or None for deficient rank."""
    d = len(rows)
    if not d:
        return ()
    matrix = [list(row) + [value] for row, value in zip(rows, rhs)]
    for j in range(d):
        pivot = next((i for i in range(j, d) if matrix[i][j]), None)
        if pivot is None:
            return None
        matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
        divisor = matrix[j][j]
        matrix[j] = [x / divisor for x in matrix[j]]
        for i in range(d):
            if i != j and matrix[i][j]:
                multiplier = matrix[i][j]
                matrix[i] = [x - multiplier*y for x, y in zip(matrix[i], matrix[j])]
    return tuple(matrix[i][-1] for i in range(d))


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), F(0))


def adjacency(n, edges):
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('Positive integer vertex count required')
    graph = [set() for _ in range(n)]
    for u, v in edges:
        if not all(isinstance(w, int) and not isinstance(w, bool) and 0 <= w < n for w in (u, v)) or u == v:
            raise ValueError('Invalid edge')
        graph[u].add(v)
        graph[v].add(u)
    return [tuple(sorted(row)) for row in graph]


def widest(graph, scores, s, t):
    best = [None]*len(graph)
    predecessor = [None]*len(graph)
    done = set()
    best[s] = scores[s]
    queue = [(-scores[s], s)]
    while queue:
        negative, u = heapq.heappop(queue)
        if u in done or -negative != best[u]:
            continue
        done.add(u)
        if u == t:
            path = [t]
            while path[-1] != s:
                path.append(predecessor[path[-1]])
            return best[t], tuple(reversed(path))
        for v in graph[u]:
            if v not in done:
                value = min(best[u], scores[v])
                if best[v] is None or value > best[v]:
                    best[v], predecessor[v] = value, u
                    heapq.heappush(queue, (-value, v))
    return None, None


def prepare(case):
    n, s, t = case['n'], case['s'], case['t']
    graph = adjacency(n, case['edges'])
    if any(not isinstance(v, int) or isinstance(v, bool) or not 0 <= v < n for v in (s, t)):
        raise ValueError('Invalid endpoint')
    box = tuple((frac(lo), frac(hi)) for lo, hi in case['box'])
    if any(lo > hi for lo, hi in box):
        raise ValueError('Reversed parameter bound')
    d = len(box)
    scores = tuple(tuple(frac(x) for x in row) for row in case['coefficients'])
    if len(scores) != n or any(len(row) != d+1 for row in scores):
        raise ValueError('One intercept and d slopes required per vertex')
    constraints = []
    for j, (lo, hi) in enumerate(box):
        normal = tuple(F(int(i == j)) for i in range(d))
        constraints += [(normal, hi), (tuple(-x for x in normal), -lo)]
    for row, rhs in case.get('inequalities', []):
        row = tuple(frac(x) for x in row)
        if len(row) != d:
            raise ValueError('Invalid half space dimension')
        constraints.append((row, frac(rhs)))
    return graph, scores, box, tuple(constraints)


def arrangement_vertices(scores, box, constraints, max_intersections=200000):
    d = len(box)
    hyperplanes = set()
    boundaries = list(constraints)
    boundaries += [(tuple(x-y for x, y in zip(a[1:], b[1:])), b[0]-a[0])
                   for a, b in combinations(scores, 2)]
    for row, rhs in boundaries:
        pivot = next((x for x in row if x), None)
        if pivot is not None:
            hyperplanes.add((tuple(x/pivot for x in row), rhs/pivot))
    hyperplanes = sorted(hyperplanes)
    total = comb(len(hyperplanes), d) if len(hyperplanes) >= d else 0
    if total > max_intersections:
        raise RuntimeError(f'Resource cap: {total} intersections exceed {max_intersections}; no partial result')
    points = set()
    independent = 0
    for picked in combinations(hyperplanes, d):
        point = solve_square([x[0] for x in picked], [x[1] for x in picked])
        if point is not None:
            independent += 1
            if all(dot(row, point) <= rhs for row, rhs in constraints):
                points.add(point)
    return tuple(sorted(points)), {'hyperplanes':len(hyperplanes), 'candidate_intersections':total,
                                  'independent_intersections':independent, 'retained_points':len(points)}


def analyze(case, max_intersections=200000):
    graph, coefficients, box, constraints = prepare(case)
    points, counters = arrangement_vertices(coefficients, box, constraints, max_intersections)
    if not points:
        return {'id':case['id'], 'status':'empty_parameter_set', 'counters':counters}
    s, t = case['s'], case['t']
    values = [tuple(row[0] + dot(row[1:], point) for row in coefficients) for point in points]
    lower = tuple(min(row[i] for row in values) for i in range(len(graph)))
    upper = tuple(max(row[i] for row in values) for i in range(len(graph)))
    fixed, fixed_path = widest(graph, lower, s, t)
    if fixed is None:
        return {'id':case['id'], 'status':'disconnected', 'counters':counters,
                'capacity_extended_value':'-infinity', 'barrier_extended_value':'+infinity'}
    evaluations = []
    for point, scores in zip(points, values):
        capacity, path = widest(graph, scores, s, t)
        evaluations.append({'theta':point, 'capacity':capacity,
                            'barrier':min(scores[s], scores[t])-capacity, 'path':path})
    extrema = {}
    for field in ('capacity', 'barrier'):
        for label, selector in [('minimum', min), ('maximum', max)]:
            witness = selector(evaluations, key=lambda row:row[field])
            extrema[f'{field}_{label}'] = witness[field]
            extrema[f'{field}_{label}_witness'] = witness
    return {'id':case['id'], 'status':'connected', 'counters':counters, 'fixed':fixed,
            'fixed_path':fixed_path, 'lower_field':lower, 'upper_field':upper,
            'guaranteed':extrema['capacity_minimum'], 'possible':extrema['capacity_maximum'],
            **extrema}
