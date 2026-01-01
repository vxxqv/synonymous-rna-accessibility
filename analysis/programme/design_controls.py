"""Exact design-identification and endpoint controls. No RNA outcomes are used."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import argparse
import json


def cube(n):
    return list(product((0, 1), repeat=n))


def design(points, order):
    columns = [s for k in range(order + 1) for s in combinations(range(len(points[0])), k)]
    return [[F(all(x[i] for i in s)) for s in columns] for x in points]


def rank(matrix):
    a = [list(map(F, row)) for row in matrix]
    r = 0
    for c in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [v / p for v in a[r]]
        for j in range(len(a)):
            if j != r:
                p = a[j][c]
                a[j] = [v - p * w for v, w in zip(a[j], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def adjacency(points):
    return [[j for j, y in enumerate(points) if sum(a != b for a, b in zip(x, y)) == 1]
            for x in points]


def bottleneck(graph, scores, start, end):
    """Max-min Dijkstra. Ties are admitted together, with no jitter."""
    best = [None] * len(graph)
    best[start] = scores[start]
    remaining = set(range(len(graph)))
    while remaining:
        eligible = [i for i in remaining if best[i] is not None]
        if not eligible:
            return None
        u = max(eligible, key=lambda i: best[i])
        remaining.remove(u)
        if u == end:
            return best[u]
        for v in graph[u]:
            candidate = min(best[u], scores[v])
            if best[v] is None or candidate > best[v]:
                best[v] = candidate
    return None


def endpoint_lifetime(graph, scores, start, end):
    b = bottleneck(graph, scores, start, end)
    if b is None:
        raise ValueError('Endpoints must be in the same unfiltered component')
    return max(F(0), min(scores[start], scores[end]) - b)


def component_sizes(graph, active):
    remaining = set(active)
    sizes = []
    while remaining:
        first = min(remaining)
        remaining.remove(first)
        stack, size = [first], 0
        while stack:
            u = stack.pop()
            size += 1
            for v in graph[u]:
                if v in remaining:
                    remaining.remove(v)
                    stack.append(v)
        sizes.append(size)
    return sizes


def integral(graph, scores):
    critical = sorted(set([F(0), F(1)] + list(scores)))
    result = F(0)
    for lo, hi in zip(critical, critical[1:]):
        mid = (lo + hi) / 2
        sizes = component_sizes(graph, [i for i, q in enumerate(scores) if q >= mid])
        if sizes:
            result += (hi - lo) * (1 - F(sum(s*s for s in sizes), sum(sizes)**2))
    return result


def mobius(points, scores):
    values = dict(zip(points, scores))
    coefficients = {}
    for mask in points:
        selected = [i for i, x in enumerate(mask) if x]
        coefficients[''.join(map(str, mask))] = sum(
            (-1)**(len(selected) - sum(x)) * values[x]
            for x in points if all(not x[i] or mask[i] for i in range(len(x))))
    return coefficients


def exhaustive_bottleneck(graph, scores, start, end):
    """Independent threshold BFS oracle rather than path optimization."""
    for t in sorted(set(scores), reverse=True):
        if scores[start] < t or scores[end] < t:
            continue
        seen, stack = {start}, [start]
        while stack:
            for v in graph[stack.pop()]:
                if scores[v] >= t and v not in seen:
                    seen.add(v)
                    stack.append(v)
        if end in seen:
            return t
    return None


def run():
    q2, q3 = cube(2), cube(3)
    # The old pair states AU, GU, GC: intercept and independent nucleotide indicators.
    old = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
    ranks = {'old_three_state_additive': rank(old),
             'Q2_additive': rank(design(q2, 1)),
             'Q3_additive': rank(design(q3, 1)),
             'Q3_pairwise': rank(design(q3, 2)),
             'Q3_all_orders': rank(design(q3, 3))}
    assert list(ranks.values()) == [3, 3, 4, 7, 8]
    g = adjacency(q3)
    a, b = q3.index((0, 0, 0)), q3.index((1, 1, 0))
    fields = {
        'constant_negative': [F(1, 2) for x in q3],
        'additive_negative': [F(x[0] + 2*x[1] + 4*x[2], 7) for x in q3],
        'pairwise_barrier': [F(1 - (x[0] - x[1])**2) for x in q3],
        'parity_three_way_barrier': [F(sum(x) % 2 == 0) for x in q3],
        'cubic_interaction_without_barrier': [F(x[0]*x[1]*x[2]) for x in q3]
    }
    expected = {'constant_negative': (0, 0, 0), 'additive_negative': (0, 0, 0),
                'pairwise_barrier': (1, F(1, 2), 0),
                'parity_three_way_barrier': (1, F(3, 4), -4),
                'cubic_interaction_without_barrier': (0, 0, 1)}
    controls = []
    for name, scores in fields.items():
        lifetime, area = endpoint_lifetime(g, scores, a, b), integral(g, scores)
        coefficients = mobius(q3, scores)
        assert (lifetime, area, coefficients['111']) == expected[name]
        assert bottleneck(g, scores, a, b) == exhaustive_bottleneck(g, scores, a, b)
        controls.append({'name': name, 'scores': scores, 'endpoint_lifetime': lifetime,
                         'integrated_fragmentation': area, 'mobius_coefficients': coefficients})
    # Same change in one non-endpoint score. One repairs the chosen endpoint path;
    # the other is matched for degree, initial score, and distance from the first endpoint.
    base = fields['parity_three_way_barrier']
    targeted, matched = list(base), list(base)
    targeted[q3.index((1, 0, 0))] = 1
    matched[q3.index((0, 0, 1))] = 1
    assert endpoint_lifetime(g, targeted, a, b) == 0
    assert endpoint_lifetime(g, matched, a, b) == 1
    assert integral(g, targeted) == integral(g, matched)
    # Exhaustive small-field oracle check, including many exact ties.
    checked = 0
    for scores in product((F(0), F(1, 2), F(1)), repeat=8):
        assert bottleneck(g, scores, a, b) == exhaustive_bottleneck(g, scores, a, b)
        checked += 1
    # Q2: a barrier implies a quantitative departure from every additive field.
    # Exact sup-norm distance to the additive subspace is |Delta|/4.
    square_checks = 0
    g2 = adjacency(q2)
    for scores in product((F(0), F(1, 2), F(1)), repeat=4):
        delta = scores[0] - scores[1] - scores[2] + scores[3]
        closed = max(F(0), min(scores[0], scores[3]) - max(scores[1], scores[2]))
        assert endpoint_lifetime(g2, scores, 0, 3) == closed
        assert closed <= abs(delta)/2
        residual_sign = [1, -1, -1, 1]
        projected = [v - delta*s/4 for v, s in zip(scores, residual_sign)]
        assert projected[0] - projected[1] - projected[2] + projected[3] == 0
        assert max(abs(x-y) for x, y in zip(scores, projected)) == abs(delta)/4
        square_checks += 1
    return {'evidence_type': 'exact finite verification plus separate analytic derivation',
            'points': q3, 'fixed_endpoints': [q3[a], q3[b]], 'ranks': ranks,
            'residual_degrees_of_freedom': {'old': 0, 'Q2_additive': 1,
                                            'Q3_additive': 4, 'Q3_pairwise': 1},
            'controls': controls, 'endpoint_intervention': {
                'base_lifetime': 1, 'targeted_lifetime': 0, 'matched_lifetime': 1,
                'targeted_integral': integral(g, targeted), 'matched_integral': integral(g, matched),
                'interpretation': 'Global fragmentation alone cannot distinguish endpoint-specific rescue.'},
            'exhaustive_tied_Q3_fields': checked, 'exhaustive_Q2_fields': square_checks,
            'non_claims': ['Finite checks do not prove a general theorem.',
                           'The design can identify statistical interactions, not physical folding causality.',
                           'No new RNA partition functions or biological validation were performed.']}


def serialize(value):
    if isinstance(value, F):
        return str(value)
    raise TypeError(type(value).__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    path = Path(parser.parse_args().out)
    result = run()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, default=serialize) + '\n', encoding='utf-8')
    print(json.dumps({'ranks': result['ranks'], 'Q3_fields_checked': result['exhaustive_tied_Q3_fields'],
                      'synthetic_controls_passed': len(result['controls'])}))
