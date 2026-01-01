"""Invented exact interval checks only; no cPRIME outcomes or primary imports."""
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path


def bounds(intervals):
    if not intervals:
        raise ValueError('At least one marker is required.')
    for lo, hi in intervals:
        if lo > hi or lo == math.inf or hi == -math.inf:
            raise ValueError('Invalid finite-response interval.')
    if len(intervals) == 1:
        return F(0), F(0)
    minimum = max(F(0), max(x[0] for x in intervals)-min(x[1] for x in intervals))
    supremum = max(hi-lo for i, (_, hi) in enumerate(intervals)
                    for j, (lo, _) in enumerate(intervals) if i != j)
    return minimum, supremum


def graph_possible(cells):
    n = len(cells[0])
    assert all(len(cell) == n for cell in cells)
    arcs = {(i, j) for cell in cells for i in range(n) for j in range(n)
            if i != j and cell[i][1] < cell[j][0]}
    remaining = set(range(n))
    while remaining:
        roots = {j for j in remaining if not any(i in remaining and k == j for i, k in arcs)}
        if not roots:
            return False
        remaining -= roots
    return True


def possible_order_mask(intervals, grid):
    """Independent finite oracle: enumerate responses, then compatible orders."""
    orders = tuple(itertools.permutations(range(len(intervals))))
    domain = [[x for x in grid if lo <= x <= hi] for lo, hi in intervals]
    mask = 0
    for values in itertools.product(*domain):
        for index, order in enumerate(orders):
            if all(values[order[i]] <= values[order[i+1]] for i in range(len(order)-1)):
                mask |= 1 << index
    return mask


def main():
    grid = tuple(F(i, 2) for i in range(5))
    intervals = tuple((a, b) for a in grid for b in grid if a <= b)
    counts = {}
    cached = {}
    for n in (1, 2, 3):
        checked = 0
        patterns = list(itertools.product(intervals, repeat=n))
        masks = []
        for pattern in patterns:
            domain = [[x for x in grid if lo <= x <= hi] for lo, hi in pattern]
            realized = [max(v)-min(v) for v in itertools.product(*domain)]
            assert bounds(pattern) == (min(realized), max(realized)), pattern
            masks.append(possible_order_mask(pattern, grid))
            checked += 1
        counts[n] = checked
        cached[n] = patterns, masks
    pair_checks = 0
    patterns, masks = cached[2]
    for i, first in enumerate(patterns):
        for j, second in enumerate(patterns):
            assert graph_possible((first, second)) == bool(masks[i] & masks[j]), (first, second)
            pair_checks += 1
    n3_checks = 0
    patterns, masks = cached[3]
    # Fixed synthetic indices, not a claim of exhaustive two-cell 3-marker pairs.
    for i in range(2048):
        a, b = (i*137+19) % len(patterns), (i*251+73) % len(patterns)
        assert graph_possible((patterns[a], patterns[b])) == bool(masks[a] & masks[b])
        n3_checks += 1
    assert bounds(((F(0), F(10)), (F(4), F(6)))) == (0, 6)
    assert bounds(((F(0), F(10)),)) == (0, 0)
    assert bounds(((F(0), math.inf),)) == (0, 0)
    assert bounds(((F(0), math.inf), (F(1), F(2)))) == (0, math.inf)
    assert bounds(((F(0), math.inf), (F(1), F(2)), (F(4), F(5)))) == (2, math.inf)
    assert bounds(((-math.inf, math.inf), (F(-1), F(1)))) == (0, math.inf)
    assert bounds(((F(0), F(1)), (F(1), F(2)))) == (0, 2)
    assert bounds(((F(0), F(0)), (F(2), F(2)))) == (2, 2)
    # Positive raw values approach the unattained supremum 10.
    approximants = [F(10)-F(1, q) for q in (1, 2, 10, 100)]
    assert all(x < 10 for x in approximants) and approximants == sorted(approximants)
    assert not graph_possible((((F(1), F(2)), (F(3), F(4))), ((F(3), F(4)), (F(1), F(2)))))
    assert graph_possible((((F(1), F(2)), (F(2), F(3))), ((F(2), F(3)), (F(1), F(2)))))
    assert graph_possible((((F(1), F(2)), (F(3), F(4))), ((-math.inf, math.inf), (-math.inf, math.inf))))
    # Identical positive controls/combined states admit opposite isolated effects.
    for alpha in (F(-1, 2), F(1, 2)):
        mu, beta, eta = F(1), F(0), -alpha
        assert mu > 0 and mu+alpha > 0 and mu+beta == 1
        assert mu+alpha+beta+eta == 1
    # Bounded interactions, exact combined-minus-marker contrasts 1 and 2.
    alpha_set = [a for a in (F(i, 4) for i in range(-8, 17))
                 if abs(F(1)-a) <= 1 and abs(F(2)-a) <= 1]
    assert min(alpha_set) == 1 and max(alpha_set) == 2
    assert all(abs(F(1)-a) <= 1 and abs(F(2)-a) <= 1 for a in (F(1), F(2)))
    # With eta=0, differing exact adjusted contrasts are inconsistent.
    assert not any(F(1) == a == F(2) for a in alpha_set)
    invalid = 0
    for pattern in [(), ((F(2), F(1)),), ((math.inf, math.inf),), ((-math.inf, -math.inf),)]:
        try:
            bounds(pattern)
        except ValueError:
            invalid += 1
    assert invalid == 4
    report = {'status':'pass','arithmetic':'fractions.Fraction with infinity sentinels, no tolerance',
              'range_patterns_by_marker_count':counts,'total_range_patterns':sum(counts.values()),
              'two_cell_two_marker_pairs_exhaustive_on_declared_grid':pair_checks,
              'two_cell_three_marker_pairs_fixed_subset':n3_checks,
              'test_grid':[str(x) for x in grid], 'interval_inventory_size':len(intervals),
              'known_answer_examples':{'nested_intervals_range':['0','6'], 'one_marker_unbounded_range':['0','0'],
                  'missing_plus_disjoint_range':['2','infinity'], 'bounded_interaction_alpha':['1','2'],
                  'unrestricted_interaction_opposite_effects':True,'strict_reversal_rejected':True,'touching_ties_accepted':True},
              'invalid_input_rejections':invalid,
              'test_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':'Synthetic finite checks supplement the general proofs; not empirical data, exhaustive real-interval enumeration, or independent external review.'}
    Path(__file__).with_name('synthetic_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Optimized Python is not allowed for assertion-based checks.')
    main()
