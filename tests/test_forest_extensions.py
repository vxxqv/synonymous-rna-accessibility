"""Independent counting identities and traversal-depth regression."""
from pathlib import Path
import math
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis/expansion'))
from graph_compiler import forest_extensions, subset_extensions


def alternating_permutations(n):
    row = [1]
    for m in range(1, n+1):
        previous = row
        row = [0]
        for k in range(1, m+1):
            row.append(row[-1]+previous[m-k])
    return row[-1]


def fence(n):
    return {(i, i+1) if i % 2 == 0 else (i+1, i) for i in range(n-1)}


def test_forest_counting_boundaries_and_independent_oracles():
    assert forest_extensions(0, set()) == 1
    assert forest_extensions(7, set()) == math.factorial(7)
    assert forest_extensions(3, {(0, 1), (1, 2), (0, 2)}) is None
    for n in range(1, 11):
        assert forest_extensions(n, fence(n)) == alternating_permutations(n)
        assert forest_extensions(n, fence(n)) == subset_extensions(n, fence(n))
    assert forest_extensions(6, {(0, 1), (0, 2), (3, 4), (3, 5)}) == 80


def test_valid_fence_does_not_use_recursive_tree_traversal():
    old_limit = sys.getrecursionlimit()
    try:
        sys.setrecursionlimit(100)
        assert forest_extensions(150, fence(150)) == alternating_permutations(150)
    finally:
        sys.setrecursionlimit(old_limit)
