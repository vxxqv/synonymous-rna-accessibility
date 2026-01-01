"""Semantic regressions for exact overlapping-code quotient and compact paths."""
from itertools import product
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.programme.overlap_v2_quotient import (
    accepts, compile_sector, endpoint_mutation_order, endpoint_path,
    frozen_variables, implication_components,
)


def test_nonempty_domains_can_have_a_contradictory_scc():
    occurrences = (((0, 1, 4), True), ((2, 3, 4), False))
    compiled = compile_sector(5, occurrences, 'LL', ((0, 2),), {0: 'A'},
                              ((), (True, False)))
    assert all(compiled['domains'])
    assert set(compiled['clauses']) == {(-2,), (-1, 2), (1,)}
    assert compiled['unsat'] is True
    assert not any(accepts(compiled, bits) for bits in product((0, 1), repeat=2))
    component, frozen = implication_components(compiled)
    assert component[0] == component[1]
    assert frozen == frozen_variables(compiled) == compiled['frozen_variables']


def test_unconstrained_within_fibre_substitutions_are_retained():
    # A and G share the purine bucket; C and U share the pyrimidine bucket.
    assert endpoint_mutation_order(2, (), '', (), {}, 'AC', 'GU') == (0, 1)
    assert endpoint_path(2, (), '', (), {}, 'AC', 'GU') == ['AC', 'GC', 'GU']


def test_bucket_changes_and_within_bucket_changes_each_occur_once():
    order = endpoint_mutation_order(3, (), '', (), {}, 'ACC', 'GAA')
    assert len(order) == len(set(order)) == 3
    assert set(order) == {0, 1, 2}
    path = endpoint_path(3, (), '', (), {}, 'ACC', 'GAA')
    assert path[0] == 'ACC' and path[-1] == 'GAA' and len(path) == 4
    assert all(sum(a != b for a, b in zip(x, y)) == 1 for x, y in zip(path, path[1:]))


def test_abstract_arg_reversal_requires_scc_invariant():
    occurrences = (((0, 1, 2), False), ((2, 1, 0), False))
    compiled = compile_sector(3, occurrences, 'RR', (), {}, ((), ()))
    assert not compiled['unsat']
    assert frozen_variables(compiled) == (0, 1)
    assert endpoint_mutation_order(3, occurrences, 'RR', (), {}, 'CGC', 'AGA') is None
    assert endpoint_path(3, occurrences, 'RR', (), {}, 'CGC', 'AGA') is None


def test_long_precedence_chain_returns_compact_reverse_order():
    k = 2000
    occurrences = [((i, k, i+1), False) for i in range(k-1)]
    occurrences.append(((k-1, k, k+1), False))
    order = endpoint_mutation_order(k+2, occurrences, 'R'*k, (), {k: 'G', k+1: 'A'},
                                    'C'*k+'GA', 'A'*k+'GA')
    # Each implication A_i => A_(i+1) makes this the unique valid order.
    assert order == tuple(reversed(range(k)))


@pytest.mark.parametrize('x,y', [('X', 'A'), ('A', 'T'), ('a', 'G')])
def test_invalid_alphabet_rejected_even_without_obligations(x, y):
    with pytest.raises(ValueError, match='valid endpoints'):
        endpoint_mutation_order(1, (), '', (), {}, x, y)
    with pytest.raises(ValueError, match='valid endpoints'):
        endpoint_path(1, (), '', (), {}, x, y)


def test_empty_instance_and_identical_endpoints():
    compiled = compile_sector(0, (), '', (), {}, ((), ()))
    assert not compiled['unsat'] and accepts(compiled, ())
    assert frozen_variables(compiled) == ()
    assert endpoint_mutation_order(0, (), '', (), {}, '', '') == ()
    assert endpoint_path(0, (), '', (), {}, '', '') == ['']
    assert endpoint_mutation_order(1, (), '', (), {0: 'G'}, 'G', 'G') == ()


@pytest.mark.parametrize('n,edges,pins,x,y', [
    (-1, (), {}, '', ''),
    (1, (), {1: 'A'}, 'A', 'A'),
    (1, ((0, 1),), {}, 'A', 'A'),
    (1, ((0, 0),), {}, 'A', 'A'),
    (1, (), {}, '', 'A'),
])
def test_no_obligations_does_not_bypass_input_bounds(n, edges, pins, x, y):
    with pytest.raises(ValueError):
        endpoint_mutation_order(n, (), '', edges, pins, x, y)
