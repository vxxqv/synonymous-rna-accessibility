"""Assumption and boundary checks for additive product threshold queries."""
from fractions import Fraction
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.product_closure import (component_closure, connection_witness,
    exact_product_bottleneck, product_connected)

P3 = ((1,), (0, 2), (1,))


def test_projection_containment_does_not_imply_component_membership():
    graphs = (P3, P3)
    scores = ((2, 0, 1),)*2
    closure = component_closure(graphs, scores, (0, 0), 2)
    assert closure.coordinates == (frozenset(range(3)),)*2
    assert not product_connected(graphs, scores, (0, 0), (2, 2), 2)
    assert connection_witness(graphs, scores, (0, 0), (2, 2), 2) is None
    assert exact_product_bottleneck(graphs, scores, (0, 0), (2, 2)) == 1


def test_rational_equality_threshold_and_identical_endpoint():
    scores = ((Fraction(1, 3), Fraction(-2, 7), Fraction(5, 11)),)
    assert product_connected((P3,), scores, (0,), (2,), Fraction(-2, 7))
    assert not product_connected((P3,), scores, (0,), (2,), Fraction(-2, 7)+Fraction(1, 1000))
    assert exact_product_bottleneck((P3,), scores, (0,), (0,)) == Fraction(1, 3)


def test_baseline_disconnection_and_singleton_factors():
    assert exact_product_bottleneck((((), ()),), ((0, 1),), (0,), (1,)) is None
    assert product_connected((((),),), ((3,),), (0,), (0,), 3)
    assert connection_witness((((),),), ((3,),), (0,), (0,), 3) == ()


def test_exact_input_and_undirected_assumptions_are_checked():
    with pytest.raises(TypeError, match='Exact integers'):
        component_closure((P3,), ((1., 0., 1.),), (0,), 0)
    with pytest.raises(ValueError, match='undirected'):
        component_closure((((1,), ()),), ((1, 1),), (0,), 0)
    with pytest.raises(ValueError, match='below threshold'):
        component_closure((P3,), ((1, 0, 1),), (1,), 1)


def test_hundred_factor_certificate_without_product_enumeration():
    dimension = 100
    graphs = (P3,)*dimension
    scores = ((2, 0, 1),)*dimension
    x = (0,)*dimension
    target = (2,)*50+(0,)*50
    moves = connection_witness(graphs, scores, x, target, dimension)
    assert moves is not None
    state = list(x)
    for coordinate, vertex in moves:
        assert vertex in graphs[coordinate][state[coordinate]]
        state[coordinate] = vertex
        assert sum(scores[i][v] for i, v in enumerate(state)) >= dimension
    assert tuple(state) == target
    assert len(moves) <= 2*sum((len(g)-1)**2 for g in graphs)
    # The all-right state has score d but every neighbor has score d-1,
    # providing an independent direct certificate of its isolation.
    isolated = (2,)*dimension
    assert sum(scores[i][v] for i, v in enumerate(isolated)) == dimension
    assert not product_connected(graphs, scores, x, isolated, dimension)
