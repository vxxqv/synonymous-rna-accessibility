"""Hand-computable counterexamples for exact connectivity query certificates."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.acquisition import (certificate, choose_query,
    clairvoyant_additional_cost, exact_minimax, minimum_unknown_separator)


DIAMOND = ((1, 2), (0, 3), (0, 3), (1, 2))


def test_missing_is_unknown_and_single_negative_does_not_disconnect():
    assert certificate(DIAMOND, (1, -1, -1, 1), 0, 3) == 'unresolved'
    assert certificate(DIAMOND, (1, 0, -1, 1), 0, 3) == 'unresolved'
    assert certificate(DIAMOND, (1, 0, 0, 1), 0, 3) == 'separated'
    assert certificate(DIAMOND, (1, 1, -1, 1), 0, 3) == 'connected'


def test_clairvoyant_certificate_is_not_worst_case_query_cost():
    state = (1, -1, -1, 1)
    assert clairvoyant_additional_cost(DIAMOND, state, (1, 0, 1, 1), 0, 3) == 1
    assert clairvoyant_additional_cost(DIAMOND, state, (1, 0, 0, 1), 0, 3) == 2
    assert exact_minimax(DIAMOND, 0, 3)(state) == 2


def test_known_present_vertices_are_protected_from_separator():
    graph = ((1,), (0, 2), (1, 3), (2,))
    assert minimum_unknown_separator(graph, (1, 1, -1, 1), 0, 3) == [2]
    assert minimum_unknown_separator(graph, (1, -1, -1, 1), 0, 3) == [1]


def test_adjacent_endpoints_are_already_certified():
    assert certificate(((1,), (0,)), (1, 1), 0, 1) == 'connected'
    assert clairvoyant_additional_cost(((1,), (0,)), (1, 1), (1, 1), 0, 1) == 0
    with pytest.raises(ValueError, match='unnecessary'):
        choose_query(((1,), (0,)), (1, 1), 0, 1, 'degree', 0, np.random.default_rng(1))


def test_invalid_observations_fail_closed():
    with pytest.raises(ValueError, match='known present'):
        certificate(DIAMOND, (-1, -1, -1, 1), 0, 3)
    with pytest.raises(ValueError, match='disagrees'):
        clairvoyant_additional_cost(DIAMOND, (1, 0, -1, 1), (1, 1, 1, 1), 0, 3)
