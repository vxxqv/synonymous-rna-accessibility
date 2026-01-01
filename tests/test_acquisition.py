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


def test_canonical_bisection_differs_from_path_midpoint_after_observations():
    path = tuple(tuple(v for v in (u-1, u+1) if 0 <= v < 9) for u in range(9))
    state = (1, -1, -1, -1, 1, -1, -1, -1, 1)
    rng = np.random.default_rng(1)
    assert choose_query(path, state, 0, 8, 'lazy_bisection', 1, rng) == 3
    assert choose_query(path, state, 0, 8, 'lazy_vertex_bisection', 1, rng) == 2
    assert choose_query(path, state, 0, 8, 'path_cut', 0, rng) == 3


def test_existing_output_stops_before_validation(tmp_path):
    from analysis.programme.acquisition_benchmark import run
    with pytest.raises(FileExistsError, match='new output directory'):
        run(tmp_path, 'validation')
    assert not list(tmp_path.iterdir())


def test_validation_hash_mismatch_fails_before_outcomes(tmp_path, monkeypatch):
    import json
    from analysis.programme import acquisition_benchmark as benchmark
    (tmp_path/'input.txt').write_text('changed source', encoding='utf-8')
    freeze = tmp_path/'freeze.json'
    freeze.write_text(json.dumps({'input_sha256': {'input.txt': '0'*64},
        'validation_seeds': list(range(2000, 2050)),
        'expected_worlds': 3800, 'expected_policy_runs': 22800}), encoding='utf-8')
    monkeypatch.setattr(benchmark, 'ROOT', tmp_path)
    monkeypatch.setattr(benchmark, 'VALIDATION_FREEZE', freeze)
    monkeypatch.setattr(benchmark, 'EXECUTION_INPUTS', ('input.txt',))
    with pytest.raises(ValueError, match='hash mismatch'):
        benchmark.verify_validation_freeze()


def test_uncertainty_uses_seed_blocks_not_individual_worlds():
    from analysis.programme.acquisition_benchmark import paired_summary
    from src.acquisition import POLICIES
    records = []
    for offset, seed in enumerate(range(2000, 2050)):
        for world in range(76):
            for policy in POLICIES:
                records.append({'world_id': f'{seed}:{world}', 'policy': policy,
                    'queries': offset if policy == 'path_cut' else 0,
                    'clairvoyant_additional': 0, 'seed': seed, 'graph': 'test',
                    'family': 'test', 'initial_condition': 'test'})
    summary = paired_summary(records, 'validation')
    interval = summary['paired_path_cut_minus_baseline']['random']['mean_difference_uncertainty']
    assert interval['standard_error'] == pytest.approx(np.sqrt(4.25))
    assert interval['degrees_of_freedom'] == 49
    assert np.mean(interval['interval']) == pytest.approx(24.5)
