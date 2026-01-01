"""Frozen synthetic query benchmark. Validation is disabled in this version.

Run from repository root:
    python analysis/programme/acquisition_benchmark.py --output PATH
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.acquisition import (POLICIES, certificate, choose_query,
    clairvoyant_additional_cost, exact_minimax, minimum_unknown_separator)
from src.barrier_identification import threshold_status


def adjacency_for(labels):
    return tuple(tuple(j for j, other in enumerate(labels) if i != j and
                       sum(a != b for a, b in zip(label, other)) == 1)
                 for i, label in enumerate(labels))


def graph_inventory():
    graphs = []
    for d in (3, 4, 6):
        labels = list(product((0, 1), repeat=d))
        graphs.append((f'Q{d}', labels, adjacency_for(labels), True))
    for d in (5, 7):
        labels = [bits for bits in product((0, 1), repeat=d)
                  if all(not (bits[i] and bits[i+1]) for i in range(d-1))]
        graphs.append((f'independent_path_{d}', labels, adjacency_for(labels), True))
    codons = ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG']
    graphs.append(('leucine', codons, adjacency_for(codons), False))
    labels = [a+b for a, b in product(codons, repeat=2)]
    graphs.append(('leucine_product', labels, adjacency_for(labels), False))
    return graphs


def nx_graph(adjacency):
    graph = nx.Graph()
    graph.add_nodes_from(range(len(adjacency)))
    graph.add_edges_from((u, v) for u, ns in enumerate(adjacency) for v in ns if u < v)
    return graph


def endpoint_pair(graph):
    distances = dict(nx.all_pairs_shortest_path_length(graph))
    pairs = [(distances[u][v], u, v) for u in sorted(graph)
             for v in sorted(graph) if u < v]
    diameter = max(d for d, _, _ in pairs)
    a, b = min((a, b) for d, a, b in pairs if d == diameter)
    return a, b, diameter


def truth_status(graph, truth, a, b):
    selected = graph.subgraph([v for v, state in enumerate(truth) if state == 1])
    return 'connected' if nx.has_path(selected, a, b) else 'separated'


def matrix_heat(adjacency):
    matrix = np.zeros((len(adjacency), len(adjacency)), dtype=float)
    for u, neighbors in enumerate(adjacency):
        matrix[u, u] = len(neighbors)
        for v in neighbors:
            matrix[u, v] = -1
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.exp(-values)) @ vectors.T


def field_values(family, labels, heat, rng):
    n = len(labels)
    if family.startswith('bernoulli_'):
        return (rng.random(n) < float(family.split('_')[1])).astype(int), None
    if family == 'heat':
        score = heat @ rng.standard_normal(n)
    else:
        coordinates = 2*np.asarray(labels, dtype=float)-1
        score = coordinates @ rng.standard_normal(coordinates.shape[1])
        if family == 'pairwise':
            for i, j in combinations(range(coordinates.shape[1]), 2):
                score += rng.standard_normal()*coordinates[:, i]*coordinates[:, j]
    threshold = float(np.median(score))
    return (score >= threshold).astype(int), threshold


def audit_prefix(adjacency, state, a, b, complete_truth_status):
    status = certificate(adjacency, state, a, b)
    lower = [float(s == 1) for s in state]
    upper = [float(s != 0) for s in state]
    # Orthogonal implementation: heap widest paths instead of Boolean BFS.
    alternative = threshold_status(adjacency, lower, upper, a, b, .5)
    if alternative != status:
        raise AssertionError(('Certificate kernels disagree', status, alternative))
    if status != 'unresolved' and status != complete_truth_status:
        raise AssertionError(('False certificate', status, complete_truth_status))
    return status


def exact_checks():
    """Enumerate every Q3 ternary observation state and every compatible truth."""
    _, labels, adjacency, _ = graph_inventory()[0]
    graph = nx_graph(adjacency)
    a, b, _ = endpoint_pair(graph)
    internal = [v for v in graph if v not in (a, b)]
    minimax = exact_minimax(adjacency, a, b)
    counts = Counter()
    for partial in product((-1, 0, 1), repeat=len(internal)):
        state = [1]*len(labels)
        for v, value in zip(internal, partial):
            state[v] = value
        status = certificate(adjacency, state, a, b)
        unknown = [v for v in internal if state[v] == -1]
        if status == 'unresolved':
            candidates = []
            for k in range(1, len(unknown)+1):
                for subset in combinations(unknown, k):
                    allowed = [v for v, value in enumerate(state)
                               if value != 0 and v not in subset]
                    if not nx.has_path(graph.subgraph(allowed), a, b):
                        candidates.append(subset)
                if candidates:
                    break
            expected = list(min(candidates))
            if minimum_unknown_separator(adjacency, state, a, b) != expected:
                raise AssertionError('Minimum unknown separator or tie rule is incorrect')
            counts['separator_comparisons'] += 1
        possible_truths = set()
        for completion in product((0, 1), repeat=len(unknown)):
            truth = state.copy()
            for v, value in zip(unknown, completion):
                truth[v] = value
            actual = truth_status(graph, truth, a, b)
            audit_prefix(adjacency, state, a, b, actual)
            possible_truths.add(actual)
            brute_cost = None
            for k in range(len(unknown)+1):
                found = False
                for subset in combinations(unknown, k):
                    trial = state.copy()
                    for v in subset:
                        trial[v] = truth[v]
                    if certificate(adjacency, trial, a, b) != 'unresolved':
                        brute_cost = k
                        found = True
                        break
                if found:
                    break
            if clairvoyant_additional_cost(adjacency, state, truth, a, b) != brute_cost:
                raise AssertionError('Clairvoyant path/cut bound disagrees with enumeration')
            counts['completion_comparisons'] += 1
        if (status == 'unresolved') != (len(possible_truths) == 2):
            raise AssertionError('Conditional certificate is not complete over all completions')
        counts['partial_states'] += 1
        counts['minimax_states_evaluated'] += int(minimax(tuple(state)) >= 0)
    initial = tuple(1 if v in (a, b) else -1 for v in graph)
    counts['q3_endpoint_only_minimax'] = minimax(initial)
    # Adjacent endpoints and already disconnected graphs are zero-query cases.
    for adj, state, pair, expected in [
            (((1,), (0,)), (1, 1), (0, 1), 'connected'),
            (((), ()), (1, 1), (0, 1), 'separated'),
            (((),), (1,), (0, 0), 'connected')]:
        if certificate(adj, state, *pair) != expected:
            raise AssertionError('Boundary case failed')
    counts['boundary_cases'] = 3
    return dict(counts)


def paired_summary(records):
    by_world = {}
    for row in records:
        by_world.setdefault(row['world_id'], {})[row['policy']] = row
    overall = {}
    for policy in POLICIES:
        rows = [world[policy] for world in by_world.values()]
        counts = [r['queries'] for r in rows]
        overall[policy] = {'median_queries': float(np.median(counts)),
            'mean_queries': float(np.mean(counts)), 'maximum_queries': max(counts),
            'sum_queries': sum(counts),
            'mean_excess_over_clairvoyant': float(np.mean([
                r['queries']-r['clairvoyant_additional'] for r in rows]))}
    paired = {}
    for baseline in POLICIES:
        if baseline == 'path_cut':
            continue
        differences = [w['path_cut']['queries']-w[baseline]['queries']
                       for w in by_world.values()]
        paired[baseline] = {'mean_difference': float(np.mean(differences)),
            'median_difference': float(np.median(differences)),
            'wins': sum(d < 0 for d in differences), 'ties': differences.count(0),
            'losses': sum(d > 0 for d in differences),
            'largest_extra_queries': max(differences),
            'largest_saved_queries': -min(differences)}
    strata = []
    keys = sorted(set((r['graph'], r['family'], r['initial_condition']) for r in records))
    for graph, family, initial in keys:
        selected = [r for r in records if
                    (r['graph'], r['family'], r['initial_condition']) ==
                    (graph, family, initial)]
        strata.append({'graph': graph, 'family': family, 'initial_condition': initial,
                      'policies': {p: {'mean': float(np.mean([r['queries'] for r in selected if r['policy'] == p])),
                                       'median': float(np.median([r['queries'] for r in selected if r['policy'] == p])),
                                       'maximum': max(r['queries'] for r in selected if r['policy'] == p)}
                                   for p in POLICIES}})
    success = all(overall['path_cut']['median_queries'] <= .9*overall[p]['median_queries']
                  and overall['path_cut']['maximum_queries'] <= overall[p]['maximum_queries']
                  for p in ('degree', 'lazy_forward', 'lazy_bisection'))
    return {'overall': overall, 'paired_path_cut_minus_baseline': paired,
            'graph_field_initial_strata': strata,
            'development_only_efficiency_criterion_met': success,
            'interpretation': 'Development simulations only. This is not validation or biological evidence.'}


def run(output):
    start = perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    checks = exact_checks()
    records = []
    worlds = []
    inventory = []
    prefixes = 0
    for gi, (name, labels, adjacency, binary) in enumerate(graph_inventory()):
        graph = nx_graph(adjacency)
        a, b, diameter = endpoint_pair(graph)
        internal = [v for v in graph if v not in (a, b)]
        heat = matrix_heat(adjacency)
        families = ['bernoulli_0.25', 'bernoulli_0.5', 'bernoulli_0.75', 'heat']
        if binary:
            families += ['additive', 'pairwise']
        inventory.append({'name': name, 'vertices': len(labels), 'edges': graph.number_of_edges(),
            'labels': labels, 'adjacency': adjacency, 'endpoints': [a, b],
            'diameter': diameter, 'families': families})
        minimax = exact_minimax(adjacency, a, b) if name == 'Q3' else None
        for fi, family in enumerate(families):
            for seed in range(1000, 1010):
                rng = np.random.default_rng(np.random.SeedSequence([seed, gi, fi, 0]))
                truth, threshold = field_values(family, labels, heat, rng)
                truth[a] = truth[b] = 1
                truth = truth.tolist()
                actual = truth_status(graph, truth, a, b)
                initial_rng = np.random.default_rng(np.random.SeedSequence([seed, gi, 0, 1]))
                initial_extra = sorted(initial_rng.permutation(internal)[:len(internal)//4].tolist())
                for initial_condition, extra in [('endpoints', []), ('quarter', initial_extra)]:
                    world_id = f'{name}:{family}:{seed}:{initial_condition}'
                    state0 = [truth[v] if v in (a, b) or v in extra else -1 for v in graph]
                    lower_bound = clairvoyant_additional_cost(adjacency, state0, truth, a, b)
                    minimax_cost = minimax(tuple(state0)) if minimax else None
                    worlds.append({'world_id': world_id, 'graph': name, 'family': family,
                        'seed': seed, 'initial_condition': initial_condition,
                        'initial_extra': extra, 'truth': truth, 'threshold': threshold,
                        'truth_status': actual, 'clairvoyant_additional': lower_bound,
                        'q3_minimax_remaining': minimax_cost})
                    for policy in POLICIES:
                        policy_start = perf_counter()
                        state = state0.copy()
                        policy_rng = np.random.default_rng(np.random.SeedSequence([seed, gi, fi, 2]))
                        queries = []
                        while True:
                            status = audit_prefix(adjacency, state, a, b, actual)
                            prefixes += 1
                            if status != 'unresolved':
                                break
                            v = choose_query(adjacency, state, a, b, policy, len(queries), policy_rng)
                            if v in (a, b) or state[v] != -1:
                                raise AssertionError('Policy queried an observed vertex')
                            state[v] = truth[v]
                            queries.append(v)
                        if not lower_bound <= len(queries) <= len(internal)-len(extra):
                            raise AssertionError('Query cost violates certificate bounds')
                        records.append({'world_id': world_id, 'graph': name, 'family': family,
                            'seed': seed, 'initial_condition': initial_condition, 'policy': policy,
                            'queries': len(queries), 'query_order': queries, 'certificate': status,
                            'clairvoyant_additional': lower_bound, 'q3_minimax_remaining': minimax_cost,
                            'policy_seconds': perf_counter()-policy_start})
        print(json.dumps({'graph_completed': name, 'worlds': len(worlds), 'policy_runs': len(records)}), flush=True)
    if len(worlds) != 760 or len(records) != 3800:
        raise AssertionError('Coverage differs from frozen contract')
    summary = paired_summary(records)
    summary.update({'phase': 'development', 'validation_outcomes_examined': False,
        'worlds': len(worlds), 'policy_runs': len(records), 'checked_query_prefixes': prefixes,
        'false_certificates': 0, 'exact_checks': checks,
        'connected_worlds': sum(w['truth_status'] == 'connected' for w in worlds),
        'separated_worlds': sum(w['truth_status'] == 'separated' for w in worlds),
        'runtime_seconds': perf_counter()-start,
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'networkx': nx.__version__},
        'input_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/acquisition.py', 'src/barrier_identification.py',
            'analysis/programme/acquisition_benchmark.py',
            'analysis/programme/acquisition_protocol_v1.json',
            'analysis/programme/acquisition_conventions_v1.json')}})
    for name, data in [('summary.json', summary), ('worlds.json', worlds),
                       ('policy_results.json', records), ('graphs.json', inventory)]:
        (output/name).write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ('worlds', 'policy_runs', 'checked_query_prefixes',
        'false_certificates', 'exact_checks', 'overall', 'development_only_efficiency_criterion_met')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.output)
