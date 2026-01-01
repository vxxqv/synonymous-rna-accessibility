"""Independent exact certificate and threshold-component audit. No production imports.

Reads only an explicitly supplied saved result and sequence source contract.
The saved counts are used to check ratio arithmetic, not to reopen count files.
"""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python is not allowed')
from collections import deque
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
from statistics import median


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(v):
    if v in ('inf', '+inf'):
        return math.inf
    if v == '-inf':
        return -math.inf
    return Q(v)


def encode(v):
    return str(v)


def design(states):
    levels = [sorted({row[j] for row in states}) for j in range(len(states[0]))]
    return [[1, *[int(row[j] == a) for j, values in enumerate(levels)
                  for a in values[1:]]] for row in states]


def certificate(cert, X, L, U):
    n, p = len(X), len(X[0])
    forced = {(i, j) for i in range(n) for j in range(n)
              if U[i] is not None and U[i] < L[j]}
    require(cert['forced_orders'] == len(forced), 'Forced order count')
    require(cert['exact_checked'] is True, 'Missing exact_checked flag')
    if cert['status'] == 'infeasible':
        weights = []
        cancellation = [Q(0)] * p
        support = []
        for item in cert['dual']:
            i, j = item['lower_vertex'], item['higher_vertex']
            require((i, j) in forced, 'Dual order not strictly forced by intervals')
            w = Q(item['weight'])
            require(w > 0, 'Nonpositive supported dual weight')
            a = [X[j][k]-X[i][k] for k in range(p)]
            require(a == item['difference'], 'Incorrect stored feature difference')
            cancellation = [v+w*z for v, z in zip(cancellation, a)]
            weights.append(w)
            support.append({'lower_vertex': i, 'higher_vertex': j,
                            'weight': str(w), 'strict_envelope_gap': str(L[j]-U[i]),
                            'difference': a})
        require(sum(weights) == 1 and all(v == 0 for v in cancellation),
                'Dual does not cancel exactly to zero with unit total weight')
        return {'status': 'infeasible', 'forced_orders': len(forced),
                'support_size': len(support), 'support': support,
                'weighted_difference_sum': list(map(str, cancellation))}
    require(cert['status'] == 'feasible', 'Unknown certificate status')
    b = list(map(Q, cert['beta']))
    require(len(b) == p, 'Beta dimension')
    t = [sum(a*z for a, z in zip(row, b)) for row in X]
    y = list(map(Q, cert['response']))
    require(list(map(Q, cert['latent'])) == t and len(y) == n, 'Latent/response values')
    require(all(t[j]-t[i] >= 1 for i, j in forced), 'Primal unit margins')
    require(all(L[i] <= y[i] and (U[i] is None or y[i] <= U[i]) for i in range(n)),
            'Response outside intervals')
    require(all(y[i] <= y[j] for i in range(n) for j in range(n) if t[i] <= t[j]),
            'Stored response not a common monotone function, including ties')
    return {'status': 'feasible', 'forced_orders': len(forced)}


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        require(0 <= a < n and 0 <= b < n and a != b, 'Invalid graph edge')
        adj[a].add(b)
        adj[b].add(a)
    return adj


def component(adj, scores, threshold, start, strict=False):
    good = (lambda v: scores[v] > threshold) if strict else (lambda v: scores[v] >= threshold)
    if not good(start):
        return set()
    reached = {start}
    queue = deque([start])
    while queue:
        for j in adj[queue.popleft()]:
            if j not in reached and good(j):
                reached.add(j)
                queue.append(j)
    return reached


def capacities(scores, adj):
    """Independent algorithm: descending threshold components, no max-min closure."""
    n = len(scores)
    result = [[-math.inf]*n for _ in range(n)]
    for threshold in sorted(set(scores), reverse=True):
        unseen = {v for v in range(n) if scores[v] >= threshold}
        while unseen:
            vertices = component(adj, scores, threshold, min(unseen))
            unseen -= vertices
            for a in vertices:
                for b in vertices:
                    if result[a][b] == -math.inf:
                        result[a][b] = threshold
    return result


def witness(a, b, scores, adj, capacity):
    if capacity == -math.inf:
        return {'capacity': '-inf', 'path': None,
                'base_component_a': sorted(component(adj, [0]*len(scores), 0, a))}
    parent = {a: None}
    queue = deque([a])
    while queue and b not in parent:
        i = queue.popleft()
        for j in sorted(adj[i]):
            if scores[j] >= capacity and j not in parent:
                parent[j] = i
                queue.append(j)
    require(b in parent, 'Capacity lacks a path')
    path, v = [], b
    while v is not None:
        path.append(v)
        v = parent[v]
    path.reverse()
    require(min(scores[v] for v in path) == capacity, 'Witness bottleneck')
    cut = component(adj, scores, capacity, a, strict=True)
    require(b not in cut, 'Strict-superlevel cut fails')
    return {'capacity': str(capacity), 'path': path,
            'strict_superlevel_component_a': sorted(cut),
            'endpoint_blocks_strict_level': scores[a] <= capacity or scores[b] <= capacity}


def pair_check(records, L, U, adj):
    low = capacities(L, adj)
    high_scores = [math.inf if v is None else v for v in U]
    high = capacities(high_scores, adj)
    pairs = list(itertools.combinations(range(len(L)), 2))
    require(len(records) == len(pairs), 'Pair coverage')
    positive, sample = 0, None
    for record, (a, b) in zip(records, pairs):
        require((record['a'], record['b']) == (a, b), 'Pair order or duplicate')
        endpoint = min(L[a], L[b])
        width = max(Q(0), endpoint-high[a][b])
        nonempty = high[a][b] < endpoint
        for key, expected in [('endpoint_lower', endpoint), ('bottleneck_lower', low[a][b]),
                              ('bottleneck_upper', high[a][b]), ('certified_interval_width', width)]:
            require(number(record[key]) == expected, 'Pair arithmetic mismatch '+key)
        require(record['certified_nonempty'] is nonempty, 'Strict nonempty criterion')
        require(record['threshold_convention'] == '(bottleneck_upper, endpoint_lower] if nonempty; acceptance >= threshold',
                'Weak threshold convention changed')
        positive += nonempty
        if sample is None or (nonempty and not sample['certified_nonempty']):
            sample = {'a': a, 'b': b, 'certified_nonempty': nonempty,
                      'endpoint_lower': str(endpoint), 'width': str(width),
                      'lower_witness': witness(a, b, L, adj, low[a][b]),
                      'upper_witness': witness(a, b, high_scores, adj, high[a][b])}
    return {'pairs_checked': len(pairs), 'certified_nonempty': positive, 'example': sample}


def run(result, contract):
    rows = result['vertices_data']
    sequences = contract['sequences']
    require(len(rows) == len(sequences) == 72 and len(set(sequences)) == 72, 'Fixed 72 domain')
    require([r['sequence'] for r in rows] == sequences, 'Sequence domain/order')
    require([r['vertex'] for r in rows] == list(range(72)), 'Vertex indices')
    states = [[s[k:k+3] for k in (23, 38, 29, 32, 41)] for s in sequences]
    require([r['codon_states'] for r in rows] == states, 'Codon state annotation')
    X = design(states)
    variable_sites = [j for j in range(63) if len({s[j] for s in sequences}) > 1]
    Xn = design([[s[j] for j in variable_sites] for s in sequences])
    edges = [(i, j) for i, j in itertools.combinations(range(72), 2)
             if sum(x != y for x, y in zip(sequences[i], sequences[j])) == 1]
    require(len(edges) == result['edges'] == 204 and result['vertices'] == 72, 'Graph size')
    require(len(X[0]) == result['codon_design_columns'] == 8 and
            len(Xn[0]) == result['nucleotide_design_columns'] == 8, 'Feature size')
    adj = adjacency(72, edges)
    wt = 'GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
    w = sequences.index(wt)
    require(all(len(r['output_counts']) == 9 and len(r['input_counts']) == 3 for r in rows), 'Count vector sizes')
    for r in rows:
        require(all(v is None or type(v) is int and v >= 0 for v in r['output_counts']+r['input_counts']), 'Saved count integrity')
    fields = []
    for br, tr in itertools.product(range(9), range(3)):
        field = []
        for row in rows:
            a, b, c, d = row['output_counts'][br], rows[w]['output_counts'][br], row['input_counts'][tr], rows[w]['input_counts'][tr]
            v = None if any(x is None for x in (a,b,c,d)) or min(b,c,d) <= 0 else Q(a*d, b*c)
            field.append(v)
        fields.append(field)
    totals = result['normalization']['TR1_aligned_totals']
    require(len(totals) == 12, 'Normalization totals shape')
    def med(row):
        values = [Q(c, totals[j]) for j,c in enumerate(row['input_counts']) if c is not None and totals[j] > 0]
        return median(values) if values else None
    medians = list(map(med, rows))
    points = []
    for i, row in enumerate(rows):
        require([None if v is None else Q(v) for v in row['scenario_ratios']] == [f[i] for f in fields], 'Saved ratio arithmetic')
        source = []
        for br in range(9):
            a, b = row['output_counts'][br], rows[w]['output_counts'][br]
            v = None if a is None or b is None or b <= 0 or totals[br+3] <= 0 or medians[i] is None or medians[w] is None or medians[i] <= 0 or medians[w] <= 0 else Q(a,b)*medians[w]/medians[i]
            source.append(v)
        require([None if v is None else Q(v) for v in row['author_relative_by_output']] == source, 'Author relative arithmetic')
        point = None if None in source else sum(source)/9
        require((None if row['point'] is None else Q(row['point'])) == point, 'Point arithmetic')
        points.append(point)
    L, U = [], []
    for i, row in enumerate(rows):
        values = [f[i] for f in fields]
        l, u = (Q(0), None) if None in values else (min(values), max(values))
        require(Q(row['lower']) == l and (None if row['upper'] is None else Q(row['upper'])) == u, 'Envelope arithmetic')
        L.append(l); U.append(u)
    require(result['finite_envelopes'] == sum(v is not None for v in U) and
            result['unknown_envelopes'] == sum(v is None for v in U), 'Envelope coverage')
    primary = certificate(result['primary_codon_null'], X, L, U)
    nucleotide = certificate(result['nucleotide_encoding_equivalence_check'], Xn, L, U)
    primary_graph = pair_check(result['bounded_pairs'], L, U, adj)
    require(primary_graph['certified_nonempty'] == result['conservative_barrier_pairs'], 'Primary graph summary')
    require(len(result['scenario_results']) == 28, 'Scenario coverage')
    summaries = []
    for k, field in enumerate(fields+[points]):
        label = f'BR{k//3+1}_TR{k%3+1}' if k < 27 else 'mean_author_per_output_relative_ES'
        record = result['scenario_results'][k]
        require(record['scenario'] == label and record['unknown_vertices'] == field.count(None), 'Scenario identity/unknown coverage')
        lower = [Q(0) if v is None else v for v in field]
        cert = certificate(record['codon_null'], X, lower, field)
        graph = pair_check(result['all_scenario_pairs'][label], lower, field, adj)
        require(record['bounded_pair_count'] == graph['certified_nonempty'], 'Scenario graph summary')
        summaries.append({'scenario': label, 'certificate': cert, 'graph': graph})
    require(set(result['all_scenario_pairs']) == {x['scenario'] for x in summaries}, 'Extra/missing scenario fields')
    return {'status': 'PASS', 'vertices':72, 'edges':204, 'ratio_values_checked':72*27,
            'author_relative_values_checked':72*9, 'primary_certificate':primary,
            'nucleotide_certificate':nucleotide, 'primary_graph':primary_graph,
            'scenario_summaries':summaries, 'certificate_count':30,
            'pair_records_checked':29*2556,
            'limitations':['Certifies saved finite-domain arithmetic and rational implications, not source extraction or biological sampling coverage.',
                           'Shared count normalizations are dependent; empirical envelopes are not confidence intervals.',
                           'Graph obstruction does not establish causal interactions; absence of certified obstruction is not a universal connection claim.',
                           'TR1-aligned totals are treated as supplied source metadata; full 3072 source count files are not reopened.']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--result',required=True,type=Path)
    parser.add_argument('--result-sha256',required=True)
    parser.add_argument('--contract',required=True,type=Path)
    parser.add_argument('--contract-sha256',required=True)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    require(not args.out.exists(), 'Refusing existing output')
    require(sha(args.result) == args.result_sha256, 'Result hash mismatch')
    require(sha(args.contract) == args.contract_sha256, 'Contract hash mismatch')
    audit=run(json.loads(args.result.read_text()),json.loads(args.contract.read_text()))
    audit['input_hashes']={'result':sha(args.result),'contract':sha(args.contract),'verifier':sha(__file__)}
    args.out.mkdir(parents=True)
    (args.out/'verification.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:audit[k] for k in ('status','certificate_count','pair_records_checked')}))


if __name__ == '__main__':
    main()
