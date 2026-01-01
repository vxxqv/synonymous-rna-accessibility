"""Bounded exact tests for overlapping standard-code constraints.

The compiler transfer is a logical derivation. Enumerations test its implementation
and probe reachability; a negative search is not a general geodesicity proof.
"""
from pathlib import Path
from itertools import product, combinations
from collections import defaultdict, deque
import argparse
import json
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'analysis/expansion'))
from feasibility_2sat import CODON_CNF, BITS, BASE, solve_cnf, evaluate, encode
from Bio.Data import CodonTable

COMPLEMENT = str.maketrans('ACGU', 'UGCA')
FORWARD = {c.replace('T', 'U'): a for c, a in CodonTable.unambiguous_dna_by_id[1].forward_table.items()}
PAIR = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}


def occurrences(n, strand, offset):
    positions = list(range(offset, n)) if strand == '+' else list(range(n-1-offset, -1, -1))
    return [(tuple(positions[k:k+3]), strand == '-') for k in range(0, len(positions)-2, 3)]


def read_codon(seq, positions, complement):
    s = ''.join(seq[p] for p in positions)
    return s.translate(COMPLEMENT) if complement else s


def labels(seq, codons):
    c = [read_codon(seq, positions, complement) for positions, complement in codons]
    if any(x not in FORWARD for x in c):
        return None
    return tuple(FORWARD[x] for x in c)


def signature(seq, codons, peptide, edges):
    serine = tuple(int(read_codon(seq, positions, complement).startswith('AG'))
                   for (positions, complement), aa in zip(codons, peptide) if aa == 'S')
    return serine + tuple(int(seq[min(i, j)] in 'AG') for i, j in edges)


def compile_overlap(n, codons, peptide, edges=(), fixed=None):
    if len(codons) != len(peptide):
        raise ValueError('One standard sense amino acid per codon occurrence')
    clauses = []
    for (positions, complement), aa in zip(codons, peptide):
        if aa not in CODON_CNF or len(positions) != 3 or len(set(positions)) != 3 or any(not 0 <= p < n for p in positions):
            raise ValueError('Invalid codon occurrence')
        for clause in CODON_CNF[aa]:
            mapped = []
            for literal in clause:
                local = abs(literal) - 1
                variable = 2*positions[local//2] + local % 2 + 1
                # Complement sends (r,d) to (1-r,1-d); reverse order is in positions.
                sign = (1 if literal > 0 else -1) * (-1 if complement else 1)
                mapped.append(sign*variable)
            clauses.append(tuple(mapped))
    for i, j in edges:
        if i == j or not (0 <= i < n and 0 <= j < n):
            raise ValueError('Invalid pairing obligation')
        ri, rj, di, dj = 2*i+1, 2*j+1, 2*i+2, 2*j+2
        clauses.extend(((ri, rj), (-ri, -rj), (-di, -dj)))
    for p, base in (fixed or {}).items():
        if not 0 <= p < n or base not in BITS:
            raise ValueError('Invalid pin')
        clauses.extend(((2*p+k+1) if bit else -(2*p+k+1),)
                       for k, bit in enumerate(BITS[base]))
    assert len(clauses) <= 66*len(codons) + 3*len(edges) + 2*len(fixed or {})
    return 2*n, tuple(clauses)


def direct(seq, codons, peptide, edges, fixed):
    return (labels(seq, codons) == tuple(peptide)
            and all(seq[i]+seq[j] in PAIR for i, j in edges)
            and all(seq[i] == b for i, b in fixed.items()))


def compiler_checks():
    rng = random.Random(73691)
    observations = []
    assignments = 0
    for k in range(18):
        n = 4 if k < 6 else 6
        strand, offset = [('+', 1), ('+', 2), ('-', 0), ('-', 1), ('-', 2), ('+', 0)][k % 6]
        codons = occurrences(n, '+', 0) + occurrences(n, strand, offset)
        peptide = tuple(rng.choice(sorted(CODON_CNF)) for _ in codons)
        edges = tuple(rng.sample(list(combinations(range(n), 2)), k % 3))
        fixed = {rng.randrange(n): rng.choice('ACGU')} if k % 2 else {}
        nvars, clauses = compile_overlap(n, codons, peptide, edges, fixed)
        solutions = []
        for x in product('ACGU', repeat=n):
            seq = ''.join(x)
            truth = direct(seq, codons, peptide, edges, fixed)
            assert truth == evaluate(clauses, encode(seq))
            assignments += 1
            if truth:
                solutions.append(seq)
        answer = solve_cnf(nvars, clauses)
        assert (answer is not None) == bool(solutions)
        if answer is not None:
            witness = ''.join(BASE[answer[2*p:2*p+2]] for p in range(n))
            assert direct(witness, codons, peptide, edges, fixed)
        observations.append({'n': n, 'second_frame': [strand, offset], 'peptide': peptide,
                             'pair_obligations': edges, 'pins': fixed, 'solutions': len(solutions)})
    # SAT, empty, contradiction and reverse-complement boundary examples.
    for n, codons, peptide, edges, fixed in [
        (0, [], [], [], {}),
        (3, [((0,1,2),False)], ['K'], [], {0:'C'}),
        (3, [((2,1,0),True)], ['K'], [], {}),
        (1, [], [], [], {0:'A'})]:
        nvars, clauses = compile_overlap(n, codons, peptide, edges, fixed)
        truth = [x for x in map(''.join, product('ACGU', repeat=n)) if direct(x,codons,peptide,edges,fixed)]
        assert (solve_cnf(nvars,clauses) is not None) == bool(truth)
    return {'assignment_truth_checks': assignments, 'random_instances': observations,
            'boundary_instances': 4}


def probe_case(n, second, pair_edge):
    codons = occurrences(n, '+', 0) + occurrences(n, *second)
    groups = defaultdict(list)
    for x in product('ACGU', repeat=n):
        seq = ''.join(x)
        peptide = labels(seq, codons)
        if peptide is not None and all(seq[i]+seq[j] in PAIR for i,j in pair_edge):
            groups[peptide].append(seq)
    components = pairs = 0
    extra = []
    detours = []
    for peptide, sequences in groups.items():
        lookup = {s:i for i,s in enumerate(sequences)}
        graph = [[] for _ in sequences]
        for u, seq in enumerate(sequences):
            for p in range(n):
                for base in 'ACGU':
                    if base != seq[p]:
                        v = lookup.get(seq[:p]+base+seq[p+1:])
                        if v is not None:
                            graph[u].append(v)
        signatures = defaultdict(list)
        unvisited = set(range(len(sequences)))
        while unvisited:
            root = min(unvisited)
            stack, component = [root], {root}
            unvisited.remove(root)
            while stack:
                for v in graph[stack.pop()]:
                    if v in unvisited:
                        unvisited.remove(v)
                        component.add(v)
                        stack.append(v)
            sig = signature(sequences[root], codons, peptide, pair_edge)
            assert all(signature(sequences[u], codons, peptide, pair_edge) == sig for u in component)
            signatures[sig].append(sequences[root])
            components += 1
            for start in sorted(component):
                distance, queue = {start: 0}, deque([start])
                while queue:
                    u = queue.popleft()
                    for v in graph[u]:
                        if v not in distance:
                            distance[v] = distance[u] + 1
                            queue.append(v)
                for end in sorted(component):
                    if end <= start:
                        continue
                    h = sum(a != b for a,b in zip(sequences[start], sequences[end]))
                    pairs += 1
                    if distance[end] != h and len(detours) < 5:
                        detours.append({'peptide': peptide, 'a': sequences[start], 'b': sequences[end],
                                        'Hamming': h, 'distance': distance[end]})
        for sig, reps in signatures.items():
            if len(reps) > 1 and len(extra) < 5:
                extra.append({'peptide': peptide, 'signature': sig, 'component_representatives': reps})
    return {'n': n, 'second_frame': second, 'pair_edges': pair_edge,
            'nonempty_fixed_peptide_sets': len(groups), 'sequences': sum(map(len,groups.values())),
            'components': components, 'connected_unordered_endpoint_pairs': pairs,
            'new_signature_obstructions': extra, 'non_Hamming_detours': detours}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    checks = compiler_checks()
    # All five relative frame choices at lengths 6 and 7, first without pairs;
    # then every single pairing edge at length 7. No favorable subset selection.
    rows = []
    for n in (6,7):
        for second in [('+',1),('+',2),('-',0),('-',1),('-',2)]:
            for edges in ([()] if n == 6 else [()] + [(e,) for e in combinations(range(n),2)]):
                rows.append(probe_case(n, second, edges))
    result = {'scope': 'Bounded exact implementation and counterexample search, not a proof of general reachability',
              'compiler_checks': checks, 'configurations': rows,
              'totals': {'configurations': len(rows),
                         'connected_unordered_endpoint_pairs': sum(x['connected_unordered_endpoint_pairs'] for x in rows),
                         'configurations_with_new_obstruction': sum(bool(x['new_signature_obstructions']) for x in rows),
                         'configurations_with_detour': sum(bool(x['non_Hamming_detours']) for x in rows)}}
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result['totals']))
