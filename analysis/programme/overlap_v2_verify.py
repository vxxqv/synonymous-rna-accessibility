"""Bounded independent direct-translation and graph audit of overlap_v2.

Only abstract standard-code constraints and synthetic short strings are used.
No source transcript, organism or intervention data are read.
"""
import argparse
from collections import Counter, defaultdict, deque
import datetime
import hashlib
import itertools
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

from Bio.Data import CodonTable
from overlap_v2_quotient import compile_sector, accepts, project, frozen_variables, endpoint_path

ROOT = Path(__file__).resolve().parents[2]
FORWARD = {c.replace('T', 'U'): a for c, a in CodonTable.unambiguous_dna_by_id[1].forward_table.items()}
PAIR = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
COMPLEMENT = dict(zip('ACGU', 'UGCA'))
COUNTS = Counter()
ROWS = []


def read_direct(sequence, occurrence):
    positions, complement = occurrence
    bases = [sequence[p] for p in positions]
    return ''.join(COMPLEMENT[b] if complement else b for b in bases)


def peptide_direct(sequence, occurrences):
    codons = [read_direct(sequence, occurrence) for occurrence in occurrences]
    return tuple(FORWARD[c] for c in codons) if all(c in FORWARD for c in codons) else None


def signature_direct(sequence, occurrences, peptide, edges):
    return (tuple(read_direct(sequence, occurrence)[:2] == 'AG'
                  for occurrence, aa in zip(occurrences, peptide) if aa == 'S'),
            tuple(sequence[p] in 'AG' for p in sorted({p for edge in edges for p in edge})))


def explicit_graph(sequences):
    lookup = {x: i for i, x in enumerate(sequences)}
    graph = [[] for _ in sequences]
    for i, sequence in enumerate(sequences):
        for p, old in enumerate(sequence):
            for base in 'ACGU':
                if base != old:
                    target = lookup.get(sequence[:p] + base + sequence[p+1:])
                    if target is not None:
                        graph[i].append(target)
    labels = {}
    for start in range(len(sequences)):
        if start in labels:
            continue
        labels[start] = start
        todo = [start]
        for u in todo:
            for v in graph[u]:
                if v not in labels:
                    labels[v] = start
                    todo.append(v)
    return graph, labels


def audit_configuration(n, occurrences, edges=(), fixed=None, kind='physical'):
    fixed = fixed or {}
    groups = defaultdict(list)
    for sequence in map(''.join, itertools.product('ACGU', repeat=n)):
        COUNTS['direct_sequence_assignments'] += 1
        if any(sequence[p] != base for p, base in fixed.items()) or any(sequence[i] + sequence[j] not in PAIR for i, j in edges):
            continue
        peptide = peptide_direct(sequence, occurrences)
        if peptide is not None:
            groups[peptide].append(sequence)
    row = Counter()
    for peptide, sequences in groups.items():
        COUNTS['nonempty_fixed_peptide_instances'] += 1
        graph, labels = explicit_graph(sequences)
        sectors = defaultdict(list)
        for i, sequence in enumerate(sequences):
            sectors[signature_direct(sequence, occurrences, peptide, edges)].append(i)
        predicted_partition = defaultdict(set)
        for signature, indices in sectors.items():
            compiled = compile_sector(n, occurrences, peptide, edges, fixed, signature)
            assert not compiled['unsat']
            generated = set()
            for bits in itertools.product((0, 1), repeat=len(compiled['variable_sites'])):
                if not accepts(compiled, bits):
                    continue
                buckets = [compiled['buckets'][p][bits[compiled['site_variables'][p]]]
                           if p in compiled['site_variables'] else compiled['buckets'][p][0]
                           for p in range(n)]
                fibre = set(map(''.join, itertools.product(*buckets)))
                for sequence in fibre:
                    assert peptide_direct(sequence, occurrences) == peptide
                    assert signature_direct(sequence, occurrences, peptide, edges) == signature
                    assert all(sequence[i] + sequence[j] in PAIR for i, j in edges)
                    assert all(sequence[p] == base for p, base in fixed.items())
                    assert project(compiled, sequence) == bits
                generated.update(fibre)
                COUNTS['complete_quotient_fibres'] += 1
                COUNTS['nontrivial_coordinate_fibres'] += len(fibre) > 1
            assert generated == {sequences[i] for i in indices}
            frozen = frozen_variables(compiled)
            COUNTS['sectors_with_SCC_variables'] += bool(frozen)
            COUNTS['sectors_with_domain_at_least_three'] += any(len(domain) >= 3 for domain in compiled['domains'])
            for i in indices:
                bits = project(compiled, sequences[i])
                predicted_partition[signature, tuple(bits[v] for v in frozen)].add(i)
            COUNTS['compiled_realizable_sectors'] += 1
        true_partition = defaultdict(set)
        for i, label in labels.items():
            true_partition[label].add(i)
        assert {frozenset(v) for v in predicted_partition.values()} == {frozenset(v) for v in true_partition.values()}
        row['components'] += len(true_partition)
        row['old_signature_extra_components'] += len(true_partition) - len(sectors)
        for start in range(len(sequences)):
            distance, todo = {start: 0}, deque([start])
            while todo:
                u = todo.popleft()
                for v in graph[u]:
                    if v not in distance:
                        distance[v] = distance[u] + 1
                        todo.append(v)
            for end in range(start+1, len(sequences)):
                same = end in distance
                if same:
                    assert distance[end] == sum(a != b for a, b in zip(sequences[start], sequences[end]))
                    COUNTS['connected_endpoint_distances'] += 1
                else:
                    COUNTS['separated_endpoint_pairs'] += 1
                if COUNTS['explicit_path_certificates'] < 1500:
                    path = endpoint_path(n, occurrences, peptide, edges, fixed, sequences[start], sequences[end])
                    assert (path is not None) == same
                    if path is not None:
                        assert path[0] == sequences[start] and path[-1] == sequences[end]
                        assert len(path)-1 == sum(a != b for a, b in zip(path[0], path[-1]))
                        for x, y in zip(path, path[1:]):
                            assert sum(a != b for a, b in zip(x, y)) == 1
                        assert all(peptide_direct(x, occurrences) == peptide for x in path)
                        assert all(all(x[i] + x[j] in PAIR for i, j in edges) for x in path)
                        assert all(all(x[p] == b for p, b in fixed.items()) for x in path)
                    COUNTS['explicit_path_certificates'] += 1
        COUNTS['full_component_partitions'] += 1
    COUNTS['configurations'] += 1
    COUNTS[kind + '_configurations'] += 1
    ROWS.append({'n': n, 'kind': kind, 'occurrences': occurrences, 'edges': edges, 'fixed': fixed,
                 'nonempty_peptides': len(groups), **row})


def physical_occurrences(n, strand, offset):
    positions = list(range(offset, n)) if strand == '+' else list(range(n-1-offset, -1, -1))
    return [(tuple(positions[i:i+3]), strand == '-') for i in range(0, len(positions)-2, 3)]


def generic_SCC_check():
    candidates = [(i,) for i in (-1, 1, -2, 2)] + [(i, j) for i in (-1, 1) for j in (-2, 2)]
    for mask in range(1 << len(candidates)):
        compiled = {'unsat': False, 'variable_sites': [0, 1], 'clauses': [c for k, c in enumerate(candidates) if mask & (1 << k)]}
        solutions = [bits for bits in itertools.product((0, 1), repeat=2) if accepts(compiled, bits)]
        graph = [[j for j, y in enumerate(solutions) if sum(a != b for a, b in zip(x, y)) == 1] for x in solutions]
        labels = {}
        for i in range(len(solutions)):
            if i in labels:
                continue
            labels[i] = i
            todo = [i]
            for u in todo:
                for v in graph[u]:
                    if v not in labels:
                        labels[v] = i
                        todo.append(v)
        frozen = frozen_variables(compiled)
        for i, x in enumerate(solutions):
            for j, y in enumerate(solutions):
                assert (labels[i] == labels[j]) == all(x[v] == y[v] for v in frozen)
        COUNTS['generic_two_variable_formulas'] += 1


def main():
    if not __debug__:
        raise RuntimeError('verification requires assertions')
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    out = (ROOT / args.outdir).resolve()
    if ROOT not in out.parents or out.exists():
        raise ValueError('fresh output directory inside repository required')
    out.mkdir(parents=True)
    source_paths = ['analysis/programme/overlap_v2_quotient.py', 'analysis/programme/overlap_v2_verify.py',
                    'analysis/programme/overlap_classification_attempt_v2.json', 'analysis/exact/independent_oracle.py']
    before = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths}
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()
    generic_SCC_check()
    for aa in sorted(set(FORWARD.values())):
        for complement in (False, True):
            for family in ((False, True) if aa == 'S' else (None,)):
                occurrence = ((0, 1, 2), complement)
                signature = ((family,) if family is not None else (), ())
                compiled = compile_sector(3, (occurrence,), (aa,), (), {}, signature)
                for sequence in map(''.join, itertools.product('ACGU', repeat=3)):
                    expected = FORWARD.get(read_direct(sequence, occurrence)) == aa
                    if family is not None:
                        expected = expected and read_direct(sequence, occurrence).startswith('AG') == family
                    inside = all(sequence[p] in domain for p, domain in enumerate(compiled['domains']))
                    actual = inside and accepts(compiled, project(compiled, sequence))
                    assert actual == expected
                    COUNTS['local_codon_family_truth_assignments'] += 1
                COUNTS['local_codon_family_relations'] += 1
    # All directed/complemented pairs of distinct-position triplets on three sites,
    # modulo exchanging the two obligations. Every simple pair graph on3 sites.
    types = [(p, complement) for p in itertools.permutations(range(3)) for complement in (False, True)]
    all_edges = list(itertools.combinations(range(3), 2))
    for first, second in itertools.combinations_with_replacement(types, 2):
        for mask in range(8):
            edges = tuple(e for i, e in enumerate(all_edges) if mask & (1 << i))
            audit_configuration(3, (first, second), edges, kind='arbitrary_triplet')
    print('all3-site ordered-triplet geometries checked', flush=True)
    # All five nonidentical physical relative frames at n4..6, no pair or
    # every one-edge obligation. This is a declared short finite universe.
    for n in (4, 5, 6):
        for strand, offset in (('+', 1), ('+', 2), ('-', 0), ('-', 1), ('-', 2)):
            occurrences = physical_occurrences(n, '+', 0) + physical_occurrences(n, strand, offset)
            for edges in [()] + [(edge,) for edge in itertools.combinations(range(n), 2)]:
                audit_configuration(n, occurrences, edges)
    print('all short physical two-frame geometries checked', flush=True)
    seed = 2026102501
    rng = random.Random(seed)
    for trial in range(60):
        n = 4 + trial % 3
        occurrences = [(tuple(rng.sample(range(n), 3)), bool(rng.randrange(2))) for _ in range(2 + trial % 4)]
        fixed = {p: rng.choice('ACGU') for p in range(n) if rng.random() < .2}
        edges = tuple(rng.sample(list(itertools.combinations(range(n), 2)), trial % 4))
        audit_configuration(n, occurrences, edges, fixed, 'synthetic_pinned')
    # Abstract reversal is not reverse-complement physical translation.
    occurrences = (((0, 1, 2), False), ((2, 1, 0), False))
    assert endpoint_path(3, occurrences, ('R', 'R'), (), {}, 'CGC', 'AGA') is None
    COUNTS['abstract_extra_invariant_fixture'] = 1
    chains = []
    for k in (2, 3, 8, 20):
        occurrences = [(tuple((i, k, i+1)), False) for i in range(k-1)] + [((k-1, k, k+1), False)]
        fixed = {k: 'G', k+1: 'A'}
        path = endpoint_path(k+2, occurrences, ('R',)*k, (), fixed, 'C'*k+'GA', 'A'*k+'GA')
        assert path is not None and len(path) == k+1
        assert all(peptide_direct(x, occurrences) == ('R',)*k for x in path)
        assert [next(p for p in range(k) if x[p] != y[p]) for x, y in zip(path, path[1:])] == list(reversed(range(k)))
        chains.append({'events': k, 'path_length': len(path)-1, 'order': list(reversed(range(k)))})
    COUNTS['unbounded_chain_fixtures'] = len(chains)
    # Domain-fibre boundaries and invalid self-pair/codon-coordinate inputs.
    audit_configuration(0, (), kind='empty')
    audit_configuration(1, (), fixed={0: 'G'}, kind='pinned_untranslated')
    for n, occurrences, peptide, edges in [(3, [((0, 0, 1), False)], ['A'], ()), (3, [], [], ((0, 0),))]:
        try:
            compile_sector(n, occurrences, peptide, edges, {}, ((), (True,) if edges else ()))
        except ValueError:
            pass
        else:
            raise AssertionError('invalid repeated position accepted')
    COUNTS['invalid_input_cases'] = 2
    payload = {'claim_id': 'OVERLAP_V2_QUOTIENT_AND_LIFTING', 'scope': 'standard code, distinct ordered codon coordinates, optional whole-codon complements, fixed bases, six canonical/wobble pair types',
               'counts': dict(COUNTS), 'configurations': ROWS, 'long_chain_fixtures': chains,
               'nonclaims': ['Finite checks do not prove universal geometry.', 'Arbitrary ordered triples need not be physical reading frames.',
                             'No nonstandard code, fitness or molecular intervention claim.']}
    result_path = out/'overlap_v2_results.json'
    result_path.write_text(json.dumps(payload, indent=2)+'\n')
    after = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths}
    if before != after:
        raise RuntimeError('source changed during verification')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    import Bio
    manifest = {'schema_version': 2, 'claim_id': payload['claim_id'], 'repository': {'commit': revision, 'dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))},
                'command': [sys.executable, 'analysis/programme/overlap_v2_verify.py', '--outdir', str(out.relative_to(ROOT))],
                'environment': {'software': [{'name': 'CPython', 'version': platform.python_version()}, {'name': 'Biopython', 'version': Bio.__version__}], 'hardware': platform.platform()},
                'mathematics': {'assertion_tested': 'Exact quotient fibres, complete SCC component partitions and Hamming distances in declared finite overlapping-code universes',
                                'coefficient_domain': 'Boolean and finite RNA alphabet, exact graph integers', 'conventions': payload['scope'],
                                'inputs': source_paths, 'bounds': dict(COUNTS), 'non_claims': payload['nonclaims']},
                'randomness': {'used': True, 'generator': 'Python MT19937', 'seed': seed},
                'run': {'started_at': started, 'runtime_seconds': time.perf_counter()-clock, 'exit_status': 0, 'status': 'completed', 'resource_limits': {'memory_bytes': None, 'cpu_seconds': None, 'max_cores': None, 'max_threads': None}},
                'input_artifacts': [{'path': p, 'sha256': before[p], 'sha256_after': after[p]} for p in source_paths],
                'outputs': [{'path': str(result_path.relative_to(ROOT)).replace('\\', '/'), 'sha256': hashlib.sha256(result_path.read_bytes()).hexdigest(), 'kind': 'result'}],
                'checks': [{'name': 'source-hashes-unchanged', 'passed': True}], 'result': 'completed',
                'declared_results': [str(result_path.relative_to(ROOT)).replace('\\', '/')], 'residual_risks': ['No hard resource limits; not proof-assistant verification.']}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'counts': dict(COUNTS), 'seconds': manifest['run']['runtime_seconds']}), flush=True)


if __name__ == '__main__':
    main()
