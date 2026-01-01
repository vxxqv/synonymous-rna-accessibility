"""Independent finite checks, with before/after hashes, for the 2-CNF compiler.

Run from repository root, providing a fresh --outdir within the repository.
This does not rerun the separate 13-million-instance graph audit.
"""
import argparse
from collections import Counter
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

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'analysis/expansion'))
sys.path.insert(0, str(ROOT / 'analysis/exact'))
from feasibility_2sat import BITS, CODON_CNF, compile_instance, encode, evaluate, find_witness, solve_cnf
from independent_oracle import CODE, CODONS, PAIRS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direct_valid(peptide, edges, fixed, sequence):
    return (len(sequence) == 3 * len(peptide)
            and all(CODE[sequence[3 * c:3 * c + 3]] == aa for c, aa in enumerate(peptide))
            and all(sequence[i] + sequence[j] in PAIRS for i, j in edges)
            and all(sequence[p] == base for p, base in fixed.items()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('assertions must remain enabled')
    out = (ROOT / args.outdir).resolve()
    if ROOT not in out.parents or out.exists():
        raise ValueError('fresh output directory inside repository required')
    out.mkdir(parents=True)
    source_paths = ['analysis/expansion/feasibility_2sat.py',
                    'analysis/expansion/verification/verify_feasibility_2sat.py',
                    'analysis/exact/independent_oracle.py']
    before = {p: digest(ROOT / p) for p in source_paths}
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    timer = time.perf_counter()
    counts = Counter()
    seed = 2026102411
    rng = random.Random(seed)
    records = []

    # Every assignment of the six bits is an RNA codon, not only sense codons.
    for aa in sorted(CODONS):
        accepted = []
        for codon in map(''.join, itertools.product('ACGU', repeat=3)):
            bits = encode(codon)
            assert evaluate(CODON_CNF[aa], bits) == (CODE[codon] == aa)
            fixed = dict(enumerate(codon))
            witness = find_witness(aa, (), fixed)
            assert (witness is not None) == (CODE[codon] == aa)
            if witness is not None:
                assert witness == codon
                accepted.append(codon)
            counts['codon_truth_table_cases'] += 1
            counts['fully_pinned_codon_cases'] += 1
            if CODE[codon] == '*':
                counts['stop_codon_rejections'] += 1
        for p in range(3):
            for base in 'ACGU':
                fixed = {p: base}
                witness = find_witness(aa, (), fixed)
                expected = any(codon[p] == base for codon in CODONS[aa])
                assert (witness is not None) == expected
                if expected:
                    assert direct_valid(aa, (), fixed, witness)
                counts['single_pin_cases'] += 1
        records.append({'amino_acid': aa, 'clauses': CODON_CNF[aa], 'accepted_codons': accepted})
    assert len(set(BITS.values())) == 4
    for x, y in itertools.product('ACGU', repeat=2):
        r, d = BITS[x]
        s, e = BITS[y]
        assert ((r != s) and not (d and e)) == (x + y in PAIRS)
        counts['ordered_pair_truth_table_cases'] += 1

    # Independent brute force also tests generic SCC assignment polarity.
    for _ in range(1000):
        n = rng.randint(0, 7)
        clauses = []
        if n:
            for _ in range(rng.randint(0, 25)):
                clauses.append(tuple(rng.choice((-1, 1)) * rng.randint(1, n)
                                     for _ in range(rng.randint(1, 2))))
        expected = next((bits for bits in itertools.product((0, 1), repeat=n)
                         if evaluate(clauses, bits)), None)
        observed = solve_cnf(n, clauses)
        assert (observed is None) == (expected is None)
        if observed is not None:
            assert evaluate(clauses, observed)
        counts['random_generic_cnf_cases'] += 1

    aas = sorted(CODONS)
    examples = []
    for m in range(1, 5):
        for trial in range(250):
            peptide = ''.join(rng.choice(aas) for _ in range(m))
            all_sequences = [''.join(codons) for codons in itertools.product(*(CODONS[a] for a in peptide))]
            all_edges = list(itertools.combinations(range(3 * m), 2))
            start = rng.choice(all_sequences)
            eligible = (all_edges if trial % 2 else
                        [(i, j) for i, j in all_edges if start[i] + start[j] in PAIRS])
            edges = tuple(edge for edge in eligible if rng.random() < 0.22)
            fixed = {p: (start[p] if trial % 3 else rng.choice('ACGU'))
                     for p in range(3 * m) if rng.random() < 0.25}
            direct = [sequence for sequence in all_sequences if direct_valid(peptide, edges, fixed, sequence)]
            witness = find_witness(peptide, edges, fixed)
            assert (witness is not None) == bool(direct)
            if witness is not None:
                assert direct_valid(peptide, edges, fixed, witness)
            n, clauses = compile_instance(peptide, edges, fixed)
            assert n == 6 * m and len(clauses) <= 66 * m + 3 * len(edges) + 2 * len(fixed)
            counts['random_whole_instances'] += 1
            counts['random_whole_satisfiable' if direct else 'random_whole_unsatisfiable'] += 1
            counts['brute_force_sequences'] += len(all_sequences)
            if len(examples) < 8:
                examples.append({'peptide': peptide, 'edges': edges, 'fixed': fixed, 'witness': witness})
    assert find_witness('', ()) == ''
    assert find_witness('L', ((0, 1), (1, 2), (0, 2))) is None  # odd orientation cycle
    assert find_witness('M', (), {0: 'C'}) is None
    assert find_witness('K', ((0, 2),)) is None  # two purines
    assert find_witness('F', ((0, 2),)) is None  # two pyrimidines
    assert find_witness('LK', ((0, 5), (5, 0))) is not None  # duplicate harmless
    witness = find_witness('S' * 200, ())
    assert direct_valid('S' * 200, (), {}, witness)  # no recursion limit dependence
    counts['explicit_boundary_cases'] = 7
    for peptide, edges, fixed in [('*', (), {}), ('A', ((0, 0),), {}),
                                   ('A', ((0, 3),), {}), ('', (), {0: 'A'}),
                                   ('A', (), {0: 'T'})]:
        try:
            compile_instance(peptide, edges, fixed)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid model input accepted')
        counts['invalid_input_cases'] += 1

    payload = {'claim': 'Exact standard-code pair-graph feasibility has a linear-size 2-CNF encoding',
               'seed': seed, 'counts': dict(counts), 'codon_clause_certificate': records,
               'examples': examples,
               'limits': ['Random whole-instance checks cover lengths 1..4, not all lengths.',
                          'The separate structural proof establishes the universal result.',
                          'Boolean bit flips do not represent nucleotide mutation adjacency.']}
    result_path = out / 'feasibility_checks.json'
    result_path.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    after = {p: digest(ROOT / p) for p in source_paths}
    assert before == after
    try:
        revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
        dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = 'unavailable', True
    manifest = {'schema_version': 2, 'claim_id': 'STANDARD_CODE_2CNF_FEASIBILITY',
                'repository': {'commit': revision, 'dirty': dirty},
                'command': [sys.executable, 'analysis/expansion/verification/verify_feasibility_2sat.py', '--outdir', str(out.relative_to(ROOT))],
                'environment': {'software': [{'name': 'CPython', 'version': platform.python_version()}], 'hardware': platform.platform()},
                'mathematics': {'assertion_tested': payload['claim'], 'coefficient_domain': 'Boolean exact finite logic',
                                'conventions': 'A=11,G=10,C=01,U=00; standard sense code; six allowed RNA pairs; arbitrary simultaneous pair graph',
                                'inputs': source_paths, 'bounds': dict(counts), 'non_claims': payload['limits']},
                'randomness': {'used': True, 'generator': 'Python random.Random MT19937', 'seed': seed},
                'run': {'started_at': started, 'runtime_seconds': time.perf_counter() - timer, 'status': 'completed', 'exit_status': 0,
                        'resource_limits': {'memory_bytes': None, 'cpu_seconds': None, 'max_cores': None, 'max_threads': None}},
                'input_artifacts': [{'path': p, 'sha256': before[p], 'sha256_after': after[p]} for p in source_paths],
                'outputs': [{'path': str(result_path.relative_to(ROOT)).replace('\\', '/'), 'sha256': digest(result_path), 'kind': 'result'}],
                'checks': [{'name': 'input-hashes-unchanged', 'passed': True}], 'result': 'completed',
                'residual_risks': ['Finite checking is not proof-assistant verification; no hard resource limits imposed.'],
                'declared_results': [str(result_path.relative_to(ROOT)).replace('\\', '/')]}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'counts': dict(counts), 'seconds': manifest['run']['runtime_seconds']}))


if __name__ == '__main__':
    main()
