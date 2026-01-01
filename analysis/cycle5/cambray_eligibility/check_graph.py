"""Independent string-pair oracle. Reads sequence-only inventory and audit outputs."""
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
import time
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent

def main():
    started = time.perf_counter()
    rows = list(csv.DictReader((ROOT / 'sequence_metadata.tsv').open(encoding='utf-8'), delimiter='\t'))
    groups = collections.defaultdict(list)
    by_id = {r['ID']: r for r in rows}
    for r in rows:
        groups[str(Seq(r['sequence']).translate())].append(r)
    edges = set()
    adjacency = collections.defaultdict(set)
    pairs_checked = 0
    for group in groups.values():
        for a, b in itertools.combinations(group, 2):
            pairs_checked += 1
            differences = []
            for i, (x, y) in enumerate(zip(a['sequence'], b['sequence'])):
                if x != y:
                    differences.append(i)
                    if len(differences) > 1:
                        break
            if len(differences) == 1:
                lo, hi = sorted((a['ID'], b['ID']))
                edges.add((lo, hi, differences[0] + 1))
                adjacency[lo].add(hi)
                adjacency[hi].add(lo)
    stored_edges = set()
    with (ROOT / 'synonymous_edges.tsv').open(encoding='utf-8') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            lo, hi = sorted((row['left_ID'], row['right_ID']))
            stored_edges.add((lo, hi, int(row['position_1based'])))
    assert edges == stored_edges
    # Distinct representation: enumerate opposite vertex pairs of all cycles of
    # length four, then retain only exact binary coordinate squares.
    cycles = set()
    for a in adjacency:
        distance_two = set().union(*(adjacency[b] for b in adjacency[a]))
        for c in distance_two:
            if a >= c:
                continue
            for b, d in itertools.combinations(adjacency[a] & adjacency[c], 2):
                ids = tuple(sorted((a, b, c, d)))
                if len(set(ids)) != 4:
                    continue
                sequences = [by_id[i]['sequence'] for i in ids]
                choices = [set(col) for col in zip(*sequences)]
                variable = [col for col in choices if len(col) > 1]
                if len(variable) == 2 and all(len(col) == 2 for col in variable):
                    cycles.add(ids)
    stored_faces = [json.loads(line) for line in (ROOT / 'binary_faces.jsonl').read_text().splitlines()]
    assert cycles == {tuple(sorted(r['IDs'])) for r in stored_faces if r['dimension'] == 2}
    # A binary cube must contain a binary square and a disjoint parallel square.
    # Try every one-coordinate translation of each found square by dictionary
    # membership, independent of the original corner-expansion algorithm.
    seq_id = {r['sequence']: r['ID'] for r in rows}
    cubes = set()
    for face in cycles:
        seqs = [by_id[i]['sequence'] for i in face]
        for pos, bases in enumerate(zip(*seqs)):
            if len(set(bases)) != 1:
                continue
            for base in set('ACGT') - {bases[0]}:
                translated = [s[:pos] + base + s[pos+1:] for s in seqs]
                if all(s in seq_id for s in translated):
                    if all(str(Seq(s).translate()) == str(Seq(seqs[0]).translate()) for s in translated):
                        cubes.add(tuple(sorted((*face, *(seq_id[s] for s in translated)))))
    assert not cubes
    seen = set()
    component_sizes = collections.Counter()
    for identifier in by_id:
        if identifier in seen:
            continue
        pending = [identifier]
        seen.add(identifier)
        size = 0
        while pending:
            current = pending.pop()
            size += 1
            for neighbor in adjacency[current] - seen:
                seen.add(neighbor)
                pending.append(neighbor)
        component_sizes[size] += 1
    expected = json.loads((ROOT / 'graph_summary.json').read_text())
    assert dict(component_sizes) == {int(k): v for k, v in expected['synonymous_component_size_counts'].items()}
    receipt = {
        'status': 'passed_distinct_algorithm_corroboration',
        'same_author_not_separate_person_review': True,
        'same_peptide_string_pairs_checked': pairs_checked,
        'exact_synonymous_edge_set_count': len(edges),
        'exact_binary_square_set_count': len(cycles),
        'binary_cube_count': len(cubes),
        'component_size_distribution_exact_match': True,
        'phenotype_values_decoded': 0,
        'elapsed_seconds': time.perf_counter() - started,
        'hashes': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in ['sequence_metadata.tsv', 'synonymous_edges.tsv', 'binary_faces.jsonl', 'graph_summary.json', 'audit_graph.py', 'check_graph.py']},
        'limits': 'Exact synonymous edges, binary squares, absence of cubes and components corroborated. Raw nonsynonymous edge count is not independently re-enumerated. Finite phenotype coverage remains unexamined.'
    }
    (ROOT / 'graph_corroboration.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
