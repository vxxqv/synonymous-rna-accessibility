"""Sequence-selected square dependencies and elementary endpoint identifiability."""
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent

def reachable(allowed, start=0, target=3):
    if start not in allowed or target not in allowed:
        return False
    seen = {start}
    pending = [start]
    while pending:
        v = pending.pop()
        for w in (v ^ 1, v ^ 2):
            if w in allowed and w not in seen:
                seen.add(w)
                pending.append(w)
    return target in seen

def main():
    rows = {r['ID']: r for r in csv.DictReader((ROOT/'sequence_metadata.tsv').open(encoding='utf-8'), delimiter='\t')}
    components = [json.loads(s) for s in (ROOT/'components.jsonl').read_text().splitlines()]
    id_component = {v:i for i,c in enumerate(components) for v in c['IDs']}
    faces = [json.loads(s) for s in (ROOT/'binary_faces.jsonl').read_text().splitlines()]
    seeds = collections.Counter(); peptides = collections.Counter(); comp = collections.Counter(); vertices = set(); contexts = collections.Counter()
    annotated = []
    for face in faces:
        ids = face['IDs']
        strings = [rows[v]['sequence'] for v in ids]
        peptide = str(Seq(strings[0]).translate())
        assert all(str(Seq(s).translate()) == peptide for s in strings)
        positions = [i for i,c in enumerate(zip(*strings)) if len(set(c)) > 1]
        assert len(positions) == 2
        context = ''.join('*' if i in positions else c for i,c in enumerate(strings[0]))
        family = {rows[v]['seed'] for v in ids}
        assert len(family) == 1
        cset = {id_component[v] for v in ids}
        assert len(cset) == 1
        seed = next(iter(family)); cid = next(iter(cset))
        seeds[seed] += 1; peptides[peptide] += 1; comp[cid] += 1; contexts[context] += 1; vertices.update(ids)
        annotated.append({**face, 'peptide':peptide, 'component_index':cid, 'component_size':components[cid]['size'], 'variable_positions_1based':[p+1 for p in positions], 'same_codon':positions[0]//3 == positions[1]//3, 'context_sha256':hashlib.sha256(context.encode()).hexdigest()})
    patterns = []
    for bits in itertools.product((False, True), repeat=4):
        a,b,c,d = bits
        trait1 = {0,3} | ({1} if a else set()) | ({2} if b else set())
        trait2 = {0,3} | ({1} if c else set()) | ({2} if d else set())
        detected = reachable(trait1) and reachable(trait2) and not reachable(trait1 & trait2)
        formula = (a or b) and (c or d) and not ((a and c) or (b and d))
        assert detected == formula
        patterns.append({'intermediate_trait_bits':bits, 'joint_only_obstruction':detected})
    assert sum(r['joint_only_obstruction'] for r in patterns) == 2
    # Joint obstruction can occur for two exact additive score fields on Q2.
    trait1 = [1,0,2,1]; trait2 = [1,2,0,1]
    assert trait1[0]+trait1[3] == trait1[1]+trait1[2]
    assert trait2[0]+trait2[3] == trait2[1]+trait2[2]
    assert reachable({i for i,x in enumerate(trait1) if x >= 1})
    assert reachable({i for i,x in enumerate(trait2) if x >= 1})
    assert not reachable({i for i in range(4) if min(trait1[i],trait2[i]) >= 1})
    summary = {
        'squares':len(faces), 'unique_square_vertices':len(vertices),
        'design_seed_families':len(seeds), 'square_counts_by_seed':dict(seeds),
        'peptide_families':len(peptides), 'square_counts_per_peptide_frequency':dict(collections.Counter(peptides.values())),
        'synonymous_components':len(comp), 'square_counts_per_component_frequency':dict(collections.Counter(comp.values())),
        'distinct_two_coordinate_contexts':len(contexts),
        'same_codon_squares':sum(r['same_codon'] for r in annotated),
        'different_codon_squares':sum(not r['same_codon'] for r in annotated),
        'squares_whose_entire_observed_component_has_four_vertices':sum(r['component_size']==4 for r in annotated),
        'all_seed1_excluded':all(not r['excluded_by_seed1_peptide_boundary'] for r in faces),
        'square_independence_limit':'Different squares may share vertices, component, peptide, design seed and experimental preparation. Counts are not counts of independent biological replicates.',
        'full_context_status':'Exact intended variable-region translation checked. Shared fixed context and codon-aligned insertion are authenticated by original methods and author code. Full physical clone genomes/backbones were not sequenced here, and Addgene full-sequence download requires login. Conditional full-ORF equality follows for any common codon-aligned fixed prefix and suffix.',
        'synthetic_truth_table':patterns,
        'additive_positive_control':{'trait1':trait1,'trait2':trait2,'threshold':1,'joint_only_obstruction':True},
        'identifiability_limit':'The joint-only endpoint obstruction is identifiable on a complete square, but it does not identify epistasis or an RNA mechanism. Exact additive opposing effects already suffice.',
        'phenotype_values_decoded':0
    }
    (ROOT/'square_dependencies.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (ROOT/'annotated_squares.json').write_text(json.dumps(annotated,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__': main()
