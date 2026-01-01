"""Sequence-only reconstruction of the declared existing FAS library.

No observed counts, phenotype values, outcome coverage or fitting inputs.
The native annotation authenticates intended coding consequence, not translation
of the experimental minigene, which was an RNA splicing assay.
"""
from pathlib import Path
import collections
import hashlib
import itertools
import json
import re
from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
NATIVE = HERE / 'native_refseq.gb'
EXPECTED_NATIVE_SHA = 'e8964c4db7d19acd5f6a1075106ec67087ad4bc5a8dc3a986f7eccb0a7f31ca7'
if hashlib.sha256(NATIVE.read_bytes()).hexdigest() != EXPECTED_NATIVE_SHA:
    raise RuntimeError('Native reference changed')
record = SeqIO.read(NATIVE, 'genbank')
wt = 'GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
ancestor = 'GATCCAGATCTAACTTGCTGTGGTTGTGTCTCCTGCTTCTCCCGATTCTAGTAATTGTTTGGG'
# Exact declared-library regex from author's 005_remove_sequencing_errors.sh.
pattern = 'GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG'
tokens = re.findall(r'\[([ACGT]+)\]|([ACGT])', pattern)
domains = [tuple(a or b) for a, b in tokens]
sequences = [''.join(x) for x in itertools.product(*domains)]
if len(sequences) != 3072 or wt not in sequences or ancestor not in sequences:
    raise RuntimeError('Library design mismatch')
position = str(record.seq).find(wt)
if position < 0 or str(record.seq).count(wt) != 1:
    raise RuntimeError('Native reference match not unique')
cds = next(f for f in record.features if f.type == 'CDS')
start, end = int(cds.location.start), int(cds.location.end)
native_cds = str(record.seq[start:end])
offset = position - start
protein = str(Seq(native_cds).translate())
if protein.rstrip('*') != cds.qualifiers['translation'][0]:
    raise RuntimeError('Independent native annotation translation mismatch')

def translated(exon):
    substituted = native_cds[:offset] + exon + native_cds[offset + len(wt):]
    return str(Seq(substituted).translate())

groups = collections.defaultdict(list)
for sequence in sequences:
    groups[translated(sequence)].append(sequence)
human = groups[protein]
hs = set(human)
edges = set()
for sequence in human:
    for i, bases in enumerate(domains):
        for base in bases:
            other = sequence[:i] + base + sequence[i + 1:]
            if other in hs and other != sequence:
                edges.add(tuple(sorted((sequence, other))))
variable = []
for i, bases in enumerate(domains):
    if len(bases) > 1:
        codon_start = ((offset + i) // 3) * 3
        changes = []
        for base in bases:
            other = wt[:i] + base + wt[i + 1:]
            changes.append({'base': base, 'same_native_protein': translated(other) == protein})
        variable.append({'exon_position_1based': i + 1,
                         'native_base': wt[i], 'declared_bases': list(bases),
                         'native_codon': native_cds[codon_start:codon_start+3],
                         'codon_position_1based': (offset + i) % 3 + 1,
                         'human_class_observed_bases': sorted({z[i] for z in human}),
                         'single_change_consequence': changes})

# All single and double substitutions in the 2016 declared local universe.
# This checks design potential only. It does not inspect measured completeness.
single_synonymous = []
for i in range(len(wt)):
    for base in 'ACGT':
        if base != wt[i] and translated(wt[:i]+base+wt[i+1:]) == protein:
            single_synonymous.append((i,base))
two_site_synonymous_pairs = 0
for (i, b), (j, c) in itertools.combinations(single_synonymous, 2):
    if i == j:
        continue
    seq = list(wt); seq[i] = b; seq[j] = c
    if translated(''.join(seq)) == protein:
        two_site_synonymous_pairs += 1

out = {
    'scope': 'Declared library sequence design and native included-isoform coding consequence only. No measured outcome coverage.',
    'native_accession': record.id, 'native_sha256': EXPECTED_NATIVE_SHA,
    'native_exon_start0': position, 'native_cds_start0': start,
    'native_exon_start_cds_offset0': offset, 'native_exon_start_codon_offset0': offset % 3,
    'native_protein_id': cds.qualifiers['protein_id'][0],
    'native_translation_annotation_matches': True,
    'declared_variants': len(sequences),
    'protein_classes': len(groups),
    'protein_class_size_histogram': dict(sorted(collections.Counter(map(len, groups.values())).items())),
    'human_reference_class_vertices': len(human),
    'human_reference_class_edges': len(edges),
    'human_class_hamming_distance_from_native_histogram': dict(sorted(collections.Counter(sum(a != b for a,b in zip(s,wt)) for s in human).items())),
    'human_class_potential_overlap_with_2016_radius2_design': sum(sum(a != b for a,b in zip(s,wt)) <= 2 for s in human),
    'human_reference_class_variable_positions': variable,
    'ancestor_same_protein_as_human': translated(ancestor) == protein,
    '2016_potential_synonymous_singles': len(single_synonymous),
    '2016_potential_complete_same_protein_two_site_squares': two_site_synonymous_pairs,
    '2016_measured_same_protein_squares': None,
    '2019_measured_same_protein_coverage': None,
    'minigene_translation_or_protein_output_authenticated': False,
}
(HERE / 'declared_human_reference_class_sequences.json').write_text(json.dumps({'scope':'Existing published library design only, not measured coverage or outcomes','native_accession':record.id,'sequences':sorted(human)},indent=2)+'\n',encoding='utf8')
(HERE / 'sequence_contract_check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(out,indent=2))
