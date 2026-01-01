"""Exact sequence-only audit of the authenticated Cambray design inventory."""
import ast
import collections
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
SOURCE_SHA = "9fcd79a97e3478c68411dee6cb853ad77e249d1d93b11bc4af58d2fcc44d4e7e"


def main():
    started = time.perf_counter()
    src = HERE / "sequence_metadata.tsv"
    if hashlib.file_digest(src.open("rb"), "sha256").hexdigest() != SOURCE_SHA:
        raise ValueError("Sequence inventory changed")
    rows = list(csv.DictReader(src.open(encoding="utf-8"), delimiter="\t"))
    if len({r['ID'] for r in rows}) != len(rows):
        raise ValueError("Duplicate identifiers")
    if any(r['ID'] != '_'.join([r['seed'], r['combi'], r['design_rep']]) for r in rows):
        raise ValueError("Identifier and design columns disagree")
    author_source = HERE.parents[1] / "literature_expansion/cambray_generate_seed.py"
    # The deposited script is Python 2. Extract only its literal code table,
    # rather than executing or parsing its unrelated generator implementation.
    table_literal = author_source.read_text(encoding='utf-8').split('code = ',1)[1].split('\n\n',1)[0]
    code = {c.upper(): a.upper() for c, a in ast.literal_eval(table_literal).items()}
    # Independent runtime table checks the extracted author code, including stops.
    from Bio.Data import CodonTable
    reference = CodonTable.unambiguous_dna_by_id[1]
    independent_code = dict(reference.forward_table)
    independent_code.update({c: '*' for c in reference.stop_codons})
    if code != independent_code:
        raise AssertionError("Author genetic code differs from standard code")
    seq_records = collections.defaultdict(list)
    seeds = collections.Counter()
    seed_combi = collections.defaultdict(collections.Counter)
    for r in rows:
        s = r['sequence']
        if len(s) != 96 or set(s) - set('ACGT'):
            raise ValueError("Noncanonical sequence or length")
        seq_records[s].append(r)
        seeds[r['seed']] += 1
        seed_combi[r['seed']][r['combi']] += 1
    seqs = sorted(seq_records)
    aa = [''.join(code[s[i:i+3]] for i in range(0,96,3)) for s in seqs]
    if any('*' in a for a in aa):
        raise ValueError("Stop codon in intended frame")
    groups = collections.defaultdict(list)
    for i, peptide in enumerate(aa):
        groups[peptide].append(i)
    trans = str.maketrans({'A':'0','C':'1','G':'2','T':'3'})
    enc = [int(s.translate(trans),4) for s in seqs]
    lookup = {v:i for i,v in enumerate(enc)}
    adjacency = [[] for s in seqs]
    raw_edges = 0
    synonymous = []
    cross_seed = 0
    for i, value in enumerate(enc):
        for pos in range(96):
            shift = 2*(95-pos)
            current = (value >> shift) & 3
            for alternate in range(current+1,4):
                candidate = value ^ ((current ^ alternate) << shift)
                j = lookup.get(candidate)
                if j is None:
                    continue
                raw_edges += 1
                if aa[i] == aa[j]:
                    synonymous.append((i,j,pos))
                    adjacency[i].append((j,pos))
                    adjacency[j].append((i,pos))
                    if {r['seed'] for r in seq_records[seqs[i]]} != {r['seed'] for r in seq_records[seqs[j]]}:
                        cross_seed += 1
        if i and i%50000==0:
            print(json.dumps({'processed':i,'unique_sequences':len(seqs),'seconds':round(time.perf_counter()-started,2)}),flush=True)
    components=[]
    seen=set()
    for start in range(len(seqs)):
        if start in seen:
            continue
        stack=[start];seen.add(start);members=[]
        while stack:
            v=stack.pop();members.append(v)
            for w,_ in adjacency[v]:
                if w not in seen:seen.add(w);stack.append(w)
        components.append(sorted(members))
    faces2=[];faces3=[];two_steps=0;endpoint_pairs=set()
    for i,value in enumerate(enc):
        neighbors=adjacency[i]
        for (a,p),(b,q) in itertools.combinations(neighbors,2):
            if p != q:
                two_steps += 1
                endpoint_pairs.add((min(a,b),max(a,b)))
        # The minimum nucleotide word uniquely roots a binary face.
        upper=[(j,p,enc[j]^value) for j,p in neighbors if enc[j]>value]
        squares={}
        for a,b in itertools.combinations(range(len(upper)),2):
            j,p,delta1=upper[a];k,q,delta2=upper[b]
            if p==q:continue
            last=lookup.get(value^delta1^delta2)
            if last is not None and aa[last]==aa[i]:
                vertices=tuple(sorted([i,j,k,last]))
                faces2.append(vertices);squares[(a,b)]=last
        for a,b,c in itertools.combinations(range(len(upper)),3):
            if (a,b) not in squares or (a,c) not in squares or (b,c) not in squares:continue
            last=lookup.get(value^upper[a][2]^upper[b][2]^upper[c][2])
            if last is not None and aa[last]==aa[i]:
                faces3.append(tuple(sorted([i,upper[a][0],upper[b][0],upper[c][0],squares[(a,b)],squares[(a,c)],squares[(b,c)],last])))
    if len(faces2)!=len(set(faces2)) or len(faces3)!=len(set(faces3)):
        raise AssertionError("Duplicate binary face")
    # A peptide's entire sampled fibre is a complete product of its observed
    # per-position nucleotide alleles iff the product cardinality equals its size.
    complete_groups=[]
    for peptide,members in groups.items():
        if len(members)<2:continue
        choices=[set(seqs[i][p] for i in members) for p in range(96)]
        total=math.prod(len(a) for a in choices)
        if total==len(members):
            complete_groups.append({'peptide':peptide,'vertices':len(members),'variable_positions_1based':[p+1 for p,a in enumerate(choices) if len(a)>1],
                                    'factor_sizes':[len(a) for a in choices if len(a)>1],'IDs':[seq_records[seqs[i]][0]['ID'] for i in members]})
    seed1peptides={aa[i] for i,s in enumerate(seqs) if any(r['seed']=='1' for r in seq_records[s])}
    excluded={i for i,s in enumerate(seqs) if aa[i] in seed1peptides}
    with (HERE/'synonymous_edges.tsv').open('x',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['left_ID','right_ID','position_1based'])
        for i,j,p in synonymous:w.writerow([seq_records[seqs[i]][0]['ID'],seq_records[seqs[j]][0]['ID'],p+1])
    with (HERE/'binary_faces.jsonl').open('x',encoding='utf-8') as f:
        for dimension,faces in [(2,faces2),(3,faces3)]:
            for vertices in faces:
                f.write(json.dumps({'dimension':dimension,'IDs':[seq_records[seqs[i]][0]['ID'] for i in vertices],
                                    'seed_ids':sorted({r['seed'] for i in vertices for r in seq_records[seqs[i]]}),
                                    'excluded_by_seed1_peptide_boundary':bool(set(vertices)&excluded)})+'\n')
    with (HERE/'components.jsonl').open('x',encoding='utf-8') as f:
        for members in components:
            if len(members)>1:
                f.write(json.dumps({'size':len(members),'IDs':[seq_records[seqs[i]][0]['ID'] for i in members],
                                   'peptide':aa[members[0]],'seed_ids':sorted({r['seed'] for i in members for r in seq_records[seqs[i]]})})+'\n')
    (HERE/'complete_observed_allele_fibres.json').write_text(json.dumps(complete_groups,indent=2)+'\n',encoding='utf-8')
    duplicate_records=[{'sequence_sha256':hashlib.sha256(s.encode()).hexdigest(),'IDs':[r['ID'] for r in rec]} for s,rec in seq_records.items() if len(rec)>1]
    (HERE/'duplicate_sequences.json').write_text(json.dumps(duplicate_records,indent=2)+'\n',encoding='utf-8')
    summary={'status':'sequence_only_graph_audit_complete','records':len(rows),'unique_ids':len(rows),'unique_96nt_sequences':len(seqs),
             'duplicate_sequence_groups':len(duplicate_records),'unique_32aa_peptides':len(groups),'stop_containing_sequences':0,
             'seed_count':len(seeds),'seed_record_counts':dict(sorted(seeds.items(),key=lambda x:int(x[0]))),
             'design_rep_counts':dict(collections.Counter(r['design_rep'] for r in rows)),
             'seed_combination_counts':{s:len(v) for s,v in seed_combi.items()},
             'raw_hamming1_edges':raw_edges,'synonymous_hamming1_edges':len(synonymous),
             'cross_seed_synonymous_edges':cross_seed,'vertices_incident_to_synonymous_edges':sum(bool(a) for a in adjacency),
             'synonymous_component_size_counts':dict(sorted(collections.Counter(map(len,components)).items())),
             'maximum_component_size':max(map(len,components)),'two_step_routes_with_hamming2_endpoints':two_steps,
             'distinct_hamming2_endpoint_pairs_with_measured_design_route':len(endpoint_pairs),
             'complete_binary_2_faces':len(faces2),'complete_binary_3_faces':len(faces3),
             'complete_nontrivial_observed_allele_peptide_fibres':len(complete_groups),
             'complete_observed_allele_fibre_size_counts':dict(sorted(collections.Counter(g['vertices'] for g in complete_groups).items())),
             'seed1_records':seeds['1'],'seed1_peptide_count':len(seed1peptides),'seed1_peptide_excluded_unique_sequences':len(excluded),
             'unexposed_boundary_synonymous_edges':sum(i not in excluded and j not in excluded for i,j,p in synonymous),
             'unexposed_boundary_complete_binary_2_faces':sum(not(set(v)&excluded) for v in faces2),
             'unexposed_boundary_complete_binary_3_faces':sum(not(set(v)&excluded) for v in faces3),
             'phenotype_values_decoded':0,'joint_phenotype_finite_coverage':'not examined before protocol',
             'scope':'Edges and faces connect deposited design sequences, not verified finite paired phenotypes. Complete binary faces hold all other nucleotide positions fixed; they are not complete unrestricted synonymous genotype spaces.',
             'sequence_metadata_sha256':SOURCE_SHA,'elapsed_seconds':time.perf_counter()-started}
    (HERE/'graph_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
