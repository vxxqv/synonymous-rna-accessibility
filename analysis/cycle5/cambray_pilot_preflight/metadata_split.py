"""Exact analysis-set peptide relatedness and fixed design-family split.

Reads sequence-only inventories. Never opens the phenotype archive.
"""
import collections
import csv
import hashlib
import json
from pathlib import Path
from Bio.Seq import Seq
from graph_rules import corner_index

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent/'cambray_eligibility'
PROTOCOL = ROOT.parent/'cambray_pilot/protocol_draft.md'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def h(payload):
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()

def main():
    inputs = ['sequence_metadata.tsv','annotated_squares.json','components.jsonl','previously_exposed_ids.json']
    assert digest(INPUT/'sequence_metadata.tsv') == '9fcd79a97e3478c68411dee6cb853ad77e249d1d93b11bc4af58d2fcc44d4e7e'
    rows = {r['ID']:r for r in csv.DictReader((INPUT/'sequence_metadata.tsv').open(encoding='utf-8'),delimiter='\t')}
    squares = json.loads((INPUT/'annotated_squares.json').read_text())
    assert len(squares) == 47
    components = [json.loads(line) for line in (INPUT/'components.jsonl').read_text().splitlines()]
    selected_cids = sorted({s['component_index'] for s in squares})
    selected_ids = sorted({v for cid in selected_cids for v in components[cid]['IDs']})
    by_seed = collections.defaultdict(set)
    for identifier in selected_ids:
        by_seed[rows[identifier]['seed']].add(str(Seq(rows[identifier]['sequence']).translate()))
    assert len(by_seed) == 32
    exposed_peptides = sorted({str(Seq(r['sequence']).translate()) for r in rows.values() if r['seed']=='1'})
    assert len(exposed_peptides) == 2803
    seeds = sorted(by_seed,key=int)
    parent = {seed:seed for seed in seeds}
    def find(seed):
        while parent[seed] != seed:
            parent[seed] = parent[parent[seed]]
            seed = parent[seed]
        return seed
    def join(a,b):
        a,b = find(a),find(b)
        if a != b:
            parent[max(a,b,key=int)] = min(a,b,key=int)
    selected_comparisons = 0
    related_pairs = []
    maximum_selected_match = 0
    for i, a in enumerate(seeds):
        for b in seeds[i+1:]:
            matches = []
            for x in sorted(by_seed[a]):
                for y in sorted(by_seed[b]):
                    selected_comparisons += 1
                    count = sum(u==v for u,v in zip(x,y))
                    maximum_selected_match = max(maximum_selected_match,count)
                    if count >= 23: # ceil(.70 * 32), exactly rather than rounded percent.
                        matches.append({'left_peptide':x,'right_peptide':y,'matching_positions':count})
            if matches:
                join(a,b)
                related_pairs.append({'left_seed':a,'right_seed':b,'witnesses':matches})
    exposed_comparisons = 0
    exposed_matches = []
    maximum_exposed_match = 0
    for seed in seeds:
        for peptide in sorted(by_seed[seed]):
            for exposed in exposed_peptides:
                exposed_comparisons += 1
                count = sum(a==b for a,b in zip(peptide,exposed))
                maximum_exposed_match = max(maximum_exposed_match,count)
                if count >= 23:
                    exposed_matches.append({'selected_seed':seed,'selected_peptide':peptide,'seed1_peptide':exposed,'matching_positions':count})
    excluded_seeds = sorted({r['selected_seed'] for r in exposed_matches},key=int)
    groups = collections.defaultdict(list)
    for seed in seeds:
        groups[find(seed)].append(seed)
    grouped = [sorted(g,key=int) for g in groups.values()]
    if excluded_seeds:
        # The draft does not specify reallocation after new relatedness exclusions.
        # Preserve all records and return the issue before making that choice.
        assignment = {}
        status = 'hold_relatedness_exclusions_require_root_split_rule'
    elif len(grouped) == 32:
        ordered = sorted(seeds,key=lambda s:h('C5-joint-v1|seed='+s))
        assignment = {s:('development' if i<16 else 'reserved_validation') for i,s in enumerate(ordered)}
        status = 'metadata_split_ready'
    else:
        ordered_groups = sorted(grouped,key=lambda g:h('C5-joint-v1|group='+','.join(g)))
        assignment = {s:('development' if i%2==0 else 'reserved_validation') for i,g in enumerate(ordered_groups) for s in g}
        status = 'metadata_split_ready'
    units = []
    for index, square in enumerate(squares):
        ids = square['IDs']
        strings = [rows[v]['sequence'] for v in ids]
        positions = [i for i,col in enumerate(zip(*strings)) if len(set(col))>1]
        assert len(positions)==2
        alleles = [sorted({s[p] for s in strings}) for p in positions]
        assert all(len(a)==2 for a in alleles)
        coordinates = {}
        for identifier, sequence in zip(ids,strings):
            bit0,bit1 = [alleles[j].index(sequence[p]) for j,p in enumerate(positions)]
            coordinates[2*bit0+bit1] = identifier
        assert len(coordinates)==4
        lexical = sorted(ids)
        withheld = lexical[corner_index(lexical)]
        withheld_coordinate = next(i for i,v in coordinates.items() if v==withheld)
        seed = square['seed_ids'][0]
        units.append({'square_index':index,'seed':seed,'arm':assignment.get(seed,'unassigned'),
            'excluded_by_relatedness':seed in excluded_seeds,'peptide':square['peptide'],
            'component_index':square['component_index'],'component_size':square['component_size'],
            'coordinate_order_00_01_10_11':[coordinates[i] for i in range(4)],
            'positions_1based':[p+1 for p in positions],
            'diagonals':[[coordinates[0],coordinates[3]],[coordinates[1],coordinates[2]]],
            'withheld_ID':withheld,'withheld_coordinate':withheld_coordinate,
            'prediction_classification_endpoints':[coordinates[withheld_coordinate^1],coordinates[withheld_coordinate^2]],
            'isolated_observed_square':square['component_size']==4})
    assert len(units)==47 and sum(len(u['diagonals']) for u in units)==94
    if assignment:
        for key in ['peptide','component_index']:
            key_arms = collections.defaultdict(set)
            for unit in units:key_arms[unit[key]].add(unit['arm'])
            assert all(len(arms)==1 for arms in key_arms.values())
        vertex_arms = collections.defaultdict(set)
        for unit in units:
            for identifier in components[unit['component_index']]['IDs']:
                vertex_arms[identifier].add(unit['arm'])
        assert all(len(arms)==1 for arms in vertex_arms.values())
    arms = {}
    for arm in ('development','reserved_validation'):
        subset = [u for u in units if u['arm']==arm]
        arm_cids = {u['component_index'] for u in subset}
        arms[arm] = {'seeds':sorted({u['seed'] for u in subset},key=int),'squares':len(subset),'diagonal_pairs':2*len(subset),
            'peptides':len({u['peptide'] for u in subset}),'observed_components':len(arm_cids),
            'square_vertices':len({v for u in subset for v in u['coordinate_order_00_01_10_11']}),
            'component_vertices':len({v for cid in arm_cids for v in components[cid]['IDs']}),
            'isolated_observed_squares':sum(u['isolated_observed_square'] for u in subset),
            'sequence_upper_bound_can_meet_coverage_minimum':len(subset)>=10 and len({u['seed'] for u in subset})>=8,
            'finite_measurement_coverage':'unexamined'}
    result = {'status':status,'input_hashes':{name:digest(INPUT/name) for name in inputs},'protocol_draft_sha256':digest(PROTOCOL),
        'selection':{'squares':47,'diagonals':94,'square_vertices':181,'observed_components':len(selected_cids),'all_component_vertices':len(selected_ids)},
        'identity_rule':{'denominator':32,'required_matches':23,'nominal_fraction':'.70','selection_scope':'Analyzed square/component peptides only; all seed1 peptides used as exposed reference. Not universal phylogenetic independence.'},
        'selected_cross_seed_peptide_comparisons':selected_comparisons,'exposed_peptide_comparisons':exposed_comparisons,
        'maximum_selected_cross_seed_matching_positions':maximum_selected_match,'maximum_seed1_matching_positions':maximum_exposed_match,
        'related_seed_pairs':related_pairs,'exposed_matches':exposed_matches,'excluded_seeds':excluded_seeds,
        'seed_groups':sorted(grouped,key=lambda g:int(g[0])),'seed_assignment':assignment,'arms':arms,'units':units,
        'hash_conventions':{'encoding':'UTF-8','seed':'C5-joint-v1|seed=<decimal ID>','group':'C5-joint-v1|group=<ascending numeric IDs separated by comma>',
                            'corner':'C5-corner-v1|<lexically sorted four IDs separated by |>; fullSHA256 as unsigned big-endian integer modulo4, indexing lexical list'},
        'phenotype_values_decoded':0,'prediction_values_calculated_on_real_data':0}
    (ROOT/'split_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','selection','selected_cross_seed_peptide_comparisons','exposed_peptide_comparisons','maximum_selected_cross_seed_matching_positions','maximum_seed1_matching_positions','excluded_seeds','arms')},indent=2))

if __name__=='__main__':main()
