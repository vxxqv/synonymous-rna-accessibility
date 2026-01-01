"""Bounded sequence/graph dependency audit. No phenotype reader or model fitting."""
if not __debug__:
    raise SystemExit('Refuse disabled assertions before input access')
import argparse
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
from Bio.Seq import Seq

EXPECTED={
    'split':'394a0c1aa23785c290b31c147d199f406cf46b269daa9159a98c419398b05a4e',
    'mask':'6e869e95b50fe21073fabfca5cda0e544b140b8e5c8a22f217b3c52394fdf3d9',
    'components':'90772dc57ab6dfdadf737e72cc0be11783245299d1ec0e9e589663a64ee1dadd',
    'sequences':'9fcd79a97e3478c68411dee6cb853ad77e249d1d93b11bc4af58d2fcc44d4e7e',
    'edges':'5d7d5f3b1160d7a8338da3f39369720fdc92c6105b47efbd902eeede48a30edd'}
RELATIONS=['shared_square_ID','shared_mutation_edge','same_original_masked_target',
    'shared_global_masked_ID_in_square','same_peptide','same_design_seed',
    'same_related_seed_group','same_observed_component','shared_component_ID']

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def check(condition,message):
    if not condition:raise ValueError(message)

def write_json(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')

def write_csv(path,rows):
    with path.open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def stats(counts):
    total=sum(counts);positive=[n for n in counts if n]
    return {'groups':len(counts),'nonempty_groups':len(positive),'total_cases':total,
        'minimum':min(counts),'maximum':max(counts),
        'maximum_fraction':max(counts)/total if total else None,
        'sum_squared_case_fractions':sum(n*n for n in counts)/(total*total) if total else None,
        'count_histogram':dict(sorted(collections.Counter(counts).items()))}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for role in EXPECTED:parser.add_argument('--'+role,required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    check(not args.output.exists(),'Output exists; preserve prior receipt')
    inputs={r:getattr(args,r) for r in EXPECTED}
    for r,p in inputs.items():check(sha(p)==EXPECTED[r],'Input hash changed: '+r)
    split=json.loads(args.split.read_text());mask=json.loads(args.mask.read_text())
    units=split['units'];mu={u['square_index']:u for u in mask['units']}
    check(len(units)==47 and len(mu)==47,'Expected all 47 original squares')
    seed_group={seed:','.join(sorted(g,key=int)) for g in split['seed_groups'] for seed in g}
    components=[json.loads(line) for line in args.components.read_text().splitlines()]
    selected_cids={u['component_index'] for u in units}
    selected_ids={v for cid in selected_cids for v in components[cid]['IDs']}
    rows={};record_count=0
    with args.sequences.open(encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            record_count+=1
            if row['ID'] in selected_ids:
                check(row['ID'] not in rows,'Duplicate selected identity');rows[row['ID']]=row
    check(record_count==244000 and set(rows)==selected_ids,'Sequence identity coverage changed')
    all_edges=set();adj={v:set() for v in selected_ids}
    with args.edges.open(encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            a,b=row['left_ID'],row['right_ID'];edge=tuple(sorted((a,b)));all_edges.add(edge)
            if a in adj:adj[a].add(b)
            if b in adj:adj[b].add(a)
    check(len(all_edges)==20031,'Global observed edge inventory changed')
    for cid in selected_cids:
        ids=set(components[cid]['IDs']);seen={next(iter(ids))};todo=list(seen)
        while todo:
            for v in adj[todo.pop()]:
                check(v in ids,'Stored component omits an adjacent vertex')
                if v not in seen:seen.add(v);todo.append(v)
        check(seen==ids,'Stored component is disconnected')
        check(len({str(Seq(rows[v]['sequence']).translate()) for v in ids})==1,'Component translation differs')
    prepared=[]
    for u in units:
        ids=set(u['coordinate_order_00_01_10_11']);check(len(ids)==4,'Repeated square corner')
        edges={tuple(sorted((a,b))) for a,b in itertools.combinations(ids,2)
            if sum(x!=y for x,y in zip(rows[a]['sequence'],rows[b]['sequence']))==1}
        check(len(edges)==4 and edges<=all_edges,'Square mutation edges do not match sequences')
        check({rows[v]['seed'] for v in ids}=={u['seed']},'Seed assignment mismatch')
        check({str(Seq(rows[v]['sequence']).translate()) for v in ids}=={u['peptide']},'Peptide assignment mismatch')
        fixed=set(mask['arms'][u['arm']]['masked_IDs']);m=mu[u['square_index']]
        check(set(m['masked_IDs_in_square'])==ids&fixed,'Mask intersections differ')
        check(m['prediction_eligible_from_metadata']==(len(ids&fixed)==1),'Mask eligibility differs')
        prepared.append(dict(u,ids=ids,edges=edges,component_ids=set(components[u['component_index']]['IDs']),
            group=seed_group[u['seed']],masked_ids=ids&fixed,prediction_eligible=m['prediction_eligible_from_metadata']))
    pairs=[];degrees={u['square_index']:{r:0 for r in RELATIONS} for u in prepared}
    for a,b in itertools.combinations(prepared,2):
        flags={
            'shared_square_ID':bool(a['ids']&b['ids']),
            'shared_mutation_edge':bool(a['edges']&b['edges']),
            'same_original_masked_target':a['withheld_ID']==b['withheld_ID'],
            'shared_global_masked_ID_in_square':bool(a['masked_ids']&b['masked_ids']),
            'same_peptide':a['peptide']==b['peptide'],'same_design_seed':a['seed']==b['seed'],
            'same_related_seed_group':a['group']==b['group'],
            'same_observed_component':a['component_index']==b['component_index'],
            'shared_component_ID':bool(a['component_ids']&b['component_ids'])}
        for r,v in flags.items():
            degrees[a['square_index']][r]+=int(v);degrees[b['square_index']][r]+=int(v)
        pairs.append({'square_a':a['square_index'],'square_b':b['square_index'],
            'arm_a':a['arm'],'arm_b':b['arm'],'cross_arm':a['arm']!=b['arm'],**flags,
            'shared_IDs':'|'.join(sorted(a['ids']&b['ids'])),
            'shared_edges':'|'.join(':'.join(e) for e in sorted(a['edges']&b['edges'])),
            'shared_masked_IDs':'|'.join(sorted(a['masked_ids']&b['masked_ids']))})
    cross={r:sum(p['cross_arm'] and p[r] for p in pairs) for r in RELATIONS}
    # Any unexpected design/sequence overlap is a blocking pre-outcome result.
    check(not any(cross.values()),'Unexpected cross-arm overlap: '+json.dumps(cross))
    summaries={};group_rows=[];target_rows=[];square_rows=[]
    for arm in ('development','reserved_validation'):
        subset=[u for u in prepared if u['arm']==arm];groups=sorted({u['group'] for u in subset},key=lambda x:tuple(map(int,x.split(','))))
        localpairs=[p for p in pairs if p['arm_a']==arm and p['arm_b']==arm]
        pred=[u for u in subset if u['prediction_eligible']]
        target_count=collections.Counter(u['withheld_ID'] for u in pred)
        for target,n in sorted(target_count.items()):
            selected=[u for u in pred if u['withheld_ID']==target]
            check(len({u['group'] for u in selected})==1,'One target crosses groups')
            target_rows.append({'arm':arm,'target_ID':target,'group':selected[0]['group'],'prediction_cases':n,
                'square_indices':'|'.join(map(str,sorted(u['square_index'] for u in selected)))})
        for g in groups:
            part=[u for u in subset if u['group']==g];pp=[u for u in part if u['prediction_eligible']]
            group_rows.append({'arm':arm,'group':g,'selected_squares':len(part),'selected_diagonals':2*len(part),
                'metadata_prediction_squares':len(pp),'distinct_prediction_targets':len({u['withheld_ID'] for u in pp}),
                'square_vertices':len(set.union(*(u['ids'] for u in part))),
                'component_vertices':len(set.union(*(u['component_ids'] for u in part))),
                'peptides':len({u['peptide'] for u in part}),'observed_components':len({u['component_index'] for u in part}),
                'isolated_observed_squares':sum(u['isolated_observed_square'] for u in part),
                'raw_square_weight':len(part)/len(subset),'primary_equal_group_weight':1/len(groups),
                'prediction_raw_case_weight':len(pp)/len(pred),
                'maximum_within_group_shared_ID_square_degree':max(degrees[u['square_index']]['shared_square_ID'] for u in part)})
        summaries[arm]={'squares':len(subset),'diagonals':2*len(subset),'groups':len(groups),
            'square_vertices':len(set.union(*(u['ids'] for u in subset))),
            'component_vertices':len(set.union(*(u['component_ids'] for u in subset))),
            'square_mutation_edges':len(set.union(*(u['edges'] for u in subset))),
            'prediction_cases':len(pred),'distinct_prediction_targets':len(target_count),
            'duplicate_target_cases':sum(n-1 for n in target_count.values()),
            'relation_pairs':{r:sum(p[r] for p in localpairs) for r in RELATIONS},
            'square_relation_max_degree':{r:max(degrees[u['square_index']][r] for u in subset) for r in RELATIONS},
            'square_count_concentration':stats([sum(u['group']==g for u in subset) for g in groups]),
            'prediction_case_concentration':stats([sum(u['group']==g for u in pred) for g in groups]),
            'distinct_target_concentration':stats([len({u['withheld_ID'] for u in pred if u['group']==g}) for g in groups]),
            'equal_group_weight':1/len(groups)}
        for u in subset:
            square_rows.append({'square_index':u['square_index'],'arm':arm,'seed':u['seed'],'group':u['group'],
                'component_index':u['component_index'],'peptide':u['peptide'],
                'masked_target':u['withheld_ID'],'prediction_eligible':u['prediction_eligible'],
                'IDs':'|'.join(sorted(u['ids'])),**{r+'_degree':degrees[u['square_index']][r] for r in RELATIONS}})
    result={'status':'metadata_overlap_audit_passed','scope':'No phenotype outcomes or archive accessed.',
        'input_hashes':EXPECTED,'source_sha256':sha(Path(__file__)),'original_squares':47,'original_diagonals':94,
        'all_square_pairs_inspected':len(pairs),'cross_arm_square_pairs':sum(p['cross_arm'] for p in pairs),
        'cross_arm_overlap_pair_counts':cross,'arms':summaries,
        'experimental_dependence':{'common_reporter_backbone':True,'common_author_study_and_library':True,
            'protein_preparations':4,'growth_preparations':3,'shared_protein_normalization':True,
            'growth_shared_initial_culture_and_reference':True,'protein_growth_preparations_paired':False,
            'source_independent_assay_split':False,'cross_arm_shared_assay_scope':'All analyzed units',
            'quantitative_assay_covariance':'Unavailable; not estimated from sequence metadata'},
        'interpretation':'Zero recorded design overlap does not imply statistical or experimental independence.'}
    args.output.mkdir(parents=True,exist_ok=False)
    write_csv(args.output/'square_pair_relations.csv',pairs)
    write_csv(args.output/'square_degrees.csv',square_rows)
    write_csv(args.output/'group_inventory.csv',group_rows)
    write_csv(args.output/'prediction_target_multiplicity.csv',target_rows)
    write_json(args.output/'audit.json',result)
    write_json(args.output/'receipt.json',{'status':result['status'],'input_hashes':EXPECTED,'source_sha256':sha(Path(__file__)),
        'files':{p.name:sha(p) for p in args.output.iterdir() if p.is_file()},'outcome_values_read':0})
    print(json.dumps({'status':result['status'],'cross_arm':cross,'arms':summaries},indent=2))

if __name__=='__main__':main()
