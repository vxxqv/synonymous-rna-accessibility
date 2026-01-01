"""Declared-design membership only for the fixed certificate vertices.

No source2016 count/phenotype/coverage file is read. Index base is zero,
matching the frozen certificate's enumeration.
"""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[2]
SOURCE=WORKSPACE/'outputs'/'synonymous-rna-accessibility'/'analysis'/'cycle6'/'fas_pilot'/'source_contract.json'
EXPECTED='f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38'
raw=SOURCE.read_bytes()
if hashlib.sha256(raw).hexdigest()!=EXPECTED:
    raise RuntimeError('Frozen sequence contract mismatch')
cfg=json.loads(raw);wt='GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
ids=[13,16,35,39,42,46,51,70]
rows=[]
for vertex in ids:
    s=cfg['sequences'][vertex]
    mutations=[{'exon_position1':k+1,'WT_base':a,'variant_base':b} for k,(a,b) in enumerate(zip(wt,s)) if a!=b]
    d=len(mutations)
    rows.append({'vertex_zero_based':vertex,'sequence':s,'distance_from_native_WT':d,
        'substitutions_relative_native_WT':mutations,
        'within_declared_2016_analyzed_radius2_design':d<=2,
        'is_native_WT':d==0,'is_single_substitution':d==1,
        'within_documented_individually_validated_single_mutant_category':d==1,
        'observed_2016_measurement_presence':None,
        'direct_2016_inclusion_skip_record_authenticated':False})
receipt={'scope':'Theoretical source-design membership, NOT observed outcome or coverage',
    'index_base':0,'source_contract_sha256':EXPECTED,'fixed_vertices':ids,'rows':rows,
    'declared_radius2_members':sum(z['within_declared_2016_analyzed_radius2_design'] for z in rows),
    'declared_radius2_nonWT_members':sum(0<z['distance_from_native_WT']<=2 for z in rows),
    'outside_declared_radius2':sum(z['distance_from_native_WT']>2 for z in rows),
    'single_mutants_among_nonWT_certificate_vertices':sum(z['is_single_substitution'] for z in rows),
    'observed_coverage_inspected':False,'outcome_tables_or_counts_read':False,
    'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'membership_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:receipt[k] for k in ['declared_radius2_members','declared_radius2_nonWT_members','outside_declared_radius2','single_mutants_among_nonWT_certificate_vertices','observed_coverage_inspected']},indent=2))
