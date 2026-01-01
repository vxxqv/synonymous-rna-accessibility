"""Metadata-only union masking and a strict three-visible-corner prediction API."""
import hashlib
import json
from pathlib import Path

from graph_rules import missing_corner_predictions

ROOT=Path(__file__).resolve().parent
V1_SPLIT_SHA='394a0c1aa23785c290b31c147d199f406cf46b269daa9159a98c419398b05a4e'


def union_mask_manifest(split):
    """Fix each arm's mask from ALL original squares before any exclusion."""
    units=split['units']
    masks={arm:sorted({u['withheld_ID'] for u in units if u['arm']==arm}) for arm in ('development','reserved_validation')}
    if set(masks['development']) & set(masks['reserved_validation']):
        raise ValueError('Cross-arm masked identity leakage')
    annotated=[]
    for unit in units:
        ids=unit['coordinate_order_00_01_10_11']
        if len(set(ids))!=4 or unit['withheld_ID'] not in ids:
            raise ValueError('Invalid original square')
        masked=sorted(set(ids)&set(masks[unit['arm']]))
        if not masked or unit['withheld_ID'] not in masked:
            raise AssertionError('Original target absent from union mask')
        eligible=len(masked)==1
        annotated.append({'square_index':unit['square_index'],'seed':unit['seed'],'arm':unit['arm'],
            'descriptive_retained':True,'descriptive_diagonals':unit['diagonals'],
            'coordinate_order_00_01_10_11':ids,'original_withheld_ID':unit['withheld_ID'],
            'masked_IDs_in_square':masked,'mask_count':len(masked),
            'prediction_eligible_from_metadata':eligible,
            'prediction_exclusion_reason':None if eligible else 'multiple_globally_masked_corners',
            'visible_prediction_IDs':sorted(set(ids)-set(masked)) if eligible else [],
            'prediction_measurement_coverage':'unexamined'})
    arms={}
    for arm,mask in masks.items():
        selected=[u for u in annotated if u['arm']==arm]
        eligible=[u for u in selected if u['prediction_eligible_from_metadata']]
        groups=[]
        for group in split['seed_groups']:
            group_units=[u for u in selected if u['seed'] in group]
            if group_units:
                groups.append({'seed_group':group,'descriptive_squares':len(group_units),
                    'prediction_squares':sum(u['prediction_eligible_from_metadata'] for u in group_units),
                    'prediction_excluded_squares':[u['square_index'] for u in group_units if not u['prediction_eligible_from_metadata']]})
        usable_groups=sum(g['prediction_squares']>0 for g in groups)
        arms[arm]={'masked_IDs':mask,'unique_masked_IDs':len(mask),
            'descriptive_squares':len(selected),'descriptive_diagonals':sum(len(u['descriptive_diagonals']) for u in selected),
            'prediction_squares_metadata_eligible':len(eligible),'prediction_squares_metadata_excluded':len(selected)-len(eligible),
            'prediction_related_seed_groups':usable_groups,'counts_by_related_seed_group':groups,
            'metadata_upper_bound_passes_10_square_8_group_minimum':len(eligible)>=10 and usable_groups>=8,
            'phenotype_coverage_and_both_class_requirement':'unexamined'}
    return {'schema':'Cambray_shared_corner_mask_v2','status':'metadata_union_mask_fixed',
        'mask_formation':'Union of every original per-square hashed target within each arm, before exclusions; never recomputed from eligible subset',
        'descriptive_squares_retained':len(annotated),'descriptive_diagonals_retained':sum(len(u['descriptive_diagonals']) for u in annotated),
        'arms':arms,'units':annotated,'no_phenotype_values_read':True}


def predict_from_visible_ids(coordinate_ids, global_mask, known_values):
    """No table access or target values are accepted.

    known_values must have exactly the three visible IDs of this square. All key
    checks precede any value access. Caller must enforce the hash-pinned arm mask
    at input extraction, freeze predictions, then separately authorize outcomes.
    """
    if len(coordinate_ids)!=4 or len(set(coordinate_ids))!=4:
        raise ValueError('Require four distinct coordinate IDs')
    masked=set(coordinate_ids)&set(global_mask)
    if len(masked)!=1:
        raise ValueError('Prediction square must contain exactly one globally masked ID')
    expected=set(coordinate_ids)-masked
    if set(known_values)!=expected or set(known_values)&set(global_mask):
        raise ValueError('Accept exactly three unmasked IDs and no other table entries')
    missing_id=next(iter(masked))
    missing=coordinate_ids.index(missing_id)
    numerical={i:known_values[identifier] for i,identifier in enumerate(coordinate_ids) if identifier in expected}
    return {'withheld_ID':missing_id,'withheld_coordinate':missing,'predictions':missing_corner_predictions(numerical,missing)}


def main():
    original=(ROOT/'preflight_manifest.json').read_bytes()
    v1=json.loads(original)
    for entry in v1['files']:
        assert hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()==entry['sha256'], entry['path']
    named=ROOT/'preflight_manifest_v1.json'
    if named.exists():
        assert named.read_bytes()==original
    else:
        named.write_bytes(original)
    payload=(ROOT/'split_manifest.json').read_bytes()
    assert hashlib.sha256(payload).hexdigest()==V1_SPLIT_SHA
    result=union_mask_manifest(json.loads(payload))
    assert result['descriptive_squares_retained']==47 and result['descriptive_diagonals_retained']==94
    result['original_split_sha256']=V1_SPLIT_SHA
    result['original_preflight_manifest_sha256']=hashlib.sha256(original).hexdigest()
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with (ROOT/'prediction_mask_v2.json').open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'descriptive_squares':47,'descriptive_diagonals':94,
        'arms':{arm:{k:v for k,v in info.items() if k not in ('masked_IDs','counts_by_related_seed_group')} for arm,info in result['arms'].items()}},indent=2))

if __name__=='__main__':main()
