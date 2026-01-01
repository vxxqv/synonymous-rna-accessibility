"""Synthetic overlap leakage tests, plus metadata invariants only."""
import hashlib
import json
from pathlib import Path
import unittest
from masking_v2 import union_mask_manifest,predict_from_visible_ids

ROOT=Path(__file__).resolve().parent


class PoisonValues(dict):
    def __getitem__(self,key):
        raise AssertionError('A rejected input accessed a value')
    def values(self):
        raise AssertionError('A rejected input accessed values')


def toy_unit(index,ids,target,seed='2',arm='development'):
    return {'square_index':index,'coordinate_order_00_01_10_11':ids,'withheld_ID':target,
        'seed':seed,'arm':arm,'diagonals':[[ids[0],ids[3]],[ids[1],ids[2]]]}


class MaskTests(unittest.TestCase):
    def test_overlapping_square_target_cannot_be_known_elsewhere(self):
        units=[toy_unit(0,['a','b','c','d'],'d'),toy_unit(1,['d','e','f','g'],'e')]
        result=union_mask_manifest({'units':units,'seed_groups':[['2']]})
        self.assertEqual(result['arms']['development']['masked_IDs'],['d','e'])
        self.assertTrue(result['units'][0]['prediction_eligible_from_metadata'])
        self.assertFalse(result['units'][1]['prediction_eligible_from_metadata'])
        self.assertEqual(result['units'][1]['prediction_exclusion_reason'],'multiple_globally_masked_corners')
        self.assertEqual(result['descriptive_squares_retained'],2)
        self.assertEqual(result['descriptive_diagonals_retained'],4)
        #Even the excluded square's target remains globally masked.
        self.assertIn('e',result['arms']['development']['masked_IDs'])
        with self.assertRaises(ValueError):
            predict_from_visible_ids(['d','e','f','g'],{'d','e'},PoisonValues(d=1,f=2,g=3))

    def test_valid_shared_target_is_deduplicated_not_excluded(self):
        units=[toy_unit(0,['a','b','c','d'],'d'),toy_unit(1,['d','e','f','g'],'d')]
        result=union_mask_manifest({'units':units,'seed_groups':[['2']]})
        self.assertEqual(result['arms']['development']['unique_masked_IDs'],1)
        self.assertEqual(result['arms']['development']['prediction_squares_metadata_eligible'],2)

    def test_strict_three_known_ids_without_any_extra_table(self):
        answer=predict_from_visible_ids(['a','b','c','d'],{'d','other_target'},{'a':10,'b':1,'c':1})
        self.assertEqual(answer['predictions']['additive'],-8)
        for keys in ({'a','b','c','d'},{'a','b','c','extra'},{'a','b','other_target'},{'a','b'}):
            with self.assertRaises(ValueError):
                predict_from_visible_ids(['a','b','c','d'],{'d','other_target'},PoisonValues.fromkeys(keys,0))
        for mask in (set(),{'c','d'}):
            with self.assertRaises(ValueError):
                predict_from_visible_ids(['a','b','c','d'],mask,PoisonValues(a=1,b=2,c=3))

    def test_cross_arm_identity_is_refused(self):
        units=[toy_unit(0,['a','b','c','d'],'d'),toy_unit(1,['d','e','f','g'],'d',seed='3',arm='reserved_validation')]
        with self.assertRaises(ValueError):union_mask_manifest({'units':units,'seed_groups':[['2'],['3']]})

    def test_real_metadata_only_exact_set_invariants(self):
        old=json.loads((ROOT/'split_manifest.json').read_text())
        new=json.loads((ROOT/'prediction_mask_v2.json').read_text())
        self.assertEqual(len(new['units']),47)
        self.assertEqual(new['descriptive_diagonals_retained'],94)
        for arm in ('development','reserved_validation'):
            raw=[u for u in old['units'] if u['arm']==arm]
            fixed_targets={u['withheld_ID'] for u in raw}
            self.assertEqual(set(new['arms'][arm]['masked_IDs']),fixed_targets)
            for unit in (u for u in new['units'] if u['arm']==arm):
                hits=set(unit['coordinate_order_00_01_10_11'])&fixed_targets
                self.assertEqual(unit['prediction_eligible_from_metadata'],len(hits)==1)
                if unit['prediction_eligible_from_metadata']:
                    self.assertEqual(set(unit['visible_prediction_IDs']),set(unit['coordinate_order_00_01_10_11'])-fixed_targets)
                    self.assertFalse(set(unit['visible_prediction_IDs'])&fixed_targets)
        original=(ROOT/'preflight_manifest.json').read_bytes()
        self.assertEqual(original,(ROOT/'preflight_manifest_v1.json').read_bytes())
        for entry in json.loads(original)['files']:
            self.assertEqual(hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest(),entry['sha256'])


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MaskTests))
    receipt={'status':'passed' if result.wasSuccessful() else 'failed','tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
        'outcome_scope':'Synthetic values and existing sequence metadata only','original_v1_artifacts_preserved':True,
        'hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('masking_v2.py','test_masking_v2.py','prediction_mask_v2.json','preflight_manifest_v1.json')}}
    with (ROOT/'synthetic_mask_receipt_v2.json').open('x',encoding='utf-8') as stream:stream.write(json.dumps(receipt,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
