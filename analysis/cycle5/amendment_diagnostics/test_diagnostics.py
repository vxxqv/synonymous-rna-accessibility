"""Synthetic-only tests. No author functions, archives or actual evaluation files."""
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import diagnostics as d


def classified(label):
    return {'classification': label, 'thresholds': [0, 0],
            'individual_connectivity': [True, True], 'joint_connectivity': label == 'jointly_connected'}


def bound(connected=True, unknown=False):
    return {'status': 'inconclusive' if unknown else ('certified_jointly_connected' if connected else 'certified_joint_only_obstruction'),
            'thresholds': [0, 0], 'individual_pessimistic': [True, True],
            'individual_optimistic': [True, True], 'joint_pessimistic': connected and not unknown,
            'joint_optimistic': connected or unknown}


def field(connected=True, unknown=False):
    return {'complete': not unknown, 'point': None if unknown else classified('jointly_connected' if connected else 'joint_only_obstruction'),
            'bounds': bound(connected, unknown)}


def fixtures():
    units = []; descriptions = []; scores = []; coverage = []
    specs = [(1, i) for i in range(3)] + [(i, 0) for i in range(2, 17)]
    for index, (group, within) in enumerate(specs):
        target = 'SYN_shared_target' if group == 1 and within < 2 else f'SYN_{index}_target'
        ids = [target] + [f'SYN_{index}_{i}' for i in range(1, 4)]
        unit = {'square_index': index, 'seed': str(group), 'arm': 'development',
                'coordinate_order_00_01_10_11': ids, 'diagonals': [[ids[0], ids[3]], [ids[1], ids[2]]],
                'withheld_ID': target}
        units.append(unit)
        eligible = group != 16
        desc = {'square_index': index, 'seed': str(group), 'group': str(group), 'eligible': eligible, 'diagonals': []}
        coverage.append({'square_index': index, 'seed': str(group), 'group': str(group), 'eligible': eligible})
        for j, endpoints in enumerate(unit['diagonals'] if eligible else []):
            local_conflict = group == 1 and (j == 0 or within != 1)
            local = classified('joint_only_obstruction' if local_conflict else 'jointly_connected')
            component = field(connected=(index % 2 == 0), unknown=(index % 3 == 0))
            alt = classified('jointly_connected') if index == 0 and j == 0 else copy.deepcopy(local)
            scenarios = []
            for protein in d.PROTEIN:
                for growth in d.GROWTH:
                    scenarios.append({'protein_column': protein, 'growth_column': growth,
                        'local': copy.deepcopy(local), 'component': copy.deepcopy(component),
                        'tolerance': [1e-10, 1e-10], 'local_alternative_tolerance': copy.deepcopy(alt),
                        'component_alternative_tolerance': copy.deepcopy(component)})
            cert = bound(connected=not (group == 1 and within == 0))
            desc['diagonals'].append({'endpoint_IDs': endpoints, 'primary': local,
                'point_tolerance': [1e-10, 1e-10], 'local_alternative_tolerance': alt,
                'observed_component': component, 'observed_component_alternative_tolerance': copy.deepcopy(component),
                'local_range_certificate': cert, 'component_range_certificate': copy.deepcopy(cert),
                'interval_tolerance': [1e-10, 1e-10], 'local_range_alternative_tolerance': copy.deepcopy(cert),
                'component_range_alternative_tolerance': copy.deepcopy(cert), 'native_shortfalls': {},
                'preparation_scenarios': scenarios, 'preparation_consistent_joint_obstruction': local_conflict})
        descriptions.append(desc)
        score = {'square_index': index, 'seed': str(group), 'group': str(group), 'status': 'scored' if eligible else 'unscorable_missing_target_or_square'}
        if eligible:
            observed = [10.0, 20.0]
            truth = group == 1 and within != 1
            error = [1, 3, 9][within] if group == 1 else 2
            score.update(observed_target=observed, true_joint_only_obstruction=truth, methods={})
            for method, multiplier in zip(d.METHODS, (1, 2, 4)):
                errors = [error*multiplier, error*multiplier*2]
                pred = truth if method == 'additive' else False
                score['methods'][method] = {'predicted_target': [y+e for y, e in zip(observed, errors)],
                    'signed_error': errors, 'absolute_error': errors, 'squared_error': [e*e for e in errors],
                    'predicted_joint_only_obstruction': pred, 'classification_error': int(pred != truth)}
        scores.append(score)
    split = {'seed_groups': [[str(i)] for i in range(1, 18)],
             'units': units + [{'square_index': 999, 'seed': '17', 'arm': 'reserved_validation',
                               'coordinate_order_00_01_10_11': ['FOREIGN_'+str(i) for i in range(4)],
                               'withheld_ID': 'FOREIGN_0', 'diagonals': [['FOREIGN_0', 'FOREIGN_3'], ['FOREIGN_1', 'FOREIGN_2']]}]}
    evaluation = {'arm': 'development', 'status': 'arm_evaluation_complete', 'interpretation_permitted': True,
                  'coverage': coverage, 'descriptions': descriptions, 'prediction_evaluation': {'records': scores},
                  'source_hashes': {'runner': d.RUNNER_SHA256}}
    return split, evaluation


class TestSyntheticDiagnostics(unittest.TestCase):
    def setUp(self):
        self.split, self.evaluation = fixtures()

    def analyze(self):
        return d.analyze(self.evaluation, self.split, 'development')

    def test_hand_calculated_losses_and_target_balance(self):
        result = self.analyze()['full_arm']
        case = result['case_balanced_prediction']; target = result['target_balanced_prediction']
        self.assertEqual((case['scored_cases'], case['contributing_groups'], case['distinct_targets']), (17, 15, 16))
        self.assertAlmostEqual(case['summaries']['additive']['continuous']['protein']['MAE'], 97/45)
        self.assertAlmostEqual(target['summaries']['additive']['continuous']['protein']['MAE'], 67/30)
        self.assertAlmostEqual(case['summaries']['additive']['continuous']['protein']['RMSE'], math.sqrt(259/45))
        self.assertAlmostEqual(target['summaries']['additive']['continuous']['protein']['RMSE'], math.sqrt(6.6))
        self.assertTrue(case['discrimination_screen_passed'])

    def test_zero_group_explicit_and_equal_group_rate(self):
        result = self.analyze()['full_arm']['descriptive']
        self.assertEqual(result['raw_local']['eligible_diagonals'], 34)
        self.assertEqual(result['raw_local']['joint_only'], 5)
        self.assertAlmostEqual(result['raw_local']['fractions']['joint_only'], 5/34)
        self.assertAlmostEqual(result['equal_contributing_group_mean_fractions']['joint_only'], 1/18)
        self.assertEqual(result['zero_eligible_groups'], ['16'])
        self.assertIsNone(result['group_denominators']['16']['fractions']['joint_only'])
        self.assertEqual(result['selected_groups'], 16)

    def test_all_numeric_leave_outs_and_ranges(self):
        result = self.analyze(); loo = result['leave_one_group_out']
        self.assertEqual([r['omitted_group'] for r in loo], list(map(str, range(1, 17))))
        self.assertEqual(loo[0]['case_balanced_prediction']['scored_cases'], 14)
        self.assertEqual(loo[0]['case_balanced_prediction']['summaries']['additive']['continuous']['protein']['MAE'], 2)
        self.assertFalse(loo[0]['case_balanced_prediction']['both_true_classes'])
        self.assertEqual(loo[-1]['case_balanced_prediction']['scored_cases'], 17)
        key = 'case_balanced_prediction/summaries/additive/continuous/protein/MAE'
        self.assertEqual(result['all_leave_one_out_minmax_ranges'][key]['all_rows'], 16)

    def test_no_ensemble_of_endpoint_labels(self):
        result = self.analyze()
        labels = [r for r in result['prediction_case_labels_preserved'] if r['withheld_ID'] == 'SYN_shared_target']
        self.assertEqual([r['true_label'] for r in labels], [True, False])
        self.assertNotEqual(labels[0]['endpoint_IDs'], labels[1]['endpoint_IDs'])

    def test_all_scenarios_and_tolerance_flips(self):
        result = self.analyze()['observed_component_diagnostics']
        self.assertEqual(result['eligible_diagonals'], 34)
        self.assertTrue(all(len(r['preparation_scenarios']) == 12 for r in result['records']))
        self.assertEqual(result['flip_counts']['local_point'], 1)
        self.assertEqual(result['flip_counts']['scenario_local'], 12)
        self.assertIn('unresolved', result['records'][0]['bypass_status_tolerance_zero'])

    def test_reviewer_component_count_fixture(self):
        expected_path = Path(__file__).parent/'reviewer_expected_component_count_fixture.json'
        self.assertEqual(hashlib.sha256(expected_path.read_bytes()).hexdigest(),
                         '19e775e320fa7c95bbb4b1a8822694ab59260f0dfa9ec35ee1eab420efb62eab')
        expected = json.loads(expected_path.read_text())
        tables = self.analyze()['observed_component_diagnostics']['count_tables']
        zero = tables['primary_point_field']['zero_tolerance']
        alt = tables['primary_point_field']['alternative_tolerance']
        def agrees(actual, wanted, denominator=34):
            self.assertEqual(actual['denominator'], denominator)
            self.assertEqual(sum(actual['counts'].values()), denominator)
            self.assertEqual({k: v for k, v in actual['counts'].items() if v}, wanted)
        agrees(zero['component_point_classifications'], expected['counts']['point_component_classification'])
        agrees(zero['component_bound_statuses'], expected['counts']['point_component_bounds'])
        agrees(alt['component_bound_statuses'], expected['counts']['point_component_bounds'])
        agrees(zero['bypass_outcomes'], expected['counts']['point_bypass'])
        agrees(alt['bypass_outcomes'], expected['counts']['alternative_point_bypass'])
        self.assertEqual(zero['bypass_outcomes']['conditional_local_joint_only_denominator'], 5)
        self.assertEqual(alt['bypass_outcomes']['conditional_local_joint_only_denominator'], 4)
        self.assertEqual((zero['available_component_points'], zero['unavailable_component_points']), (22, 12))
        self.assertEqual(len(tables['named_preparation_scenarios']), 12)
        self.assertFalse(tables['pooled_across_scenarios'])
        for key, wanted in expected['per_scenario'].items():
            row = tables['named_preparation_scenarios'][key]
            agrees(row['zero_tolerance']['bypass_outcomes'], wanted['bypass'])
            agrees(row['alternative_tolerance']['bypass_outcomes'], wanted['alternative_bypass'])
            agrees(row['zero_tolerance']['component_bound_statuses'], wanted['component_bounds'])
            agrees(row['alternative_tolerance']['component_bound_statuses'], wanted['component_bounds'])
            agrees(row['zero_tolerance']['component_point_classifications'], expected['counts']['point_component_classification'])
        for tolerance in ('zero_tolerance', 'alternative_tolerance'):
            row = tables['empirical_range_field'][tolerance]
            agrees(row['local_certificate_statuses'], expected['counts']['range_local_bounds'])
            agrees(row['component_certificate_statuses'], expected['counts']['range_component_bounds'])
            self.assertEqual(row['component_certificate_statuses']['counts']['inconclusive'], 0)

    def test_empty_component_count_tables(self):
        tables = d.aggregate_component_counts([])
        self.assertEqual(len(tables['named_preparation_scenarios']), 12)
        points = list(tables['primary_point_field'].values())
        points += [v for scenario in tables['named_preparation_scenarios'].values()
                   for k, v in scenario.items() if k in ('zero_tolerance', 'alternative_tolerance')]
        for table in points:
            self.assertEqual(table['eligible_diagonals'], 0)
            for key in ('local_classifications', 'component_point_classifications', 'component_bound_statuses', 'bypass_outcomes'):
                self.assertEqual(table[key]['denominator'], 0)
                self.assertTrue(all(n == 0 for n in table[key]['counts'].values()))
            self.assertTrue(all(v is None for v in table['bypass_outcomes']['conditional_fractions'].values()))
        for table in tables['empirical_range_field'].values():
            self.assertEqual(table['eligible_diagonals'], 0)
            self.assertEqual(sum(table['component_certificate_statuses']['counts'].values()), 0)
        for row in self.evaluation['descriptions']:
            row.update(eligible=False, diagonals=[])
        for row in self.evaluation['coverage']:
            row['eligible'] = False
        for row in self.evaluation['prediction_evaluation']['records']:
            row['status'] = 'missing_known_replicates'
        result = self.analyze()
        self.assertEqual(result['observed_component_diagnostics']['count_tables'], tables)

    def test_unavailable_component_point_and_bounds_are_not_zeros_or_failures(self):
        local = classified('joint_only_obstruction')
        row = d.component_count_table([(local, {'complete': False, 'point': None, 'bounds': None}),
                                       (None, None)])
        self.assertEqual(row['eligible_diagonals'], 2)
        self.assertEqual(row['component_point_classifications']['counts']['point_unavailable'], 2)
        self.assertEqual(row['component_bound_statuses']['counts']['certificate_unavailable'], 2)
        self.assertEqual(row['bypass_outcomes']['counts']['observed_component_bypass_unresolved'], 1)
        self.assertEqual(row['bypass_outcomes']['counts']['local_classification_unavailable'], 1)
        self.assertEqual(row['bypass_outcomes']['counts']['not_a_local_joint_only_obstruction'], 0)
        self.assertEqual(row['bypass_outcomes']['conditional_local_joint_only_denominator'], 1)
        records = self.analyze()['observed_component_diagnostics']['records']
        for record in records:
            for key in ('local_range_certificate', 'component_range_certificate',
                        'local_range_alternative_tolerance', 'component_range_alternative_tolerance'):
                record[key] = None
        tables = d.aggregate_component_counts(records)
        for table in tables['empirical_range_field'].values():
            self.assertEqual(table['local_certificate_statuses']['counts']['certificate_unavailable'], 34)
            self.assertEqual(table['component_certificate_statuses']['counts']['certificate_unavailable'], 34)
            self.assertEqual(table['component_certificate_statuses']['counts']['certified_joint_only_obstruction'], 0)

    def test_v1_numeric_and_record_outputs_remain_identical(self):
        path = Path(__file__).parent/'history_v1/diagnostics.py'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         '9e0307ad72b6dd0be771dee9d886bb87e63d09ed85a7ceea5409b2e3f90f897d')
        spec = importlib.util.spec_from_file_location('preserved_diagnostics_v1', path)
        old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
        before = old.analyze(self.evaluation, self.split, 'development'); after = self.analyze()
        for key in before:
            if key != 'observed_component_diagnostics':
                self.assertEqual(before[key], after[key])
        for key in before['observed_component_diagnostics']:
            self.assertEqual(before['observed_component_diagnostics'][key], after['observed_component_diagnostics'][key])

    def test_zero_baseline_cannot_pass(self):
        for row in self.evaluation['prediction_evaluation']['records']:
            if row['status'] == 'scored':
                row['methods']['adjacent_mean'].update(predicted_target=row['observed_target'][:], signed_error=[0, 0], absolute_error=[0, 0], squared_error=[0, 0])
        result = self.analyze()['full_arm']['case_balanced_prediction']
        self.assertFalse(result['discrimination_screen_passed'])
        self.assertIsNone(result['MAE_comparisons'][0]['relative_gain'])

    def test_no_scored_cases_or_coverage(self):
        for row in self.evaluation['prediction_evaluation']['records']:
            row['status'] = 'missing_known_replicates'; row.pop('methods', None)
        result = self.analyze()['full_arm']['case_balanced_prediction']
        self.assertEqual(result['scored_cases'], 0)
        self.assertFalse(result['discrimination_screen_passed'])
        self.assertTrue(all(v['cases'] == 0 for v in result['groups'].values()))

    def test_one_class_is_not_discrimination_success(self):
        for row in self.evaluation['descriptions']:
            for diagonal in row['diagonals']:
                diagonal['primary'] = classified('jointly_connected')
                diagonal['preparation_consistent_joint_obstruction'] = False
                for scenario in diagonal['preparation_scenarios']:
                    scenario['local'] = classified('jointly_connected')
        for row in self.evaluation['prediction_evaluation']['records']:
            if row['status'] == 'scored':
                row['true_joint_only_obstruction'] = False
                for method in d.METHODS:
                    row['methods'][method]['predicted_joint_only_obstruction'] = False
                    row['methods'][method]['classification_error'] = 0
        result = self.analyze()['full_arm']['case_balanced_prediction']
        self.assertFalse(result['both_true_classes'])
        self.assertFalse(result['discrimination_screen_passed'])
        self.assertIsNone(result['summaries']['additive']['pooled_case_balanced_accuracy'])

    def test_nine_cases_are_a_coverage_failure(self):
        for row in self.evaluation['prediction_evaluation']['records'][9:]:
            row['status'] = 'missing_known_replicates'; row.pop('methods', None)
        result = self.analyze()['full_arm']['case_balanced_prediction']
        self.assertEqual(result['scored_cases'], 9)
        self.assertIn('prediction_coverage_below_10_squares_8_groups', result['failure_reasons'])

    def test_cross_arm_injection_and_incomplete_inventory(self):
        e = copy.deepcopy(self.evaluation); e['coverage'][0]['square_index'] = 999
        with self.assertRaises(ValueError): d.analyze(e, self.split, 'development')
        self.evaluation['coverage'].pop()
        with self.assertRaises(ValueError): self.analyze()

    def test_loss_corruption_shared_target_change_and_scenario_loss(self):
        self.evaluation['prediction_evaluation']['records'][0]['methods']['additive']['squared_error'][0] += 1
        with self.assertRaises(ValueError): self.analyze()
        self.split, self.evaluation = fixtures()
        self.evaluation['prediction_evaluation']['records'][1]['observed_target'][0] += 1
        with self.assertRaises(ValueError): self.analyze()
        self.split, self.evaluation = fixtures()
        self.evaluation['descriptions'][0]['diagonals'][0]['preparation_scenarios'].pop()
        with self.assertRaises(ValueError): self.analyze()

    def test_cli_guards_hold_and_input_immutability(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            root = Path(tmp); sp = root/'split.json'; ep = root/'evaluation.json'
            sp.write_text(json.dumps(self.split)); self.evaluation['source_hashes']['split'] = d.digest(sp.read_bytes())
            ep.write_text(json.dumps(self.evaluation)); h = d.digest(ep.read_bytes())
            before = (ep.read_bytes(), sp.read_bytes())
            with self.assertRaises(ValueError): d.run(ep, '0'*64, sp, 'development', root/'bad_hash')
            with self.assertRaises(ValueError): d.run(ep, h, sp, 'reserved_validation', root/'bad_arm')
            self.assertFalse((root/'bad_arm').exists())
            self.assertEqual(d.run(ep, h, sp, 'development', root/'passed'), 'supplemental_diagnostics_complete')
            with self.assertRaises(FileExistsError): d.run(ep, h, sp, 'development', root/'passed')
            self.assertEqual(before, (ep.read_bytes(), sp.read_bytes()))
            self.evaluation.update(status='held_for_author_mean_provenance', interpretation_permitted=False, holds=['SYN_HOLD'])
            self.evaluation.pop('descriptions'); self.evaluation.pop('prediction_evaluation')
            ep.write_text(json.dumps(self.evaluation)); h = d.digest(ep.read_bytes())
            self.assertEqual(d.run(ep, h, sp, 'development', root/'held'), 'provenance_hold_no_interpretation')
            held = json.loads((root/'held/diagnostics.json').read_text())
            self.assertNotIn('full_arm', held)

    def test_optimization_refusal_before_access(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            output = Path(tmp)/'no_output'
            for flag in ('-O', '-OO'):
                result = subprocess.run([sys.executable, flag, '-B', str(Path(d.__file__).resolve()),
                    '--evaluation', 'NOT_EXIST', '--evaluation-sha256', '0'*64, '--split', 'NOT_EXIST',
                    '--allow-arm', 'development', '--output', str(output)], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Refuse -O and -OO', result.stderr)
                self.assertFalse(output.exists())

    def test_hash_schema_and_duplicate_json_refusal(self):
        with self.assertRaises(ValueError): d.load('{"x":1,"x":2}')
        with self.assertRaises(ValueError): d.load('{"x":NaN}')
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            root = Path(tmp); sp = root/'split.json'; ep = root/'eval.json'
            sp.write_text(json.dumps(self.split)); self.evaluation['source_hashes']['split'] = '0'*64
            ep.write_text(json.dumps(self.evaluation))
            with self.assertRaises(ValueError): d.run(ep, d.digest(ep.read_bytes()), sp, 'development', root/'bad_split')
            self.evaluation['source_hashes']['split'] = d.digest(sp.read_bytes())
            self.evaluation['source_hashes']['runner'] = '0'*64; ep.write_text(json.dumps(self.evaluation))
            with self.assertRaises(ValueError): d.run(ep, d.digest(ep.read_bytes()), sp, 'development', root/'bad_runner')


if __name__ == '__main__':
    unittest.main(verbosity=2)
