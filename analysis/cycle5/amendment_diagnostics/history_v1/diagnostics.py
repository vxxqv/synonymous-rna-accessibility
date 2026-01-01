"""Post-unseal descriptive supplements. Never reads a raw archive or hidden arm.

This script has not been authorized on real evaluation outputs. No author
numeric functions are imported. Input evaluation bytes and arm are explicit.
"""
if not __debug__:
    raise RuntimeError('Refuse -O and -OO before input access or output creation')

import argparse
import hashlib
import json
import math
from pathlib import Path

RUNNER_SHA256 = '5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5'
METHODS = ('additive', 'adjacent_mean', 'three_corner_mean')
TRAITS = ('protein', 'growth')
PROTEIN = tuple('clean.lin.prot.rep' + v for v in ('One', 'Two', 'Three', 'Four'))
GROWTH = ('growth.R60one', 'growth.R60two', 'growth.R60three')
METRICS = ('joint_only', 'preparation_consistent', 'range_certificate')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON key: ' + key)
        result[key] = value
    return result


def load(data):
    def bad_constant(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(data, object_pairs_hook=unique_object, parse_constant=bad_constant)


def finite(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('Expected a finite numerical value')
    return value


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)


def average(values):
    return math.fsum(values) / len(values) if values else None


def group_order(group):
    return tuple(map(int, group.split(',')))


def metadata(split, arm):
    group_for_seed = {}
    for members in split['seed_groups']:
        if not members or members != sorted(members, key=int):
            raise ValueError('Noncanonical group')
        group = ','.join(members)
        for seed in members:
            if seed in group_for_seed:
                raise ValueError('Repeated seed')
            group_for_seed[seed] = group
    units = {}
    identity_arm = {}; group_arm = {}
    for unit in split['units']:
        if unit['square_index'] in units:
            raise ValueError('Duplicate square metadata')
        units[unit['square_index']] = unit
        group = group_for_seed[unit['seed']]
        if group in group_arm and group_arm[group] != unit['arm']:
            raise ValueError('Related group crosses arms')
        group_arm[group] = unit['arm']
        for identifier in unit['coordinate_order_00_01_10_11']:
            if identifier in identity_arm and identity_arm[identifier] != unit['arm']:
                raise ValueError('Genotype crosses arms')
            identity_arm[identifier] = unit['arm']
    selected = {i: u for i, u in units.items() if u['arm'] == arm}
    groups = sorted({group_for_seed[u['seed']] for u in selected.values()}, key=group_order)
    if len(groups) != 16:
        raise ValueError('Expected exactly16 selected related groups in this frozen design')
    if not selected:
        raise ValueError('Empty arm')
    for unit in selected.values():
        ids = unit['coordinate_order_00_01_10_11']
        if len(ids) != 4 or len(set(ids)) != 4 or unit['withheld_ID'] not in ids:
            raise ValueError('Invalid square or target identity')
    return selected, groups, group_for_seed


def check_identity(record, units, group_for_seed):
    index = record['square_index']
    if index not in units:
        raise ValueError('Foreign-arm or unknown square')
    unit = units[index]
    if record['seed'] != unit['seed'] or record['group'] != group_for_seed[unit['seed']]:
        raise ValueError('Seed/group mismatch')
    return unit


def by_unique_square(records, units, group_for_seed):
    result = {}
    for record in records:
        check_identity(record, units, group_for_seed)
        index = record['square_index']
        if index in result:
            raise ValueError('Duplicate square output')
        result[index] = record
    if set(result) != set(units):
        raise ValueError('Incomplete metadata-selected square inventory')
    return result


def local_summary(descriptions, groups):
    per_group = {g: {'eligible_squares': 0, 'eligible_diagonals': 0,
                     **{m: 0 for m in METRICS}} for g in groups}
    for row in descriptions:
        stats = per_group[row['group']]
        stats['eligible_squares'] += int(row['eligible'])
        for diagonal in row['diagonals']:
            stats['eligible_diagonals'] += 1
            stats['joint_only'] += diagonal['primary']['classification'] == 'joint_only_obstruction'
            stats['preparation_consistent'] += diagonal['preparation_consistent_joint_obstruction']
            stats['range_certificate'] += diagonal['local_range_certificate']['status'] == 'certified_joint_only_obstruction'
    for stats in per_group.values():
        n = stats['eligible_diagonals']
        stats['fractions'] = {m: stats[m] / n if n else None for m in METRICS}
    totals = {key: sum(v[key] for v in per_group.values())
              for key in ('eligible_squares', 'eligible_diagonals', *METRICS)}
    n = totals['eligible_diagonals']
    totals['fractions'] = {m: totals[m] / n if n else None for m in METRICS}
    contributors = [v for v in per_group.values() if v['eligible_diagonals']]
    means = {m: average([v['fractions'][m] for v in contributors]) for m in METRICS}
    return {'raw_local': totals, 'group_denominators': per_group,
            'selected_groups': len(groups), 'contributing_groups': len(contributors),
            'zero_eligible_groups': [g for g in groups if not per_group[g]['eligible_diagonals']],
            'equal_contributing_group_mean_fractions': means,
            'coverage_hurdle_passed': totals['eligible_squares'] >= 10 and len(contributors) >= 8,
            'descriptive_existence_milestone': totals['eligible_squares'] >= 10 and len(contributors) >= 8
                and totals['preparation_consistent'] > 0}


def score_summary(records, groups, target_balanced=False):
    scored = [r for r in records if r['status'] == 'scored']
    per_group = {g: [r for r in scored if r['group'] == g] for g in groups}
    def group_loss(rows, function):
        if not target_balanced:
            return average([function(r) for r in rows])
        targets = sorted({r['withheld_ID'] for r in rows})
        return average([average([function(r) for r in rows if r['withheld_ID'] == target])
                        for target in targets])
    summaries = {}
    for method in METHODS:
        continuous = {}
        for k, trait in enumerate(TRAITS):
            gm = {g: group_loss(rows, lambda r: r['methods'][method]['absolute_error'][k])
                  for g, rows in per_group.items()}
            gs = {g: group_loss(rows, lambda r: r['methods'][method]['squared_error'][k])
                  for g, rows in per_group.items()}
            mse = average([v for v in gs.values() if v is not None])
            continuous[trait] = {'MAE': average([v for v in gm.values() if v is not None]),
                'RMSE': math.sqrt(mse) if mse is not None else None,
                'group_MAE': gm, 'group_MSE': gs}
        gc = {g: group_loss(rows, lambda r: r['methods'][method]['classification_error'])
              for g, rows in per_group.items()}
        confusion = {key: 0 for key in ('TP', 'TN', 'FP', 'FN')}
        for row in scored:
            truth = row['true_joint_only_obstruction']
            pred = row['methods'][method]['predicted_joint_only_obstruction']
            confusion[('TP' if pred else 'FN') if truth else ('FP' if pred else 'TN')] += 1
        pos = confusion['TP'] + confusion['FN']; neg = confusion['TN'] + confusion['FP']
        summaries[method] = {'continuous': continuous,
            'group_balanced_classification_error': average([v for v in gc.values() if v is not None]),
            'group_classification_error': gc, 'pooled_case_confusion': confusion,
            'pooled_case_balanced_accuracy': (confusion['TP']/pos + confusion['TN']/neg)/2 if pos and neg else None}
    both = len({r['true_joint_only_obstruction'] for r in scored}) == 2
    contributing = sum(bool(rows) for rows in per_group.values())
    reasons = []
    if len(scored) < 10 or contributing < 8:
        reasons.append('prediction_coverage_below_10_squares_8_groups')
    if not both:
        reasons.append('both_true_classes_not_observed')
    comparisons = []
    for baseline in METHODS[1:]:
        for trait in TRAITS:
            a = summaries['additive']['continuous'][trait]['MAE']
            b = summaries[baseline]['continuous'][trait]['MAE']
            passed = b is not None and b > 0 and a <= .9*b
            comparisons.append({'baseline': baseline, 'trait': trait, 'additive_MAE': a,
                'baseline_MAE': b, 'absolute_MAE_difference': a-b if b is not None else None,
                'relative_gain': (b-a)/b if b is not None and b > 0 else None,
                'passes_10_percent': passed})
            if not passed:
                reasons.append('MAE_hurdle_failed:'+baseline+':'+trait)
        a = summaries['additive']['group_balanced_classification_error']
        b = summaries[baseline]['group_balanced_classification_error']
        if b is None or a > b:
            reasons.append('classification_error_worse_or_unavailable:'+baseline)
    return {'weighting': 'group_target_case_loss' if target_balanced else 'group_case_loss',
        'scored_cases': len(scored), 'selected_groups': len(groups), 'contributing_groups': contributing,
        'distinct_targets': len({r['withheld_ID'] for r in scored}), 'both_true_classes': both,
        'groups': {g: {'cases': len(rows), 'distinct_targets': len({r['withheld_ID'] for r in rows})}
                   for g, rows in per_group.items()},
        'summaries': summaries, 'MAE_comparisons': comparisons,
        'discrimination_screen_passed': not reasons, 'failure_reasons': reasons}


def component_detail(field):
    return {'complete': field['complete'], 'point': field['point'], 'bounds': field['bounds']}


def bypass(local, component):
    if local['classification'] != 'joint_only_obstruction':
        return 'not_a_local_joint_only_obstruction'
    if component['bounds'].get('joint_pessimistic') is True:
        return 'certified_observed_component_bypass'
    if component['bounds'].get('joint_optimistic') is False:
        return 'observed_component_optimistic_disconnection'
    return 'observed_component_bypass_unresolved'


def diagonal_diagnostics(descriptions):
    rows = []; flips = {'local_point': 0, 'component_point_or_bound': 0,
                        'local_range': 0, 'component_range': 0,
                        'scenario_local': 0, 'scenario_component_point_or_bound': 0}
    def state(field):
        return (field['point']['classification'] if field['point'] is not None else None,
                field['bounds']['status'])
    for square in descriptions:
        for diagonal in square['diagonals']:
            flips['local_point'] += diagonal['primary']['classification'] != diagonal['local_alternative_tolerance']['classification']
            flips['component_point_or_bound'] += state(diagonal['observed_component']) != state(diagonal['observed_component_alternative_tolerance'])
            flips['local_range'] += diagonal['local_range_certificate']['status'] != diagonal['local_range_alternative_tolerance']['status']
            flips['component_range'] += diagonal['component_range_certificate']['status'] != diagonal['component_range_alternative_tolerance']['status']
            scenarios = []
            pairs = set()
            for scenario in diagonal['preparation_scenarios']:
                pair = (scenario['protein_column'], scenario['growth_column'])
                if pair in pairs:
                    raise ValueError('Duplicate preparation scenario')
                pairs.add(pair)
                flips['scenario_local'] += scenario['local']['classification'] != scenario['local_alternative_tolerance']['classification']
                flips['scenario_component_point_or_bound'] += state(scenario['component']) != state(scenario['component_alternative_tolerance'])
                scenarios.append({**scenario, 'bypass_status_tolerance_zero': bypass(scenario['local'], scenario['component']),
                    'bypass_status_alternative_tolerance': bypass(scenario['local_alternative_tolerance'], scenario['component_alternative_tolerance'])})
            if pairs != {(p, g) for p in PROTEIN for g in GROWTH}:
                raise ValueError('Incomplete12-scenario inventory')
            preparation_consistent = all(s['local']['classification'] == 'joint_only_obstruction' for s in scenarios)
            if type(diagonal['preparation_consistent_joint_obstruction']) is not bool or diagonal['preparation_consistent_joint_obstruction'] != preparation_consistent:
                raise ValueError('Preserved preparation-consistency flag mismatch')
            rows.append({'square_index': square['square_index'], 'seed': square['seed'], 'group': square['group'],
                'endpoint_IDs': diagonal['endpoint_IDs'], 'local_point': diagonal['primary'],
                'local_alternative_tolerance': diagonal['local_alternative_tolerance'],
                'point_tolerance': diagonal['point_tolerance'],
                'observed_component': component_detail(diagonal['observed_component']),
                'observed_component_alternative_tolerance': component_detail(diagonal['observed_component_alternative_tolerance']),
                'bypass_status_tolerance_zero': bypass(diagonal['primary'], diagonal['observed_component']),
                'bypass_status_alternative_tolerance': bypass(diagonal['local_alternative_tolerance'], diagonal['observed_component_alternative_tolerance']),
                'local_range_certificate': diagonal['local_range_certificate'],
                'component_range_certificate': diagonal['component_range_certificate'],
                'interval_tolerance': diagonal['interval_tolerance'],
                'local_range_alternative_tolerance': diagonal['local_range_alternative_tolerance'],
                'component_range_alternative_tolerance': diagonal['component_range_alternative_tolerance'],
                'native_shortfalls': diagonal['native_shortfalls'], 'preparation_scenarios': scenarios})
    return {'eligible_diagonals': len(rows), 'flip_counts': flips, 'records': rows,
        'scope': 'Preserved evaluated observed-component bounds only; no claim about unmeasured synonymous space.'}


def flat_numeric(value, path=()):
    if isinstance(value, dict):
        for key, item in value.items():
            if key not in ('groups', 'group_denominators', 'group_MAE', 'group_MSE', 'group_classification_error'):
                yield from flat_numeric(item, path+(key,))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            yield from flat_numeric(item, path+(str(i),))
    elif value is None or isinstance(value, (int, float, bool)):
        yield '/'.join(path), value


def analyze(evaluation, split, arm):
    units, groups, gf = metadata(split, arm)
    coverage = by_unique_square(evaluation['coverage'], units, gf)
    desc = by_unique_square(evaluation['descriptions'], units, gf)
    records = by_unique_square(evaluation['prediction_evaluation']['records'], units, gf)
    scored = []; case_labels = []; target_values = {}
    for index, row in desc.items():
        unit = units[index]
        if type(row['eligible']) is not bool or row['eligible'] != coverage[index]['eligible']:
            raise ValueError('Coverage/descriptive mismatch')
        expected = {tuple(v) for v in unit['diagonals']}
        actual = [tuple(v['endpoint_IDs']) for v in row['diagonals']]
        if (row['eligible'] and (len(actual) != 2 or set(actual) != expected)) or (not row['eligible'] and actual):
            raise ValueError('Eligible diagonal inventory mismatch')
    for index, record in records.items():
        if record['status'] != 'scored':
            continue
        if not desc[index]['eligible']:
            raise ValueError('Scored ineligible square')
        unit = units[index]; target = unit['withheld_ID']
        position = unit['coordinate_order_00_01_10_11'].index(target)
        endpoints = [unit['coordinate_order_00_01_10_11'][j] for j in (position^1, position^2)]
        truth = record['true_joint_only_obstruction']
        if type(truth) is not bool or len(record['observed_target']) != 2 or set(record['methods']) != set(METHODS):
            raise ValueError('Invalid scored schema')
        actual_diagonal = next(d for d in desc[index]['diagonals'] if set(d['endpoint_IDs']) == set(endpoints))
        if truth != (actual_diagonal['primary']['classification'] == 'joint_only_obstruction'):
            raise ValueError('Scored endpoint truth differs from preserved description')
        if target in target_values and target_values[target] != record['observed_target']:
            raise ValueError('Shared target observation changed between cases')
        target_values[target] = record['observed_target']
        for method in METHODS:
            values = record['methods'][method]
            if type(values['predicted_joint_only_obstruction']) is not bool:
                raise ValueError('Nonboolean classification')
            if values['classification_error'] != int(values['predicted_joint_only_obstruction'] != truth):
                raise ValueError('Classification loss mismatch')
            if type(values['classification_error']) is not int:
                raise ValueError('Classification loss must be an integer indicator')
            for key in ('predicted_target', 'signed_error', 'absolute_error', 'squared_error'):
                if len(values[key]) != 2:
                    raise ValueError('Wrong trait dimension')
            for k in range(2):
                error = finite(values['predicted_target'][k])-finite(record['observed_target'][k])
                for key, expected in [('signed_error', error), ('absolute_error', abs(error)), ('squared_error', error*error)]:
                    if not close(finite(values[key][k]), expected):
                        raise ValueError('Preserved loss arithmetic mismatch')
        row = {**record, 'withheld_ID': target, 'classification_endpoint_IDs': endpoints}
        scored.append(row)
        case_labels.append({'square_index': index, 'group': row['group'], 'withheld_ID': target,
            'endpoint_IDs': endpoints, 'true_label': truth,
            'method_labels': {m: row['methods'][m]['predicted_joint_only_obstruction'] for m in METHODS}})
    full = {'descriptive': local_summary(list(desc.values()), groups),
            'case_balanced_prediction': score_summary(scored, groups),
            'target_balanced_prediction': score_summary(scored, groups, True)}
    leave_out = []
    for group in groups:
        kept = [g for g in groups if g != group]
        sr = [r for r in scored if r['group'] != group]
        dr = [r for r in desc.values() if r['group'] != group]
        leave_out.append({'omitted_group': group, 'descriptive': local_summary(dr, kept),
            'case_balanced_prediction': score_summary(sr, kept),
            'target_balanced_prediction': score_summary(sr, kept, True)})
    ranges = {}
    for row in leave_out:
        for key, value in flat_numeric({k: v for k, v in row.items() if k != 'omitted_group'}):
            ranges.setdefault(key, []).append(value)
    ranges = {key: {'minimum': min(known) if known else None, 'maximum': max(known) if known else None,
                    'defined_rows': len(known), 'undefined_rows': len(values)-len(known), 'all_rows': len(leave_out)}
              for key, values in ranges.items() for known in [[v for v in values if v is not None]]}
    return {'status': 'supplemental_diagnostics_complete', 'interpretation_permitted': True, 'arm': arm,
        'full_arm': full, 'prediction_case_labels_preserved': case_labels,
        'all_prediction_record_statuses': [{'square_index': i, 'status': r['status'], 'reason': r.get('reason', r.get('missing_columns'))}
                                           for i, r in sorted(records.items())],
        'leave_one_group_out': leave_out, 'all_leave_one_out_minmax_ranges': ranges,
        'observed_component_diagnostics': diagonal_diagnostics(list(desc.values())),
        'scope': 'Supplement only. Preserves primary evaluation; no new model, ensemble, favorable selection or confidence interval.'}


def run(evaluation_path, evaluation_sha256, split_path, arm, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError('Output must be a new path')
    eb = Path(evaluation_path).read_bytes()
    if digest(eb) != evaluation_sha256:
        raise ValueError('Evaluation SHA256 mismatch')
    evaluation = load(eb)
    if evaluation.get('arm') != arm:
        raise ValueError('Arm access refused')
    sb = Path(split_path).read_bytes(); split = load(sb)
    if evaluation.get('source_hashes', {}).get('split') != digest(sb):
        raise ValueError('Split SHA256 differs from evaluation provenance')
    if evaluation.get('source_hashes', {}).get('runner') != RUNNER_SHA256:
        raise ValueError('Unreviewed evaluation schema/runner hash')
    metadata(split, arm)
    if evaluation.get('interpretation_permitted') is not True:
        result = {'status': 'provenance_hold_no_interpretation', 'interpretation_permitted': False,
                  'arm': arm, 'source_status': evaluation.get('status'), 'holds': evaluation.get('holds', [])}
    elif evaluation.get('status') != 'arm_evaluation_complete':
        raise ValueError('Unexpected evaluation status')
    else:
        result = analyze(evaluation, split, arm)
    result['input_hashes'] = {'evaluation': digest(eb), 'split': digest(sb),
                            'diagnostic_script': digest(Path(__file__).read_bytes())}
    result['numeric_identity_check_tolerance'] = {'absolute': 1e-12, 'relative': 1e-12,
        'purpose': 'Verify stored loss arithmetic only; never a biological threshold.'}
    output.mkdir(parents=True, exist_ok=False)
    with (output/'diagnostics.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False, sort_keys=True)
        stream.write('\n')
    return result['status']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation', required=True)
    parser.add_argument('--evaluation-sha256', required=True)
    parser.add_argument('--split', required=True)
    parser.add_argument('--allow-arm', choices=('development', 'reserved_validation'), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(run(args.evaluation, args.evaluation_sha256, args.split, args.allow_arm, args.output))
