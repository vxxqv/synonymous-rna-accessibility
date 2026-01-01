"""Read only the authorized development summaries; aggregate existing labels, not models."""
import collections
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(test, message):
    if not test:
        raise RuntimeError(message)


def main():
    root = Path.cwd()
    owned = root / 'work/joint_function_cycle5/development_interpretation_review'
    sources = {
        'evaluation': ('outputs/synonymous-rna-accessibility/results/cycle5/development/evaluation/evaluation.json',
                       '756ff513c7c4f628110b07de19c0608ef5b1854da548b0227bca44db75e8f8fa'),
        'diagnostics': ('outputs/synonymous-rna-accessibility/results/cycle5/development/diagnostics/diagnostics.json',
                        'c5a311f6470481ca6c0de474d7e80fdd348e100334822def01e2a342655f66db'),
        'protocol': ('outputs/synonymous-rna-accessibility/analysis/cycle5/cambray_pilot/protocol_v1.md',
                     '0f51ccc808559ea69df3bdf294520d81babb67bf14ca92ce2abb6fcc665760ac'),
        'measurement_limits': ('work/joint_function_cycle5/measurement_limits/report.md',
                               '1171ef7242e299214c8449ee0222b5f50f1a982bdb73b1b36169aab77d7d3740'),
    }
    inventory = {}
    for name, (rel, expected) in sources.items():
        p = root / rel
        actual = sha(p)
        require(actual == expected, 'Source hash mismatch: ' + name)
        inventory[name] = {'path': rel, 'sha256': actual, 'bytes': p.stat().st_size}
    e = json.loads((root / sources['evaluation'][0]).read_text())
    d = json.loads((root / sources['diagnostics'][0]).read_text())
    require(e['arm'] == d['arm'] == 'development', 'Arm mismatch')
    require(e['interpretation_permitted'] and d['interpretation_permitted'], 'Interpretation hold')
    require(d['input_hashes']['evaluation'] == inventory['evaluation']['sha256'], 'Diagnostic source mismatch')
    require(e['source_hashes']['protocol'] == inventory['protocol']['sha256'], 'Protocol provenance mismatch')

    diagonals = []
    for square in e['descriptions']:
        for row in square['diagonals']:
            diagonals.append({'square_index': square['square_index'], 'group': square['group'],
                              'seed': square['seed'], 'isolated_observed_square': square['isolated_observed_square'],
                              **row})
    labels = collections.Counter(r['primary']['classification'] for r in diagonals)
    range_counts = collections.Counter(r['local_range_certificate']['status'] for r in diagonals)
    point_bounds = collections.Counter(r['observed_component']['bounds']['status'] for r in diagonals)
    groups = sorted({r['group'] for r in diagonals}, key=int)
    consistent = sum(r['preparation_consistent_joint_obstruction'] for r in diagonals)
    require(len(diagonals) == 24 and len(groups) == 10, 'Coverage mismatch')
    require(dict(labels) == {'jointly_connected': 20, 'single_trait_disconnection': 2,
                              'joint_only_obstruction': 2}, 'Point count mismatch')
    require(consistent == 1, 'Series-consistent count mismatch')
    require(range_counts['certified_joint_only_obstruction'] == 0, 'Range count mismatch')
    require(e['primary_summary']['eligible_diagonals'] == len(diagonals), 'Primary denominator mismatch')
    require(e['primary_summary']['preparation_consistent_diagonals'] == consistent, 'Primary consistency mismatch')
    scenario_counts = {}
    for row in diagonals:
        require(len(row['preparation_scenarios']) == 12, 'Incomplete scenario set')
        for scenario in row['preparation_scenarios']:
            key = scenario['protein_column'] + '|' + scenario['growth_column']
            scenario_counts.setdefault(key, collections.Counter())[scenario['local']['classification']] += 1
    require(len(scenario_counts) == 12, 'Scenario cardinality mismatch')
    tables = d['observed_component_diagnostics']['count_tables']
    for key, counts in scenario_counts.items():
        actual = tables['named_preparation_scenarios'][key]['zero_tolerance']['local_classifications']['counts']
        require(all(actual[k] == v for k, v in counts.items()), 'Scenario counts mismatch: ' + key)
        require(sum(counts.values()) == 24, 'Scenario denominator mismatch')
    prediction_statuses = collections.Counter(r['status'] for r in e['prediction_evaluation']['records'])
    require(dict(prediction_statuses) == {'scored': 11, 'missing_known_replicates': 11,
                                        'excluded_shared_corner_mask': 1}, 'Prediction status mismatch')
    require(not e['prediction_evaluation']['discrimination_screen_passed'], 'Unexpected prediction pass')
    require(len(d['leave_one_group_out']) == 16, 'Incomplete group-deletion set')
    require(all(not r['case_balanced_prediction']['discrimination_screen_passed']
                for r in d['leave_one_group_out']), 'Unexpected deletion pass')
    require(all(v == 0 for v in d['observed_component_diagnostics']['flip_counts'].values()),
            'Unexpected tolerance flip')
    require(d['full_arm']['case_balanced_prediction']['summaries'] ==
            d['full_arm']['target_balanced_prediction']['summaries'], 'Weighting identity mismatch')

    loo = []
    for r in d['leave_one_group_out']:
        v = r['case_balanced_prediction']
        loo.append({'omitted_group': r['omitted_group'], 'descriptive': r['descriptive'],
                    'prediction': {k: v[k] for k in ['scored_cases', 'distinct_targets', 'contributing_groups',
                                                       'both_true_classes', 'discrimination_screen_passed',
                                                       'failure_reasons', 'MAE_comparisons']}})
    ledger = {
        'scope': 'Complete authorized development inventory; no reserved or raw-source access.',
        'coverage_all_23': e['coverage'], 'descriptive_diagonals_all_24': diagonals,
        'groups_all_16': d['full_arm']['descriptive']['group_denominators'],
        'prediction_statuses_all_23': d['all_prediction_record_statuses'],
        'full_arm_metrics': d['full_arm'], 'leave_one_group_out_all_16': loo,
        'component_scenario_range_tables': tables,
        'all_leave_one_out_minmax_ranges': d['all_leave_one_out_minmax_ranges'],
    }
    checks = {
        'status': 'development_only_interpretation_review_complete',
        'verdict': 'Unchanged reserved descriptive question remains useful subject to root decision and report freeze.',
        'source_inventory': inventory,
        'independent_preserved_label_aggregations': {
            'eligible_diagonals': len(diagonals), 'eligible_groups': groups,
            'primary_labels': dict(labels), 'range_statuses': dict(range_counts),
            'component_point_bounds': dict(point_bounds), 'series_consistent_diagonals': consistent,
            'named_scenario_local_counts': {k: dict(v) for k, v in scenario_counts.items()},
            'prediction_statuses': dict(prediction_statuses),
            'mean_check_statuses': dict(collections.Counter(r['status'] for r in e['mean_diagnostics'])),
            'all_16_group_deletion_prediction_screens_failed': True,
            'all_recorded_tolerance_flips_zero': True,
            'target_and_case_balanced_summary_identity': True,
        },
        'exposure_ledger': {
            'newly_authorized': 'All preserved DEVELOPMENT outcomes, exclusions, margins and prediction errors in the two summaries.',
            'prior_exposed': 'Source Methods and published aggregates, source-limits reports, frozen protocol; historical seed-1 exclusion unchanged.',
            'not_opened': ['raw phenotype ZIP', 'reserved missingness', 'reserved outcomes', 'reserved class balance'],
            'no_actions': ['fitting', 'threshold/model change', 'new literature', 'manuscript edit', 'commit', 'publication', 'reserved execution'],
        },
        'review_limits': ['Separately authored with inherited context, not external/blinded review.',
                          'Aggregates preserved records; does not reconstruct assay source or calibrated uncertainty.',
                          'The reviewer previously authored the independently reviewed supplemental helper.'],
    }
    for name, obj in [('development_review_ledger.json', ledger), ('checks.json', checks)]:
        target = owned / name
        require(not target.exists(), 'Refusing output overwrite: ' + name)
        target.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    for name, (rel, expected) in sources.items():
        require(sha(root / rel) == expected, 'Input changed during review: ' + name)
    print(json.dumps({'status': 'passed', 'sources_unchanged': True, 'diagonals': len(diagonals),
                      'outputs': {n: sha(owned / n) for n in ['review.md', 'checks.json',
                                                             'development_review_ledger.json']}}))


if __name__ == '__main__':
    main()
