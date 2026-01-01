"""Frozen exact arithmetic on previously reconstructed COMT assay summaries."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import time

MAPS = ('code', 'geo')
ASSAYS = ('R', 'D', 'F')
REVERSED = ('reversal_R_up_T_down', 'reversal_R_down_T_up')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def number(value):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError('Finite numeric stored input required')
    return Q(value)

def classify(rlo, rhi, tlo, thi, ties=False):
    if rlo > 0 and thi < 0: return REVERSED[0]
    if rhi < 0 and tlo > 0: return REVERSED[1]
    if rlo > 0 and tlo > 0: return 'same_direction_a_higher'
    if rhi < 0 and thi < 0: return 'same_direction_b_higher'
    return 'unresolved_with_ties' if ties else 'unresolved_no_ties'

def inspect_group(group):
    if group is None: return {'reason': 'missing_group'}
    samples, values = group['samples'], group['values']
    if not isinstance(samples, list) or not all(type(s) is str for s in samples):
        raise ValueError('Invalid sample identities')
    if len(samples) != len(set(samples)) or len(samples) != len(values):
        raise ValueError('Invalid sample/value alignment')
    mask = [x is not None for x in values]
    for v in values:
        if v is not None: number(v)
    if type(group['technical_finite']) is not int or group['technical_finite'] != sum(mask):
        raise ValueError('Finite-count mismatch')
    value = group['median']
    if sum(mask) == 0:
        if value is not None: raise ValueError('Median without finite contributors')
        return {'reason': 'no_finite_group', 'samples': samples, 'mask': mask}
    if value is None: raise ValueError('Missing stored median with finite contributors')
    return {'median': number(value), 'samples': samples, 'mask': mask}

def group_map(row, assay):
    groups = row['assays'][assay]['groups']
    if not isinstance(groups, list): raise ValueError('Group list required')
    result = {}
    for group in groups:
        g = group['bio']
        if type(g) is not int or g not in range(1, 5) or g in result:
            raise ValueError('Duplicate or invalid group')
        result[g] = group
    return result

def pair_map(a, b):
    z, provenance, issues = {}, {}, []
    for assay in ASSAYS:
        aa, bb = group_map(a, assay), group_map(b, assay)
        z[assay], provenance[assay] = [], []
        for g in range(1, 5):
            av, bv = inspect_group(aa.get(g)), inspect_group(bb.get(g))
            reason = av.get('reason') or bv.get('reason')
            if reason is None and av['samples'] != bv['samples']: reason = 'technical_sample_list_mismatch'
            if reason is None and av['mask'] != bv['mask']: reason = 'technical_contribution_mask_mismatch'
            z[assay].append(None if reason else av['median'] - bv['median'])
            provenance[assay].append({'group': g, 'samples_a': av.get('samples'), 'samples_b': bv.get('samples'),
                                      'finite_mask_a': av.get('mask'), 'finite_mask_b': bv.get('mask'), 'issue': reason})
            if reason: issues.append({'assay': assay, 'group': g, 'reason': reason})
    cases = []
    for i, j, k in itertools.product(range(4), repeat=3):
        r = None if z['R'][i] is None or z['D'][j] is None else z['R'][i] - z['D'][j]
        t = None if z['D'][j] is None or z['F'][k] is None else z['D'][j] - z['F'][k]
        cases.append({'R_group': i+1, 'D_group': j+1, 'F_group': k+1,
                      'DeltaR': None if r is None else str(r), 'DeltaT': None if t is None else str(t)})
    complete = not issues
    if complete:
        rr, tt = [Q(c['DeltaR']) for c in cases], [Q(c['DeltaT']) for c in cases]
        vals = (min(rr), max(rr), min(tt), max(tt))
        ties = any(v == 0 for v in rr + tt)
        bounds = dict(zip(('R_lower', 'R_upper', 'T_lower', 'T_upper'), map(str, vals)))
        category = classify(*vals, ties)
        limiting = {name: [idx for idx, value in enumerate(rr if name.startswith('R') else tt)
                          if value == Q(bound)] for name, bound in bounds.items()}
    else:
        bounds = None; category = 'unresolved_missing_or_misaligned'; ties = None; limiting = None
    return {'complete': complete, 'category': category, 'issues': issues,
            'same_assay_differences': {a: [None if x is None else str(x) for x in zz] for a, zz in z.items()},
            'source_groups': provenance, 'joint_cases': cases, 'bounds': bounds,
            'any_exact_zero': ties, 'extremum_case_indices_zero_based': limiting}

def pair_result(meta, inputs):
    permap = {m: pair_map(inputs[m][meta['a']], inputs[m][meta['b']]) for m in MAPS}
    complete = all(v['complete'] for v in permap.values())
    if complete:
        bounds = {n: (min if n.endswith('lower') else max)(Q(v['bounds'][n]) for v in permap.values())
                  for n in ('R_lower', 'R_upper', 'T_lower', 'T_upper')}
        category = classify(*[bounds[n] for n in ('R_lower', 'R_upper', 'T_lower', 'T_upper')],
                            any(v['any_exact_zero'] for v in permap.values()))
        margin = (min(bounds['R_lower'], -bounds['T_upper']) if category == REVERSED[0]
                  else min(-bounds['R_upper'], bounds['T_lower']) if category == REVERSED[1] else Q(0))
        radius = str(margin / 4)
        bounds = {k: str(v) for k, v in bounds.items()}
    else:
        bounds = None; category = 'unresolved_missing_or_misaligned'; radius = None
    return {**meta, 'maps': permap, 'complete': complete, 'category': category, 'combined_bounds': bounds,
            'unstructured_log_radius_exclusive': radius,
            'map_categories_differ': len({v['category'] for v in permap.values()}) != 1}

def validate(data, contract):
    if set(data) != set(MAPS): raise ValueError('Exactly both source maps required')
    targets = {r['hgvs_nt']: r for r in contract['targets']}
    members = {r['sample']: r for r in contract['members'].values()}
    indexed = {}
    for mp in MAPS:
        rows = data[mp]
        indexed[mp] = {r['hgvs_nt']: r for r in rows}
        if len(rows) != len(indexed[mp]) or set(indexed[mp]) != set(targets):
            raise ValueError('Target coverage or uniqueness mismatch')
        for ident, row in indexed[mp].items():
            expected = targets[ident]
            for key in ('site', 'ROI', 'native_codon_site', 'is_SNV', 'edit_key'):
                if row[key] != expected[key]: raise ValueError('Target metadata mismatch: ' + key)
            for assay in ASSAYS:
                for g, obj in group_map(row, assay).items():
                    inspect_group(obj)
                    for sample in obj['samples']:
                        s = members.get(sample)
                        if s is None or (s['trait'], s['ROI'], s[mp + '_bio']) != (assay, row['ROI'], g):
                            raise ValueError('Unrecognized or misassigned sample identity')
    return indexed

def summaries(records, targets):
    output = {}
    for subset in ('all', 'SNV'):
        rr = records if subset == 'all' else [r for r in records if r['both_SNV']]
        vv = targets if subset == 'all' else [r for r in targets if r['is_SNV']]
        site_rows = []
        for site in sorted({int(v['native_codon_site']) for v in vv}):
            sr = [r for r in rr if r['site'] == site]
            nr = sum(r['category'] in REVERSED for r in sr)
            site_rows.append({'site': site, 'variant_count': sum(int(v['native_codon_site']) == site for v in vv),
                              'pair_count': len(sr), 'complete_pairs': sum(r['complete'] for r in sr),
                              'uniform_reversals': nr, 'category_counts': dict(sorted(Counter(r['category'] for r in sr).items())),
                              'reversal_fraction': str(Q(nr, len(sr))) if sr else None})
        informative = [s for s in site_rows if s['pair_count']]
        output[subset] = {'variants': len(vv), 'sites': len(site_rows), 'pairs': len(rr),
                          'sites_with_pairs': len(informative), 'singleton_sites': len(site_rows)-len(informative),
                          'complete_pairs': sum(r['complete'] for r in rr),
                          'joint_case_records': 128 * len(rr),
                          'category_counts': dict(sorted(Counter(r['category'] for r in rr).items())),
                          'uniform_reversal_pairs': sum(r['category'] in REVERSED for r in rr),
                          'sites_with_uniform_reversal': sum(s['uniform_reversals'] > 0 for s in site_rows),
                          'equal_informative_site_reversal_fraction': str(sum((Q(s['reversal_fraction']) for s in informative), Q(0)) / len(informative)) if informative else None,
                          'map_category_disagreements': sum(r['map_categories_differ'] for r in rr),
                          'pair_counts_by_ROI': dict(sorted(Counter(r['ROI'] for r in rr).items())),
                          'sites_detail': site_rows}
    return output

def execute(freeze_path, input_path, inventory_path, contract_path, out):
    if not __debug__: raise RuntimeError('Optimized Python unsupported')
    out = Path(out)
    if out.exists(): raise FileExistsError('Fresh output directory required before input decoding')
    freeze_path = Path(freeze_path)
    freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
    if freeze['status'] != 'FROZEN_SYNTHETIC_PASS_ROOT_AUTHORIZED_EXPLORATORY_RUN':
        raise ValueError('Frozen authorization required')
    required = {'analyze.py', 'verify.py', 'test_synthetic.py', 'proposed_protocol.md', 'synthetic_receipt.json'}
    if not required <= set(freeze['owned_files']): raise ValueError('Incomplete implementation freeze')
    for name, digest in freeze['owned_files'].items():
        if Path(name).name != name: raise ValueError('Unsafe freeze member')
        if sha(freeze_path.parent/name) != digest: raise ValueError('Frozen source changed: ' + name)
    if sha(__file__) != freeze['owned_files']['analyze.py']: raise ValueError('Executing source differs from freeze')
    paths = {'input': Path(input_path), 'inventory': Path(inventory_path), 'contract': Path(contract_path)}
    before = {k: sha(p) for k, p in paths.items()}
    if before != freeze['input_hashes']: raise ValueError('Input hash mismatch before outcome decoding')
    out.mkdir(parents=True)
    start = time.perf_counter()
    dump(out/'run_started.json', {'freeze_sha256': sha(freeze_path), 'input_hashes': before, 'status': 'STARTED'})
    try:
        data = json.loads(paths['input'].read_text(encoding='utf-8'))
        inventory = json.loads(paths['inventory'].read_text(encoding='utf-8'))
        contract = json.loads(paths['contract'].read_text(encoding='utf-8'))
        inputs = validate(data, contract)
        if len(inventory['pairs']) != 134: raise ValueError('Expected complete 134-pair inventory')
        expected = {(a['hgvs_nt'], b['hgvs_nt']) for a, b in itertools.combinations(sorted(contract['targets'], key=lambda r: r['hgvs_nt'].encode('utf-8')), 2)
                    if a['native_codon_site'] == b['native_codon_site']}
        if {(r['a'], r['b']) for r in inventory['pairs']} != expected: raise ValueError('Incomplete or extra pair inventory')
        records = [pair_result(p, inputs) for p in inventory['pairs']]
        result = {'schema': 'COMT_local_pair_order_v1', 'status': 'COMPLETED_EXPLORATORY_FINITE_COMPARISON',
                  'assertion_tested': 'Fixed local allele pair and fixed opposite orientation for all128 declared map/preparation combinations.',
                  'units': 'natural_log_source_normalized_proxy_difference',
                  'exactness': 'Exact rational arithmetic of stored finite binary floats; not exact biological quantities.',
                  'summaries': summaries(records, contract['targets']), 'pairs': records,
                  'nonclaims': ['independent biological replication', 'confidence intervals', 'physical protein amount', 'causal compensation', 'fitness trade-off', 'accessibility barrier'],
                  'frozen_prediction_result_unchanged': True}
        dump(out/'result.json', result)
        after = {k: sha(p) for k, p in paths.items()}
        if before != after: raise RuntimeError('Inputs changed during run')
        receipt = {'status': 'COMPLETED', 'freeze_sha256': sha(freeze_path), 'input_hashes_before': before,
                   'input_hashes_after': after, 'result_sha256': sha(out/'result.json'),
                   'pairs': len(records), 'joint_cases': sum(len(m['joint_cases']) for r in records for m in r['maps'].values()),
                   'python': platform.python_version(), 'elapsed_seconds': time.perf_counter()-start,
                   'argv': {'freeze': str(freeze_path), 'input': str(input_path), 'inventory': str(inventory_path), 'contract': str(contract_path), 'out': str(out)}}
        dump(out/'receipt.json', receipt)
        return result['summaries']
    except Exception as e:
        dump(out/'failure.json', {'status': 'FAILED_STOP', 'exception': type(e).__name__, 'message': str(e)})
        raise

if __name__ == '__main__':
    if not __debug__: raise SystemExit('Optimized Python unsupported')
    p = argparse.ArgumentParser()
    for name in ('freeze', 'input', 'inventory', 'contract', 'out'): p.add_argument('--'+name, required=True)
    a = p.parse_args()
    summary = execute(a.freeze, a.input, a.inventory, a.contract, a.out)
    print(json.dumps({s: {k:v for k,v in r.items() if k != 'sites_detail'} for s,r in summary.items()}, indent=2))
