"""Separate COMT downstream reference. No author scientific function is imported.

RF row masks and selective extraction are explicitly conditioned inputs.
Only explicit, hash-authorized files are read by the command line.
"""
if not __debug__:
    raise RuntimeError('Optimized Python is prohibited before input access')

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

MODELS = ('pooled_median', 'ROI_median', 'RNA_PAVA', 'RNA_linear', 'RNA_load_linear')
BASELINES = MODELS[:-1]
NUMERIC_RELATIVE_TOLERANCE = 2e-12
PROTOCOL_SHA256 = '97aea7e56ddf9c74a2fcbf4a17ce2bbd9280c8417d2534997835878b2a368d61'
INTEGRATION_SHA256 = 'c480ca92ce36b6e95849fceaea3e6bf49d14bcd7842b1e5fe1efd543ff33f4d2'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rational(x):
    if isinstance(x, bool) or x is None:
        raise ValueError('Boolean/None is not a numerical input')
    if isinstance(x, float) and not math.isfinite(x):
        raise ValueError('Nonfinite input')
    return x if isinstance(x, F) else F(str(x))


def median(values):
    s = sorted(values)
    if not s:
        raise ValueError('Empty median')
    q, r = divmod(len(s), 2)
    return s[q] if r else (s[q - 1] + s[q]) / 2


def weighted_median(values, weights):
    d = defaultdict(F)
    for y, w in zip(values, weights, strict=True):
        if w <= 0:
            raise ValueError('Nonpositive weight')
        d[y] += w
    levels = sorted(d)
    if not levels:
        raise ValueError('Empty weighted median')
    half, mass = sum(d.values(), F()) / 2, F()
    for i, y in enumerate(levels):
        mass += d[y]
        if mass > half:
            return y
        if mass == half:
            return (y + levels[i + 1]) / 2
    raise AssertionError('Unreachable median')


def weights_for_sites(rows):
    counts = Counter(r['site'] for r in rows)
    return [F(1, counts[r['site']]) for r in rows]


def isotonic(x, y, w):
    """Min-max characterization, not the production block-merging routine.

    z_i = max_{a<=i} min_{b>=i} weighted_mean(y_a,...,y_b), after ties.
    """
    table = defaultdict(lambda: [F(), F()])
    for xi, yi, wi in zip(x, y, w, strict=True):
        table[xi][0] += wi
        table[xi][1] += wi * yi
    levels = sorted(table)
    if not levels:
        raise ValueError('Empty isotonic fit')
    cw, cy = [F()], [F()]
    for a in levels:
        cw.append(cw[-1] + table[a][0])
        cy.append(cy[-1] + table[a][1])
    fitted = []
    for i in range(len(levels)):
        candidates = []
        for a in range(i + 1):
            candidates.append(min((cy[b + 1] - cy[a]) / (cw[b + 1] - cw[a])
                                  for b in range(i, len(levels))))
        fitted.append(max(candidates))
    return levels, fitted


def interpolate(fit, x):
    a, y = fit
    i = bisect_left(a, x)
    if i == 0:
        return y[0]
    if i == len(a):
        return y[-1]
    return y[i - 1] + (x - a[i - 1]) * (y[i] - y[i - 1]) / (a[i] - a[i - 1])


def constrained_linear(x, z, y, w):
    """Exact least squares by residualizing RNA against load, then KKT clamp.

    Returns (intercept, RNA coefficient, load coefficient, centered rank).
    The rank-deficient extrapolation convention is fixed by the protocol.
    """
    total = sum(w, F())
    if total <= 0:
        raise ValueError('Empty fit')
    mx, mz, my = [sum(wi * vi for wi, vi in zip(w, v)) / total for v in (x, z, y)]
    xc, zc, yc = [[v - m for v in vs] for vs, m in ((x, mx), (z, mz), (y, my))]
    dot = lambda a, b: sum(wi * ai * bi for wi, ai, bi in zip(w, a, b))
    zz, xx = dot(zc, zc), dot(xc, xc)
    if zz:
        load_x, load_y = dot(zc, xc) / zz, dot(zc, yc) / zz
        u = [a - load_x * b for a, b in zip(xc, zc)]
        uu = dot(u, u)
        a = max(F(), dot(u, yc) / uu) if uu else F()
        b = load_y - a * load_x
        rank = 2 if uu else 1
    elif xx:
        a, b, rank = max(F(), dot(xc, yc) / xx), F(), 1
    else:
        a, b, rank = F(), F(), 0
    return my - a * mx - b * mz, a, b, rank


def fit_and_predict(training, assessment):
    # This API deliberately refuses assessment targets.
    expected = {'hgvs_nt', 'site', 'ROI', 'R', 'L'}
    if any(set(r) != expected for r in assessment):
        raise ValueError('Assessment schema must not contain a response')
    if not training:
        raise ValueError('No training rows')
    w = weights_for_sites(training)
    x, z, y = [[rational(r[name]) for r in training] for name in ('R', 'L', 'T')]
    pooled = weighted_median(y, w)
    roi = {}
    for k in {r['ROI'] for r in training}:
        idx = [i for i, r in enumerate(training) if r['ROI'] == k]
        roi[k] = weighted_median([y[i] for i in idx], [w[i] for i in idx])
    if any(r['ROI'] not in roi for r in assessment):
        raise ValueError('Missing assessment ROI in training')
    iso = isotonic(x, y, w)
    lin = constrained_linear(x, [F()] * len(x), y, w)
    joint = constrained_linear(x, z, y, w)
    result = {}
    for row in assessment:
        xr, zl = rational(row['R']), rational(row['L'])
        result[row['hgvs_nt']] = {
            'pooled_median': pooled, 'ROI_median': roi[row['ROI']],
            'RNA_PAVA': interpolate(iso, xr),
            'RNA_linear': lin[0] + lin[1] * xr,
            'RNA_load_linear': joint[0] + joint[1] * xr + joint[2] * zl,
        }
    return result, {'weights': w, 'RNA_linear': lin, 'RNA_load_linear': joint,
                    'training_ids': [r['hgvs_nt'] for r in training],
                    'training_sites': sorted({r['site'] for r in training}),
                    'R_range': (min(x), max(x)), 'L_range': (min(z), max(z))}


def complete(row, load):
    return all(row['contrasts'][k]['complete'] for k in ('R', 'T', load))


def cross_validate(rows, load):
    byid = {r['hgvs_nt']: r for r in rows}
    if len(byid) != len(rows):
        raise ValueError('Duplicate variant')
    sites = defaultdict(set)
    for r in rows:
        if type(r['assessment_fold']) is not int or r['assessment_fold'] not in range(5):
            raise ValueError('Invalid frozen assessment fold')
        sites[r['site']].add(r['assessment_fold'])
    if any(len(v) != 1 for v in sites.values()):
        raise ValueError('A codon site crosses folds')
    available = [r for r in rows if complete(r, load)]
    output, fits = {}, {}
    unavailable = {r['hgvs_nt'] for r in rows if not complete(r, load)}
    for fold in range(5):
        tr = [r for r in available if r['assessment_fold'] != fold]
        te = [r for r in available if r['assessment_fold'] == fold]
        if not te:
            fits[fold] = {'status': 'no_complete_assessment'}
            continue
        train = [dict(hgvs_nt=r['hgvs_nt'], site=r['site'], ROI=r['ROI'],
                      R=r['contrasts']['R']['center'], L=r['contrasts'][load]['center'],
                      T=r['contrasts']['T']['center']) for r in tr]
        test = [{k: v for k, v in r.items() if k != 'T'} for r in
                [dict(hgvs_nt=r['hgvs_nt'], site=r['site'], ROI=r['ROI'],
                      R=r['contrasts']['R']['center'], L=r['contrasts'][load]['center']) for r in te]]
        try:
            pred, fit = fit_and_predict(train, test)
        except ValueError as e:
            # Record the entire frozen failed fold. Never drop only high-error alleles.
            unavailable.update(r['hgvs_nt'] for r in te)
            fits[fold] = {'status': 'fit_failure', 'reason': str(e)}
            continue
        output.update(pred)
        fits[fold] = {'status': 'predicted', **fit}
    return output, fits, unavailable


def score(pred, rows):
    byid = {r['hgvs_nt']: r for r in rows}
    errors, sites = {}, defaultdict(list)
    for vid, pp in pred.items():
        r = byid[vid]
        errors[vid] = {m: abs(pp[m] - rational(r['contrasts']['T']['center'])) for m in MODELS}
        sites[r['site']].append(vid)
    per_site = {s: {m: sum(errors[v][m] for v in vv) / len(vv) for m in MODELS} for s, vv in sites.items()}
    mae = {m: sum(p[m] for p in per_site.values()) / len(per_site) if per_site else None for m in MODELS}
    return {'variants': len(errors), 'sites': len(sites), 'site_MAE': mae,
            'per_site_errors': per_site, 'variant_errors': errors}


def criterion(scores):
    coverage = scores['variants'] >= 100 and scores['sites'] >= 40
    v = scores['site_MAE']
    defined = all(v[k] is not None and v[k] > 0 for k in BASELINES)
    flags = [v['RNA_load_linear'] is not None and v[k] is not None and v[k] > 0
             and 10 * v['RNA_load_linear'] <= 9 * v[k] for k in BASELINES]
    return {'coverage_pass': coverage, 'relative_comparison_defined': defined,
            'gain_pass_each': flags, 'pass': coverage and defined and all(flags)}


def sign_table(rows):
    values = [-1, 0, 1, 'unknown']
    out = {str(a) + '/' + str(b): 0 for a in values for b in values}
    for r in rows:
        labels = []
        for trait in ('R', 'T'):
            a = r['contrasts'][trait]
            if not a['complete']:
                labels.append('unknown')
            else:
                labels.append(1 if a['lower'] > 0 else -1 if a['upper'] < 0 else 0)
        out[str(labels[0]) + '/' + str(labels[1])] += 1
    return out


def downstream(contrasts):
    out = {'cases': {}}
    for mapping in ('code', 'geo'):
        rows = contrasts[mapping]
        out['cases'][mapping] = {}
        for mode, load in (('heavy', 'L'), ('mid', 'L_mid')):
            pred, fits, missing = cross_validate(rows, load)
            scores = score(pred, rows)
            sensitivities = {}
            for name, condition in [('SNV', lambda r: r['is_SNV']),
                                    ('ROI_1', lambda r: int(r['ROI']) == 1),
                                    ('ROI_2', lambda r: int(r['ROI']) == 2)]:
                selected = {r['hgvs_nt'] for r in rows if condition(r)}
                sensitivities[name] = score({k: v for k, v in pred.items() if k in selected}, rows)
            out['cases'][mapping][mode] = {'predictions': pred, 'fits': fits,
                'input_centers': {r['hgvs_nt']: {'R': r['contrasts']['R']['center'],
                    'L': r['contrasts'][load]['center'], 'fold': r['assessment_fold']} for r in rows},
                'scores': scores, 'sensitivities': sensitivities, 'criterion': criterion(scores),
                'unavailable': sorted(missing),
                'unavailable_details': {r['hgvs_nt']: {
                    'site': r['site'], 'ROI': r['ROI'],
                    'missing_contrasts': [k for k in ('R', 'T', load) if not r['contrasts'][k]['complete']],
                    'failed_design': complete(r, load)} for r in rows if r['hgvs_nt'] in missing},
                'fixed_denominators': {'all_variants': len(rows), 'all_sites': len({r['site'] for r in rows}),
                    'SNVs': sum(bool(r['is_SNV']) for r in rows), 'SNV_sites': len({r['site'] for r in rows if r['is_SNV']})},
                'RNA_T_sign_confusion': sign_table(rows)}
    out['primary_operational_expansion_pass'] = all(out['cases'][m]['heavy']['criterion']['pass'] for m in ('code', 'geo'))
    return out


def raw_variant_id(r):
    return ':'.join(str(r[k]) for k in ('#CHROM', 'POS', 'REF', 'ALT'))


def edit_key(r, grouped_rows=None):
    ref, alt, pos = r['REF'], r['ALT'], int(r['POS'])
    if ',' in ref or ',' in alt:
        if grouped_rows is None:
            raise ValueError('Comma allele needs original per-base records')
        changes = sorted({(int(v['POS_NT']), v['REF_NT'], v['ALT_NT'])
                          for v in grouped_rows if v['REF_NT'] != v['ALT_NT']})
        return ';'.join(str(p) + ':' + a + '>' + b for p, a, b in changes)
    if len(ref) != len(alt) or not ref:
        raise ValueError('Only equal-length source substitutions are supported')
    return ';'.join(str(pos + i) + ':' + a + '>' + b
                    for i, (a, b) in enumerate(zip(ref, alt)) if a != b)


def normalize(rows, statuses, contract):
    """Exact independent collapse/count algebra, conditional on supplied RF decisions."""
    if isinstance(statuses, dict):
        statuses = [{'record_id': k, **v} for k, v in statuses.items()]
    sm = {s['record_id']: s for s in statuses}
    if len(sm) != len(statuses) or len({r['record_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate status/row record ID')
    if set(sm) != {r['record_id'] for r in rows}:
        raise ValueError('RF status coverage differs from rows')
    members = {v['sample']: v for v in contract['members'].values()}
    if len(members) != len(contract['members']):
        raise ValueError('Ambiguous sample metadata')
    groups = defaultdict(list)
    for row in rows:
        status = sm[row['record_id']]
        if status['VAR_ID'] != raw_variant_id(row):
            raise ValueError('Mask identity mismatch')
        if type(status['frequency_pass']) is not bool or type(status['globally_rejected']) is not bool:
            raise ValueError('Mask booleans must be explicit')
        for field in ('DP', 'CAO'):
            if rational(row[field]) < 0:
                raise ValueError('Negative count/depth')
        if status['frequency_pass'] and not status['globally_rejected']:
            key = (row['Sample'], raw_variant_id(row))
            groups[key].append(row)
    collapsed, denominators = {}, defaultdict(list)
    for key, rs in groups.items():
        first = rs[0]
        if len({int(r['POS_NT']) for r in rs}) != len(rs):
            raise ValueError('Duplicate per-base source row')
        dp = median([rational(r['DP']) for r in rs])
        cao = median([rational(r['CAO']) for r in rs])
        rois = {int(str(r.get('ROI', members[r['Sample']]['ROI'])).split('_')[-1]) for r in rs}
        if len(rois) != 1:
            raise ValueError('Conflicting ROI within sample/variant')
        v = {'sample': key[0], 'ROI': next(iter(rois)), 'VAR_ID': key[1],
             'position': int(first['POS']), 'DP': dp, 'CAO': cao, 'edit_key': edit_key(first, rs),
             'source_record_ids': [r['record_id'] for r in rs]}
        collapsed[key] = v
        denominators[(key[0], v['ROI'], v['position'])].append(v)
    for alleles in denominators.values():
        depths = {v['DP'] for v in alleles}
        if len(depths) != 1:
            raise ValueError('Conflicting collapsed position depths')
        wt = next(iter(depths)) - sum(v['CAO'] for v in alleles)
        if wt < 0:
            raise ValueError('Negative reconstructed reference count')
        for v in alleles:
            v['WT'] = wt
            v['q'] = (v['CAO'] + F(1, 2)) / (wt + F(1, 2))
    if any(v['sample'] not in members or members[v['sample']]['ROI'] != v['ROI'] for v in collapsed.values()):
        raise ValueError('Source sample/ROI not authenticated by metadata contract')
    controls = defaultdict(list)
    for v in collapsed.values():
        if members[v['sample']]['trait'] == 'NC':
            controls[(v['ROI'], v['VAR_ID'])].append(v['q'])
    target = {(int(t['ROI']), t['edit_key']): t for t in contract['targets']}
    if len(target) != len(contract['targets']):
        raise ValueError('Nonunique target sequence key')
    measurements = {}
    for v in collapsed.values():
        t = target.get((v['ROI'], v['edit_key']))
        if t is None or members[v['sample']]['trait'] == 'NC':
            continue
        cc = controls.get((v['ROI'], v['VAR_ID']), [])
        background = median(cc) if cc else F()
        corrected = v['q'] - background
        key = (t['hgvs_nt'], v['sample'])
        if key in measurements:
            raise ValueError('Multiple raw variants map to same sample/target')
        measurements[key] = {'hgvs_nt': key[0], 'Sample': key[1], 'q': v['q'],
            'corrected_q': corrected, 'nc_absence_default_used': not bool(cc),
            'log_q': math.log(float(corrected)) if corrected > 0 else None,
            'status': 'positive' if corrected > 0 else 'nonpositive_after_source_NC_subtraction'}
    # Retain absent target/sample combinations explicitly, rather than silently
    # narrowing the denominator to measurements that survived preprocessing.
    for t in contract['targets']:
        for sample, meta in members.items():
            if meta['trait'] != 'NC' and meta['ROI'] == t['ROI']:
                key = (t['hgvs_nt'], sample)
                measurements.setdefault(key, {'hgvs_nt': key[0], 'Sample': sample,
                    'log_q': None, 'status': 'absent_after_declared_source_filter_or_not_in_source'})
    return {'measurements': measurements, 'collapsed': collapsed, 'denominators': denominators}


def contrast_summaries(measurements, contract):
    members = {v['sample']: v for v in contract['members'].values()}
    byid = defaultdict(dict)
    iterable = measurements.values() if isinstance(measurements, dict) else measurements
    for r in iterable:
        if r['Sample'] in byid[r['hgvs_nt']]:
            raise ValueError('Repeated target/sample')
        byid[r['hgvs_nt']][r['Sample']] = r['log_q']
    output = {'code': [], 'geo': []}
    for mapping in output:
        for target in contract['targets']:
            assays, details = {}, {}
            for trait in ('D', 'R', 'F', 'H', 'M'):
                groups, gd = [], []
                for b in (1, 2, 3, 4):
                    selected = sorted(sample for sample, meta in members.items()
                        if meta['ROI'] == target['ROI'] and meta['trait'] == trait and meta[mapping + '_bio'] == b)
                    recorded = [byid[target['hgvs_nt']].get(sample) for sample in selected]
                    vals = [v for v in recorded if v is not None]
                    groups.append(median(vals) if vals else None)
                    gd.append(dict(samples=selected, values=recorded, technical_finite=len(vals)))
                assays[trait] = groups
                details[trait] = gd
            cs = {}
            for name, (a, b) in {'R': ('R', 'D'), 'T': ('D', 'F'), 'L': ('H', 'R'), 'L_mid': ('M', 'R')}.items():
                x, y = assays[a], assays[b]
                valid = all(v is not None for v in x + y)
                if valid:
                    lo, hi = min(x) - max(y), max(x) - min(y)
                    cs[name] = {'complete': True, 'center': median(x) - median(y),
                        'lower': lo, 'upper': hi, 'resolved_sign': 1 if lo > 0 else -1 if hi < 0 else 0}
                else:
                    cs[name] = dict(complete=False, center=None, lower=None, upper=None, resolved_sign='unknown')
            output[mapping].append({**target, 'contrasts': cs, 'reference_preparations': assays,
                                    'reference_preparation_details': details})
    return output


def assert_close(a, b, location):
    if a is None or b is None:
        if a is not b:
            raise AssertionError(location + ': missingness mismatch')
        return
    if isinstance(a, (bool, str)) or isinstance(b, (bool, str)):
        if a != b:
            raise AssertionError(location + ': label mismatch')
        return
    if not math.isclose(float(a), float(b), rel_tol=NUMERIC_RELATIVE_TOLERANCE,
                        abs_tol=NUMERIC_RELATIVE_TOLERANCE):
        raise AssertionError(location + ': numeric mismatch ' + str(a) + ' != ' + str(b))


def author_fraction(obj):
    # Read only a separately authorized stored result. No author function executes.
    if not isinstance(obj, dict) or set(obj) != {'exact_stored_input_arithmetic', 'float'}:
        raise ValueError('Unexpected author arithmetic schema')
    result = F(obj['exact_stored_input_arithmetic'])
    assert_close(result, obj['float'], 'stored fraction/float')
    return result


def compare_downstream(reference, predictions, evaluation):
    checked = Counter()
    for mapping in ('code', 'geo'):
        for mode in ('heavy', 'mid'):
            rr = reference['cases'][mapping][mode]
            pp = predictions[mapping][mode]
            ee = evaluation['cases'][mapping][mode]
            expected_ids = set(rr['predictions'])
            actual_rows = {p['hgvs_nt']: p for p in pp['predictions']}
            if len(actual_rows) != len(pp['predictions']) or set(actual_rows) != expected_ids:
                raise AssertionError('Prediction coverage/identity mismatch')
            for vid, pred in rr['predictions'].items():
                for model in MODELS:
                    actual = author_fraction(actual_rows[vid]['predictions'][model])
                    if pred[model] != actual:
                        raise AssertionError('Exact prediction mismatch ' + mapping + '/' + mode + '/' + vid + '/' + model)
                    checked['predictions_exact'] += 1
            actual_folds = {f['fold']: f for f in pp['folds']}
            if set(actual_folds) != set(rr['fits']):
                raise AssertionError('Fold coverage mismatch')
            for fold, fit in rr['fits'].items():
                af = actual_folds[fold]
                if fit['status'] == 'predicted':
                    if af['status'] != 'predicted' or fit['training_ids'] != af['training_ids']:
                        raise AssertionError('Training membership or fold status mismatch')
                    if fit['training_sites'] != af['training_sites']:
                        raise AssertionError('Training-site mismatch')
                    if fit['weights'] != [author_fraction(w) for w in af['weights']]:
                        raise AssertionError('Training-site weights mismatch')
                    if af['training_count'] != len(fit['training_ids']) or af['assessment_count'] != sum(p['fold'] == fold for p in actual_rows.values()):
                        raise AssertionError('Predicted-fold counts mismatch')
                    for model in ('RNA_linear', 'RNA_load_linear'):
                        actual = af['fit'][model]
                        exact = tuple(author_fraction(actual[k]) for k in ('intercept', 'rna_slope', 'load_slope'))
                        if fit[model][:3] != exact or fit[model][3] != actual['centered_rank']:
                            raise AssertionError('Coefficient/rank mismatch')
                        checked['fit_coefficients_and_rank_exact'] += 1
                elif fit['status'] == 'fit_failure':
                    if af['status'] != 'failed_training_design':
                        raise AssertionError('Failed-fold status mismatch')
                elif af['status'] != 'no_complete_assessment_rows':
                    raise AssertionError('Empty-assessment fold status mismatch')
            au = {r['hgvs_nt']: r for r in pp['unavailable']}
            eu = {r['hgvs_nt']: r for r in ee['unavailable']}
            if len(au) != len(pp['unavailable']) or len(eu) != len(ee['unavailable']):
                raise AssertionError('Duplicate unavailable records')
            if set(au) != set(rr['unavailable']) or au != eu:
                raise AssertionError('Unavailable target coverage/status differs')
            for record in au.values():
                if record['reason'] != 'incomplete_required_assay_preparations' and not record['reason'].startswith('failed_training_design:'):
                    raise AssertionError('Unrecognized unavailability reason')
                detail = rr['unavailable_details'][record['hgvs_nt']]
                if str(record['site']) != str(detail['site']) or str(record['ROI']).replace('ROI_', '') != str(detail['ROI']):
                    raise AssertionError('Unavailable metadata mismatch')
                if detail['failed_design'] != record['reason'].startswith('failed_training_design:'):
                    raise AssertionError('Missingness versus failed-design mismatch')
                if not detail['failed_design'] and record['missing_contrasts'] != detail['missing_contrasts']:
                    raise AssertionError('Missing-assay identities mismatch')
            for vid, recorded in actual_rows.items():
                source_input = rr['input_centers'][vid]
                assert_close(source_input['R'], recorded['R'], 'assessment RNA input identity')
                assert_close(source_input['L'], recorded['L'], 'assessment load input identity')
                if source_input['fold'] != recorded['fold']:
                    raise AssertionError('Assessment fold identity mismatch')
                fit = rr['fits'][recorded['fold']]
                x, z = rational(recorded['R']), rational(recorded['L'])
                rx = not fit['R_range'][0] <= x <= fit['R_range'][1]
                rz = not fit['L_range'][0] <= z <= fit['L_range'][1]
                flags = {'RNA_PAVA': {'rna_outside': rx},
                         'RNA_linear': {'rna_outside': rx, 'load_outside': rz},
                         'RNA_load_linear': {'rna_outside': rx, 'load_outside': rz}}
                if recorded['range_flags'] != flags:
                    raise AssertionError('Extrapolation flags mismatch')
                checked['range_flag_sets'] += 1
            for key in ('scores',):
                compare_scores(rr[key], ee[key], checked)
            for subset in ('SNV', 'ROI_1', 'ROI_2'):
                compare_scores(rr['sensitivities'][subset], ee['sensitivities'][subset], checked)
            if rr['criterion'] != ee['criterion']:
                raise AssertionError('Operational criterion mismatch')
            if rr['RNA_T_sign_confusion'] != ee['RNA_T_sign_confusion']:
                raise AssertionError('Sign confusion mismatch')
            if rr['fixed_denominators'] != ee['fixed_denominators']:
                raise AssertionError('Fixed metadata denominators mismatch')
            errors = {e['hgvs_nt']: e for e in ee['variant_errors']}
            if len(errors) != len(ee['variant_errors']) or set(errors) != expected_ids:
                raise AssertionError('Error coverage mismatch')
            for vid, error in rr['scores']['variant_errors'].items():
                for model in MODELS:
                    if error[model] != author_fraction(errors[vid]['absolute_errors'][model]):
                        raise AssertionError('Exact error mismatch')
                    checked['variant_errors_exact'] += 1
            checked['cases'] += 1
    if reference['primary_operational_expansion_pass'] != evaluation['primary_operational_expansion_pass']:
        raise AssertionError('Primary rule mismatch')
    return dict(checked)


def compare_contrasts(independent, stored):
    checked = 0
    for mapping in ('code', 'geo'):
        actual = {r['hgvs_nt']: r for r in stored[mapping]}
        expected = {r['hgvs_nt']: r for r in independent[mapping]}
        if len(actual) != len(stored[mapping]) or set(actual) != set(expected):
            raise AssertionError('Contrast target completeness mismatch')
        for vid, row in expected.items():
            for trait in ('R', 'T', 'L', 'L_mid'):
                aa, bb = row['contrasts'][trait], actual[vid]['contrasts'][trait]
                if aa['complete'] != bb['complete']:
                    raise AssertionError('Preparation coverage mismatch')
                for key in ('center', 'lower', 'upper'):
                    assert_close(aa[key], bb[key], mapping + '/' + vid + '/' + trait + '/' + key)
                if aa['complete'] and aa['resolved_sign'] != bb['resolved_sign']:
                    raise AssertionError('Source sign mismatch')
                checked += 1
            if 'assays' in actual[vid]:
                for assay, expected_groups in row['reference_preparations'].items():
                    recorded = actual[vid]['assays'][assay]
                    groups = {g['bio']: g for g in recorded['groups']}
                    if set(groups) != {1, 2, 3, 4}:
                        raise AssertionError('Assay preparation group coverage mismatch')
                    for b in (1, 2, 3, 4):
                        assert_close(expected_groups[b - 1], groups[b]['median'], 'technical-log median')
                        detail = row['reference_preparation_details'][assay][b - 1]
                        ag = groups[b]
                        if len(ag['samples']) != len(set(ag['samples'])) or set(ag['samples']) != set(detail['samples']):
                            raise AssertionError('Preparation source-sample lineage differs')
                        av = dict(zip(ag['samples'], ag['values'], strict=True))
                        for sample, value in zip(detail['samples'], detail['values'], strict=True):
                            assert_close(value, av[sample], 'source preparation log value')
                        if ag['technical_finite'] != detail['technical_finite']:
                            raise AssertionError('Finite technical preparation count mismatch')
                        checked += 1
    return checked


def compare_normalization(independent, stored):
    actual = {(r['hgvs_nt'], r['Sample']): r for r in stored['measurements']}
    expected = independent['measurements']
    if len(actual) != len(stored['measurements']) or set(actual) != set(expected):
        raise AssertionError('Normalized measurement coverage mismatch')
    checked = 0
    for key, r in expected.items():
        s = actual[key]
        if 'q' not in r:
            if 'q' in s or s.get('log_q') is not None or s['status'] != r['status']:
                raise AssertionError('Absent source measurement status mismatch')
            checked += 1
            continue
        for name in ('q', 'corrected_q'):
            if r[name] != F(s[name]):
                raise AssertionError('Exact normalized count ratio mismatch ' + str(key) + '/' + name)
            checked += 1
        assert_close(r['log_q'], s['log_q'], 'positive-log or missingness')
        if r['nc_absence_default_used'] != s['nc_absence_default_used']:
            raise AssertionError('Negative-control absence status mismatch')
        checked += 2
    actual_denominators = {(r['Sample'], int(r['POS'])): r for r in stored['denominator_lineage']}
    if len(actual_denominators) != len(stored['denominator_lineage']):
        raise AssertionError('Duplicate denominator lineage')
    expected_denominators = {(s, p): v for (s, roi, p), v in independent['denominators'].items()}
    if set(actual_denominators) != set(expected_denominators):
        raise AssertionError('Background denominator coverage mismatch')
    for key, alleles in expected_denominators.items():
        a = actual_denominators[key]
        if F(a['DP']) != alleles[0]['DP'] or F(a['WT']) != alleles[0]['WT']:
            raise AssertionError('Exact DP/WT mismatch')
        if F(a['sum_retained_CAO']) != sum(v['CAO'] for v in alleles):
            raise AssertionError('Exact total alternate count mismatch')
        source_alleles = {v['VAR_ID']: F(v['CAO']) for v in a['alleles']}
        if source_alleles != {v['VAR_ID']: v['CAO'] for v in alleles}:
            raise AssertionError('Exact collapsed allele counts mismatch')
        checked += 4
    return checked


def read_tsv(path, statuses=False):
    with Path(path).open('r', newline='', encoding='utf8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if statuses:
        for row in rows:
            for k in ('frequency_pass', 'globally_rejected'):
                if row[k] not in ('TRUE', 'FALSE'):
                    raise ValueError('Unexpected mask boolean')
                row[k] = row[k] == 'TRUE'
    return rows


def validate_contrast_metadata(contrasts, contract):
    targets = {t['hgvs_nt']: t for t in contract['targets']}
    if len(targets) != len(contract['targets']):
        raise ValueError('Duplicated metadata target')
    if len(targets) != 125 or len({r['site'] for r in targets.values()}) != 53:
        raise ValueError('The frozen 125-variant/53-site target universe is required')
    if sum(bool(r['is_SNV']) for r in targets.values()) != 113 or len({r['site'] for r in targets.values() if r['is_SNV']}) != 52:
        raise ValueError('The fixed 113-SNV/52-site sensitivity universe is required')
    if set(contrasts) != {'code', 'geo'}:
        raise ValueError('Both and only the frozen grouping maps are required')
    for mapping in contrasts:
        rows = {r['hgvs_nt']: r for r in contrasts[mapping]}
        if len(rows) != len(contrasts[mapping]) or set(rows) != set(targets):
            raise ValueError('Fixed metadata denominator changed')
        for vid, row in rows.items():
            for k in ('site', 'ROI', 'assessment_fold', 'is_SNV'):
                if row[k] != targets[vid][k]:
                    raise ValueError('Fixed target metadata changed: ' + k)
            for name in ('R', 'T', 'L', 'L_mid'):
                c = row['contrasts'][name]
                if type(c['complete']) is not bool:
                    raise ValueError('Unrecognized completeness flag')
                if c['complete']:
                    lo, mid, hi = [rational(c[k]) for k in ('lower', 'center', 'upper')]
                    if not lo <= mid <= hi:
                        raise ValueError('Invalid declared contrast interval')


def compare_scores(ref, obj, counts):
    if (ref['variants'], ref['sites']) != (obj['variants'], obj['sites']):
        raise AssertionError('Score denominator mismatch')
    for m in MODELS:
        actual = None if obj['site_MAE'][m] is None else author_fraction(obj['site_MAE'][m])
        if ref['site_MAE'][m] != actual:
            raise AssertionError('Exact site-balanced MAE mismatch')
        counts['scores_exact'] += 1
    if set(ref['per_site_errors']) != set(obj['per_site_errors']):
        raise AssertionError('Per-site error coverage mismatch')
    for s, vv in ref['per_site_errors'].items():
        for m in MODELS:
            if vv[m] != author_fraction(obj['per_site_errors'][s][m]):
                raise AssertionError('Exact per-site error mismatch')
            counts['per_site_errors_exact'] += 1


def jsonable(value):
    if isinstance(value, F):
        return {'exact': str(value), 'float': float(value)}
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [jsonable(v) for v in value]
    return value


def run(authorization, out):
    out = Path(out).resolve()
    if out.exists():
        raise FileExistsError('Fresh output required')
    auth_path = Path(authorization).resolve()
    auth = json.loads(auth_path.read_text(encoding='utf8'))
    if auth.get('schema') != 'COMT_independent_numeric_authorization_v1' or auth.get('status') != 'explicitly_authorized':
        raise ValueError('Explicit independent authorization required')
    if auth.get('reference_sha256') != sha(__file__):
        raise ValueError('Independent code not hash pinned')
    role_set = set(auth['inputs'])
    required = {'protocol', 'integration_addendum', 'contract', 'contrasts', 'predictions', 'evaluation'}
    raw = {'rows_1', 'rows_2', 'row_status_1', 'row_status_2'}
    if role_set not in (required, required | raw | {'normalization'}):
        raise ValueError('Unexpected/incomplete input-role schema')
    resolved = {}
    for role, spec in auth['inputs'].items():
        if set(spec) != {'path', 'sha256'}:
            raise ValueError('Invalid pinned file specification')
        path = Path(spec['path']).resolve()
        if out == path or out in path.parents or path in out.parents:
            raise ValueError('Output must not contain/replace an input or its ancestor')
        if sha(path) != spec['sha256']:
            raise ValueError('Input hash mismatch: ' + role)
        resolved[role] = path
    if sha(resolved['protocol']) != PROTOCOL_SHA256 or sha(resolved['integration_addendum']) != INTEGRATION_SHA256:
        raise ValueError('The scientific specification differs from the held pre-outcome version')
    # All input hashes pass before any scientific body is parsed.
    contrasts = json.loads(resolved['contrasts'].read_text(encoding='utf8'))
    contract = json.loads(resolved['contract'].read_text(encoding='utf8'))
    validate_contrast_metadata(contrasts, contract)
    raw_checks = {}
    if raw <= role_set:
        rows, masks = [], []
        for roi in (1, 2):
            rows.extend(read_tsv(resolved['rows_' + str(roi)]))
            masks.extend(read_tsv(resolved['row_status_' + str(roi)], statuses=True))
        normalized = normalize(rows, masks, contract)
        reconstructed_contrasts = contrast_summaries(normalized['measurements'], contract)
        raw_checks['contrast_records'] = compare_contrasts(reconstructed_contrasts, contrasts)
        raw_checks['normalization_fields'] = compare_normalization(normalized,
            json.loads(resolved['normalization'].read_text(encoding='utf8')))
        # The raw reconstruction is checked first; model inputs are the pinned stored
        # finite contrast numbers, as required by the protocol's arithmetic convention.
    reference = downstream(contrasts)
    # Derive predictions before opening author answers.
    predictions = json.loads(resolved['predictions'].read_text(encoding='utf8'))
    evaluation = json.loads(resolved['evaluation'].read_text(encoding='utf8'))
    checked = compare_downstream(reference, predictions, evaluation)
    out.mkdir(parents=True)
    (out / 'reconstruction.json').write_text(json.dumps(jsonable(reference), indent=2) + '\n', encoding='utf8')
    receipt = {'status': 'PASS_INDEPENDENT_DOWNSTREAM_COMPARISON', 'checks': checked,
        'authorization_sha256': sha(auth_path), 'reference_sha256': sha(__file__),
        'input_hashes': {r: sha(p) for r, p in resolved.items()},
        'numerical_tolerance': NUMERIC_RELATIVE_TOLERANCE, 'raw_checks': raw_checks,
        'scientific_thresholds_changed': False,
        'RF_training_independently_reconstructed': False,
        'primary_operational_expansion_pass': reference['primary_operational_expansion_pass']}
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf8')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--authorization', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.authorization, args.out), indent=2))
