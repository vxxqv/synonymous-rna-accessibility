"""Frozen PERSIST transfer comparison. Real stages require a root execution freeze.

Importing this module reads no source table, computes no sequence feature and runs
no folding. Unit tests use synthetic data. The deposited AUP field is never read.
"""
from collections import Counter, deque
from contextlib import contextmanager
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import tempfile

import numpy as np

ALPHABET = 'ACGU'
STOPS = {'UAA', 'UAG', 'UGA'}
CODONS = tuple(''.join(c) for c in itertools.product(ALPHABET, repeat=3) if ''.join(c) not in STOPS)
COMPOSITION_NAMES = tuple('mono_'+x for x in ALPHABET)+tuple('di_'+a+b for a in ALPHABET for b in ALPHABET)+('barcode_GC',)
CODON_NAMES = tuple('codon_'+c for c in CODONS)
MUP_NAME = 'MUP_ViennaRNA_2_7_2_Turner2004_37C_deposited_full'
MODELS = ('development_median', 'development_mean', 'composition_ridge',
          'codon_composition_ridge', 'nearest_CDS', 'codon_composition_MUP_ridge')
PRIMARY = MODELS[-1]
LAMBDAS = (0.1, 1.0, 10.0, 100.0)
MD_SCALARS = ('backbone_length','backtrack','backtrack_type','betaScale','circ','circ_alpha0',
              'circ_penalty','compute_bpp','cv_fact','dangles','energy_set','gquad','helical_rise',
              'logML','max_bp_span','min_loop_size','nc_fact','noGU','noGUclosure','noLP',
              'nonstandards','oldAliEn','pf_smooth','ribo','salt','saltDPXInit','saltDPXInitFact',
              'saltMLLower','saltMLUpper','sfact','special_hp','temperature','uniq_ML','window_size')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def text_sha(value):
    return hashlib.sha256(value.encode('ascii')).hexdigest()


def write_new(path, payload):
    path = Path(path)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write('\n')


def finite(values, label):
    a = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f'Nonfinite required {label}')
    return a


def normalize_rna(sequence):
    if not isinstance(sequence, str) or not sequence or any(c not in 'ACGTU' for c in sequence):
        raise ValueError('Expected nonempty uppercase unambiguous deposited nucleotide string')
    return sequence.replace('T', 'U')


def snapshot_md(md, RNA):
    return {'scalars':{k:getattr(md, k) for k in MD_SCALARS},
            'pair_matrix':[list(row) for row in md.pair], 'reverse_pair_types':list(md.rtype),
            'alias':[RNA.shortP_getitem(md.alias, i) for i in range(RNA.MAXALPHA+1)],
            'option_string':md.option_string()}


def configured_md(settings):
    import RNA
    if RNA.__version__ != settings['version'] or RNA.pf_float_precision() != settings['pf_float_precision_code']:
        raise ValueError('ViennaRNA version or PF precision differs from frozen configuration')
    if not RNA.params_load_RNA_Turner2004():
        raise RuntimeError('Turner2004 parameter load failed')
    md = RNA.md()
    for k, v in settings['effective_md']['scalars'].items():
        setattr(md, k, v)
    if snapshot_md(md, RNA) != settings['effective_md']:
        raise ValueError('Effective model configuration mismatch')
    return RNA, md


def mup_from_upper_bpp(matrix, n, tolerance=1e-10):
    p = finite(matrix, 'base-pair probability matrix')
    if p.shape != (n+1, n+1) or n < 1:
        raise ValueError('Expected ViennaRNA one-based (N+1) square probability matrix')
    if np.any(np.abs(p[0, :]) > tolerance) or np.any(np.abs(p[:, 0]) > tolerance):
        raise ValueError('Nonzero unused matrix index zero')
    if np.any(np.abs(np.tril(p)) > tolerance):
        raise ValueError('Expected upper triangle only, with zero diagonal')
    upper = p[1:, 1:]
    if np.any(upper < -tolerance) or np.any(upper > 1+tolerance):
        raise ValueError('Probability outside declared numerical tolerance')
    sums = np.sum(upper, axis=0)+np.sum(upper, axis=1)
    if np.any(sums < -tolerance) or np.any(sums > 1+tolerance):
        raise ValueError('Per-nucleotide pairing probability outside tolerance')
    pair_total = math.fsum(float(p[i, j]) for i in range(1, n+1) for j in range(i+1, n+1))
    mup = 1-2*pair_total/n
    alternative = math.fsum(1-float(v) for v in sums)/n
    if abs(mup-alternative) > 1e-12 or not -tolerance <= mup <= 1+tolerance:
        raise ValueError('MUP identity or range check failed')
    return {'value':mup, 'row_formula_value':alternative, 'probability_minimum':float(upper.min()),
            'probability_maximum':float(upper.max()), 'maximum_pairing_sum':float(sums.max()),
            'clipped':False, 'matrix_shape':list(p.shape)}


def fold_one(sequence, settings):
    """One fixed calculation. Called on real rows only inside an authorized stage."""
    RNA, md = configured_md(settings)
    rna = normalize_rna(sequence)
    fc = RNA.fold_compound(rna, md, RNA.OPTION_MFE | RNA.OPTION_PF)
    _, mfe = fc.mfe()
    if not math.isfinite(mfe) or abs(mfe) >= RNA.INF/100:
        raise ValueError('Invalid MFE rescaling input')
    # The installed binding exposes pf_scale as read-only. Use its supported
    # MFE-driven rescaling method once, with no setter, fallback or retry.
    fc.exp_params_rescale(float(mfe))
    scale = float(fc.exp_params.pf_scale)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError('Invalid partition-function scale')
    _, ensemble = fc.pf()
    if not math.isfinite(ensemble) or abs(ensemble) >= RNA.INF/100:
        raise ValueError('Invalid partition-function result')
    report = mup_from_upper_bpp(fc.bpp(symmetric=False), len(rna), settings['probability_tolerance'])
    report.update({'length':len(rna), 'normalized_sequence_sha256':text_sha(rna),
                   'mfe_for_rescaling_only':float(mfe), 'ensemble_energy_diagnostic_only':float(ensemble),
                   'effective_pf_scale':scale, 'rescaling_calls':1})
    return report


def sequence_features(row, mup):
    rna = normalize_rna(row['RNA.sequence'])
    cds = normalize_rna(row['CDS.sequence.longest'])
    barcode = normalize_rna(row['Barcode'])
    if len(rna) < 2 or len(cds) % 3 or not cds.startswith('AUG') or cds[-3:] not in STOPS:
        raise ValueError('Invalid fixed feature input')
    sense = [cds[i:i+3] for i in range(0, len(cds)-3, 3)]
    if not sense or any(c in STOPS for c in sense):
        raise ValueError('Invalid sense CDS codon inventory')
    mono, di, codon = Counter(rna), Counter(rna[i:i+2] for i in range(len(rna)-1)), Counter(sense)
    comp = [mono[c]/len(rna) for c in ALPHABET]+[di[a+b]/(len(rna)-1) for a in ALPHABET for b in ALPHABET]
    comp += [(barcode.count('G')+barcode.count('C'))/len(barcode)]
    result = comp+[codon[c]/len(sense) for c in CODONS]+[float(mup)]
    finite(result, 'sequence predictors')
    return result


def hamming(a, b):
    if len(a) != len(b):
        raise ValueError('CDS Hamming distance requires equal length')
    return sum(x != y for x, y in zip(a, b))


def groups_at_distance(ids, sequences, threshold=10):
    graph = {i:[] for i in ids}
    for i, j in itertools.combinations(ids, 2):
        if hamming(sequences[i], sequences[j]) <= threshold:
            graph[i].append(j)
            graph[j].append(i)
    unseen, groups = set(ids), []
    for start in sorted(ids):
        if start not in unseen:
            continue
        unseen.remove(start)
        group, queue = [start], deque([start])
        while queue:
            for nxt in sorted(graph[queue.popleft()]):
                if nxt in unseen:
                    unseen.remove(nxt)
                    group.append(nxt)
                    queue.append(nxt)
        groups.append(sorted(group))
    return groups


def fit_ridge(x, y, penalty):
    x, y = finite(x, 'ridge features'), finite(y, 'development targets')
    if x.ndim != 2 or y.shape != (len(x),) or len(x) < 1 or penalty not in LAMBDAS:
        raise ValueError('Invalid ridge dimensions or penalty')
    active = np.flatnonzero(np.max(x, axis=0) != np.min(x, axis=0))
    means = np.mean(x[:, active], axis=0)
    scales = np.std(x[:, active], axis=0, ddof=0)
    if np.any(scales <= 0) or not np.all(np.isfinite(scales)):
        raise ValueError('Unrepairable training feature scale')
    z = (x[:, active]-means)/scales
    intercept = float(np.mean(y))
    beta = np.linalg.solve(z.T@z+penalty*np.eye(len(active)), z.T@(y-intercept))
    finite(beta, 'ridge coefficients')
    return {'penalty':penalty, 'active_columns':active.tolist(), 'means':means.tolist(),
            'scales':scales.tolist(), 'coefficients':beta.tolist(), 'intercept':intercept,
            'training_count':len(y), 'feature_count':x.shape[1], 'standardization_ddof':0}


def predict_ridge(model, x):
    x = finite(x, 'prediction features')
    if x.ndim != 2 or x.shape[1] != model['feature_count']:
        raise ValueError('Prediction feature shape mismatch')
    z = (x[:, model['active_columns']]-np.asarray(model['means']))/np.asarray(model['scales'])
    return finite(model['intercept']+z@np.asarray(model['coefficients']), 'predictions')


def select_penalty(losses, development_y):
    if set(losses) != set(LAMBDAS) or any(not math.isfinite(v) or v < 0 for v in losses.values()):
        raise ValueError('Invalid complete CV loss grid')
    y = finite(development_y, 'development targets')
    scale = float(np.max(y)-np.min(y))
    if not math.isfinite(scale):
        raise ValueError('Unrepairable development response range')
    tolerance = 1e-12*scale
    best = min(losses.values())
    tied = [lam for lam in LAMBDAS if (losses[lam] == best if scale == 0 else losses[lam] <= best+tolerance)]
    return max(tied), {'response_range':scale, 'tolerance':tolerance, 'best_mae':best, 'tied_penalties':tied}


def tune_ridge(x, y, ids, groups):
    if sorted(i for g in groups for i in g) != sorted(ids) or len(set(ids)) != len(ids) or len(groups) < 2:
        raise ValueError('Invalid nondegenerate development partition')
    x, y = finite(x, 'development matrix'), finite(y, 'development targets')
    index = {key:i for i, key in enumerate(ids)}
    losses, predictions, folds = {}, {}, {}
    for penalty in LAMBDAS:
        oof = np.full(len(ids), np.nan)
        current = []
        for group in groups:
            test = [index[i] for i in group]
            train = [j for j in range(len(ids)) if j not in test]
            fitted = fit_ridge(x[train], y[train], penalty)
            oof[test] = predict_ridge(fitted, x[test])
            current.append({'held_ids':group, 'training_ids':[ids[i] for i in train], 'preparation':fitted})
        finite(oof, 'out-of-fold predictions')
        losses[penalty] = float(np.mean(np.abs(oof-y)))
        predictions[str(penalty)] = oof.tolist()
        folds[str(penalty)] = current
    selected, tie = select_penalty(losses, y)
    return fit_ridge(x, y, selected), {'selected_penalty':selected, 'mae_by_penalty':{str(k):v for k,v in losses.items()},
                                     'tie':tie, 'oof_predictions':predictions, 'folds':folds}


def fit_six(feature_table, development_y, protocol, cds):
    train, holdout = protocol['development_ids'], protocol['holdout_ids']
    if set(development_y) != set(train) or set(feature_table) != set(train+holdout):
        raise ValueError('Exactly frozen training targets and complete feature rows required')
    y = finite([development_y[i] for i in train], 'development targets')
    x = finite([feature_table[i] for i in train], 'development matrix')
    test = finite([feature_table[i] for i in holdout], 'holdout predictors')
    if x.shape[1] != len(COMPOSITION_NAMES)+len(CODON_NAMES)+1:
        raise ValueError('Frozen feature count differs')
    prediction = {MODELS[0]:np.full(len(holdout), np.median(y)), MODELS[1]:np.full(len(holdout), np.mean(y))}
    fitted, tuning = {}, {}
    for name, columns in [(MODELS[2],len(COMPOSITION_NAMES)),(MODELS[3],len(COMPOSITION_NAMES)+len(CODON_NAMES)),(PRIMARY,x.shape[1])]:
        fitted[name], tuning[name] = tune_ridge(x[:, :columns], y, train, protocol['development_groups'])
        prediction[name] = predict_ridge(fitted[name], test[:, :columns])
    nearest, neighbors = [], {}
    for key in holdout:
        distances = [hamming(cds[key], cds[i]) for i in train]
        minimum = min(distances)
        ties = [i for i, dist in enumerate(distances) if dist == minimum]
        nearest.append(float(np.mean(y[ties])))
        neighbors[key] = {'distance':minimum, 'training_ids':[train[i] for i in ties]}
    prediction[MODELS[4]] = np.asarray(nearest)
    return {'development_ids':train, 'holdout_ids':holdout, 'models':list(MODELS),
            'predictions':{name:finite(prediction[name], 'predictions').tolist() for name in MODELS},
            'fitted_ridges':fitted, 'tuning':tuning, 'nearest_neighbors':neighbors,
            'development_targets':{i:float(development_y[i]) for i in train}}


def read_selected_targets(path, selected_ids, endpoint, observer=None):
    """Decode only requested numeric cells. All other data cells remain bytes."""
    selected = set(selected_ids)
    if len(selected) != len(selected_ids):
        raise ValueError('Duplicate requested IDs')
    found = {}
    with Path(path).open('rb') as stream:
        names = stream.readline().rstrip(b'\r\n').decode('utf-8').split('\t')
        if len(names) != len(set(names)) or names.count('SequenceID') != 1 or names.count(endpoint) != 1:
            raise ValueError('Ambiguous endpoint or identifier header')
        id_col, value_col = names.index('SequenceID'), names.index(endpoint)
        for raw_line in stream:
            cells = raw_line.rstrip(b'\r\n').split(b'\t')
            if len(cells) <= id_col:
                raise ValueError('Malformed identifier row')
            key = cells[id_col].decode('ascii')
            if key not in selected:
                continue
            if key in found:
                raise ValueError(f'Duplicate selected identifier: {key}')
            if len(cells) != len(names):
                raise ValueError(f'Malformed selected row: {key}')
            if observer:
                observer(key, endpoint)
            try:
                value = float(cells[value_col].decode('ascii'))
            except (ValueError, UnicodeDecodeError) as exc:
                raise ValueError(f'Invalid selected target: {key}') from exc
            if not math.isfinite(value):
                raise ValueError(f'Nonfinite selected target: {key}')
            found[key] = value
    if set(found) != selected:
        raise ValueError('Missing frozen target identifiers')
    return {key:found[key] for key in selected_ids}


def ranks(a):
    a = np.asarray(a)
    order = np.argsort(a, kind='stable')
    result = np.empty(len(a), dtype=float)
    j = 0
    while j < len(a):
        k = j+1
        while k < len(a) and a[order[k]] == a[order[j]]:
            k += 1
        result[order[j:k]] = (j+k-1)/2+1
        j = k
    return result


def rank_agreement(a, b):
    x, y = ranks(a), ranks(b)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.dot(x-x.mean(), y-y.mean())/math.sqrt(np.sum((x-x.mean())**2)*np.sum((y-y.mean())**2)))


def gain_decision(maes):
    primary = maes[PRIMARY]
    gains = {name:None if maes[name] == 0 else 1-primary/maes[name] for name in MODELS[:-1]}
    passed = all(maes[name] > 0 and primary <= 0.9*maes[name] for name in MODELS[:-1])
    return {'at_least_ten_percent_better_than_every_baseline':passed, 'fractional_gains':gains,
            'zero_error_baselines':[name for name in MODELS[:-1] if maes[name] == 0]}


def evaluate(prediction_record, targets):
    ids = prediction_record['holdout_ids']
    if set(targets) != set(ids) or len(ids) < 2:
        raise ValueError('Complete frozen holdout targets required')
    y = finite([targets[i] for i in ids], 'holdout targets')
    errors, summaries = {}, {}
    for name in MODELS:
        pred = finite(prediction_record['predictions'][name], 'fixed predictions')
        if pred.shape != y.shape:
            raise ValueError('Prediction shape mismatch')
        errors[name] = np.abs(pred-y)
        summaries[name] = {'mae':float(errors[name].mean()), 'rmse':float(np.sqrt(np.mean((pred-y)**2))),
                           'spearman_rank_agreement':rank_agreement(pred, y)}
    influence = []
    for j, key in enumerate(ids):
        mae = {name:float(np.mean(np.delete(error, j))) for name, error in errors.items()}
        influence.append({'omitted_id':key, 'mae':mae, **gain_decision(mae)})
    ranges = {}
    for name in MODELS[:-1]:
        values = [r['fractional_gains'][name] for r in influence if r['fractional_gains'][name] is not None]
        ranges[name] = {'minimum':min(values) if values else None, 'maximum':max(values) if values else None,
                        'undefined_omissions':len(influence)-len(values)}
    return {'models':summaries, **gain_decision({name:summaries[name]['mae'] for name in MODELS}),
            'rows':[{'SequenceID':key, 'target':float(y[j]), 'predictions':{name:prediction_record['predictions'][name][j] for name in MODELS},
                     'absolute_errors':{name:float(errors[name][j]) for name in MODELS}} for j,key in enumerate(ids)],
            'leave_one_holdout_row_out':influence, 'influence_gain_ranges':ranges,
            'influence_interpretation':'Dependent sensitivity diagnostics, not confidence intervals',
            'estimand':'Deposited cellular degradation coefficient, original numeric units'}


def load_metadata(protocol, root):
    raw = json.loads((root/protocol['metadata_path']).read_text(encoding='utf-8'))
    by_id = {}
    for row in raw:
        if row['SequenceID'] in by_id:
            raise ValueError('Duplicate metadata identifier')
        by_id[row['SequenceID']] = row
    result = {}
    for entry in protocol['cohort']:
        key = entry['SequenceID']
        if key not in by_id:
            raise ValueError('Missing frozen sequence metadata')
        row = by_id[key]
        if any(text_sha(row[field]) != entry[digest] for field,digest in
               [('RNA.sequence','deposited_full_sha256'),('CDS.sequence.longest','CDS_sha256'),('Barcode','barcode_sha256')]):
            raise ValueError(f'Sequence metadata hash mismatch: {key}')
        if text_sha(normalize_rna(row['RNA.sequence'])) != entry['normalized_full_sha256']:
            raise ValueError('Normalized full-sequence hash mismatch')
        result[key] = row
    return result


def verify_execution_freeze(freeze_path, root, protocol_path, stage):
    freeze = json.loads(Path(freeze_path).read_text(encoding='utf-8'))
    if freeze.get('status') != 'authorized_by_root' or stage not in freeze.get('authorized_stages', []):
        raise ValueError('Root execution freeze does not authorize this stage')
    required = [str(Path(protocol_path).resolve().relative_to(root)), str(Path(__file__).resolve().relative_to(root))]
    hashes = freeze.get('input_sha256', {})
    if any(p not in hashes for p in required):
        raise ValueError('Execution freeze omits protocol or implementation')
    for relative, expected in hashes.items():
        if sha256(root/relative) != expected:
            raise ValueError(f'Execution freeze input hash mismatch: {relative}')
    protocol = json.loads(Path(protocol_path).read_text(encoding='utf-8'))
    if platform.python_version() != protocol['software']['Python'] or np.__version__ != protocol['software']['NumPy']:
        raise ValueError('Python or NumPy version differs from the frozen configuration')
    for relative, expected in protocol['source_sha256'].items():
        if sha256(root/relative) != expected:
            raise ValueError(f'Protocol source hash mismatch: {relative}')
    return protocol


def verify_prediction_commit(prediction_path, freeze_path):
    frozen = json.loads(Path(freeze_path).read_text(encoding='utf-8'))
    expected = frozen['prediction_sha256']
    if sha256(prediction_path) != expected:
        raise ValueError('Prediction hash mismatch before holdout access')
    commit, git_path = frozen['commit'], frozen['git_path']
    if not commit or ':' in git_path or git_path.startswith('/'):
        raise ValueError('Invalid prediction commit reference')
    content = subprocess.check_output(['git','-C',frozen['repository'],'show',f'{commit}:{git_path}'])
    if hashlib.sha256(content).hexdigest() != expected:
        raise ValueError('Committed predictions do not match frozen file')
    return frozen


def require_frozen_artifact(path, freeze_path, root):
    relative = str(Path(path).resolve().relative_to(root))
    frozen = json.loads(Path(freeze_path).read_text(encoding='utf-8'))
    expected = frozen.get('input_sha256', {}).get(relative)
    if expected is None or sha256(path) != expected:
        raise ValueError('Generated feature artifact is absent from the root training freeze')


def validate_feature_artifact(data, protocol, rows):
    expected_names = list(COMPOSITION_NAMES+CODON_NAMES+(MUP_NAME,))
    if data['ids'] != protocol['all_ids'] or data['names'] != expected_names or set(data['values']) != set(protocol['all_ids']):
        raise ValueError('Feature artifact cohort or column mismatch')
    for key in protocol['all_ids']:
        values = finite(data['values'][key], 'frozen feature row')
        if values.shape != (len(expected_names),):
            raise ValueError('Feature row has wrong shape')
        report = data['diagnostics'][key]
        rna = normalize_rna(rows[key]['RNA.sequence'])
        if report['normalized_sequence_sha256'] != text_sha(rna) or report['length'] != len(rna):
            raise ValueError('Feature diagnostic sequence mismatch')
        if report['rescaling_calls'] != 1 or report['clipped'] or not math.isfinite(report['value']) or not -1e-10 <= report['value'] <= 1+1e-10:
            raise ValueError('Feature calculation did not satisfy the frozen contract')
        if sequence_features(rows[key], report['value']) != values.tolist():
            raise ValueError('Feature vector differs from its frozen input calculation')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('features','train','validate','_fold_worker'))
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--freeze', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--features')
    parser.add_argument('--predictions')
    parser.add_argument('--prediction-freeze')
    parser.add_argument('--sequence-id')
    args = parser.parse_args()
    root, output = Path(args.root).resolve(), Path(args.output)
    logical_stage = 'features' if args.stage == '_fold_worker' else args.stage
    protocol = verify_execution_freeze(args.freeze, root, args.protocol, logical_stage)
    if args.stage == '_fold_worker':
        rows = load_metadata(protocol, root)
        if args.sequence_id not in rows:
            raise ValueError('Worker requested an unfrozen sequence')
        write_new(output, fold_one(rows[args.sequence_id]['RNA.sequence'], protocol['vienna']))
        return
    output.mkdir(parents=True, exist_ok=False)
    try:
        if args.stage == 'features':
            rows = load_metadata(protocol, root)
            values, diagnostics = {}, {}
            for key in protocol['all_ids']:
                destination = output/(key+'_fold.json')
                command = [sys.executable, str(Path(__file__).resolve()), '_fold_worker', '--protocol',args.protocol,
                           '--freeze',args.freeze,'--root',str(root),'--output',str(destination),'--sequence-id',key]
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                if result.returncode or result.stdout.strip() or result.stderr.strip():
                    write_new(output/(key+'_worker_failure.json'), {'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
                    raise RuntimeError(f'Folding failure or numerical warning for {key}; no retry')
                diagnostics[key] = json.loads(destination.read_text())
                values[key] = sequence_features(rows[key], diagnostics[key]['value'])
            write_new(output/'features.json', {'ids':protocol['all_ids'], 'names':list(COMPOSITION_NAMES+CODON_NAMES+(MUP_NAME,)),
                       'values':values, 'diagnostics':diagnostics, 'protocol_sha256':sha256(args.protocol)})
        elif args.stage == 'train':
            if not args.features:
                raise ValueError('Frozen complete feature artifact required')
            require_frozen_artifact(args.features, args.freeze, root)
            data = json.loads(Path(args.features).read_text())
            if data['ids'] != protocol['all_ids'] or data['protocol_sha256'] != sha256(args.protocol) or data['names'] != list(COMPOSITION_NAMES+CODON_NAMES+(MUP_NAME,)):
                raise ValueError('Feature cohort or protocol mismatch')
            rows = load_metadata(protocol, root)
            validate_feature_artifact(data, protocol, rows)
            y = read_selected_targets(root/protocol['target_table_path'], protocol['development_ids'], protocol['endpoint'])
            fitted = fit_six(data['values'], y, protocol, {i:normalize_rna(r['CDS.sequence.longest']) for i,r in rows.items()})
            fitted.update({'protocol_sha256':sha256(args.protocol),'features_sha256':sha256(args.features),
                           'implementation_sha256':sha256(__file__), 'holdout_targets_read':False})
            write_new(output/'predictions.json', fitted)
        else:
            if not args.predictions or not args.prediction_freeze:
                raise ValueError('Committed prediction freeze required before holdout access')
            verified = verify_prediction_commit(args.predictions, args.prediction_freeze)
            fitted = json.loads(Path(args.predictions).read_text())
            if fitted['protocol_sha256'] != sha256(args.protocol) or fitted['holdout_ids'] != protocol['holdout_ids'] or fitted['implementation_sha256'] != sha256(__file__):
                raise ValueError('Prediction protocol or cohort mismatch')
            y = read_selected_targets(root/protocol['target_table_path'], protocol['holdout_ids'], protocol['endpoint'])
            result = evaluate(fitted, y)
            result.update({'prediction_freeze':verified, 'protocol_sha256':sha256(args.protocol)})
            write_new(output/'validation.json', result)
    except Exception as exc:
        write_new(output/'failure.json', {'stage':args.stage, 'exception':type(exc).__name__, 'message':str(exc),
                                       'cohort_changed':False, 'retry_or_imputation':False})
        raise


if __name__ == '__main__':
    main()
