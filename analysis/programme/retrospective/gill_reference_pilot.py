"""Frozen same-genotype RNA-library prediction pilot with selective XLSX access.

Synthetic mode opens no workbook. Train never accesses a column-D cell value,
including its header. Validate checks saved prediction hashes before opening D.
XML is lexically parsed, but only permitted cells have their values interpreted.
No author code is executed. No confidence interval or significance test is made.
"""
from __future__ import annotations

import argparse
from collections import deque
from decimal import Decimal, InvalidOperation
import hashlib
from heapq import heappop, heappush
from itertools import product
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import numpy as np

NS = {'x': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
TAG = '{' + NS['x'] + '}'
HEADERS = {'A': 'Sequence', 'B': 'CDS_1 Barcode Count', 'C': 'CDS_2 Barcode Count',
           'D': 'CDS_3 Barcode Count', 'E': 'CDS_8 Barcode Count',
           'F': 'CDS_9 Barcode Count', 'G': 'CDS_10 Barcode Count'}
TRAIN_COLUMNS = ('B', 'C', 'E', 'F', 'G')
TEST_COLUMNS = ('D', 'G')
PSEUDOCOUNTS = (0.25, 0.5, 1.0)
BASELINES = ('constant', 'gc', 'additive', 'isotonic_additive', 'heat', 'training_mean')
MODELS = ('constant', 'gc', 'additive', 'isotonic_additive', 'pairwise_half', 'heat', 'training_mean')
HEAT_TIME = math.log(2)/4
REPO = Path(__file__).resolve().parents[3]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def numerical_settings():
    return {'models': list(MODELS), 'baselines': list(BASELINES),
            'pseudocounts': list(PSEUDOCOUNTS), 'primary_pseudocount': 0.5,
            'pairwise_shrinkage': 0.5, 'heat_time': HEAT_TIME,
            'laplacian': 'unnormalized', 'equal_vertex_weights': True,
            'fit_intercept_then_reference_anchor': True,
            'matrix_vector_order': 'ordered Python math.fsum of float64 products',
            'isotonic_ties': 'exact equality of ordered additive predictions, no tolerance',
            'isotonic_interpolation': 'linear between knots, clipped outside',
            'success_fraction': 0.9, 'minimum_baseline_must_be_positive': True,
            'primary_vertices': 32, 'diagnostic_nonreference_vertices': 31,
            'duplicate_RNA_rule': 'stop if selected count vectors are identical; identity unresolved',
            'uncertainty': 'no CI or p-value; fixed finite sensitivity scenarios'}


def design(keys):
    keys = list(keys)
    if len(keys) != 32 or len(set(keys)) != 32 or keys != sorted(keys):
        raise ValueError('Exactly 32 unique keys in frozen lexicographic order required')
    if any(not isinstance(k, str) or len(k) != 8 or set(k)-set('ACGT') for k in keys):
        raise ValueError('Each key must be an eight-base DNA string')
    levels = [sorted({key[p] for key in keys}) for p in range(8)]
    if [len(v) for v in levels] != [2, 1, 1, 4, 1, 1, 4, 1]:
        raise ValueError('Frozen 2 x 4 x 4 design at positions 1,4,7 required')
    if [''.join(v) for v in product(*levels)] != keys:
        raise ValueError('Incomplete Cartesian key coverage')
    dims = (2, 4, 4)
    parts = {}
    for mask in product((0, 1), repeat=3):
        matrix = np.ones((1, 1))
        for bit, size in zip(mask, dims):
            mean = np.ones((size, size))/size
            matrix = np.kron(matrix, np.eye(size)-mean if bit else mean)
        parts[mask] = matrix
    projectors = [sum((p for mask, p in parts.items() if sum(mask) == order), np.zeros((32, 32)))
                  for order in range(4)]
    if [int(np.trace(p)) for p in projectors] != [1, 7, 15, 9]:
        raise AssertionError('Incorrect factorial projector dimensions')
    for i, a in enumerate(projectors):
        for j, b in enumerate(projectors):
            if not np.array_equal(a@b, a if i == j else np.zeros((32, 32))):
                raise AssertionError('Projectors are not exactly orthogonal')
    graph = tuple(tuple(j for j, y in enumerate(keys) if sum(a != b for a, b in zip(x, y)) == 1)
                  for x in keys)
    if sum(map(len, graph)) != 224 or {len(v) for v in graph} != {7}:
        raise AssertionError('Incorrect frozen graph')
    heat = sum((math.exp(-HEAT_TIME*sum(d for bit, d in zip(mask, dims) if bit))*p
                for mask, p in parts.items()), np.zeros((32, 32)))
    return {'keys': keys, 'projectors': projectors, 'additive': projectors[0]+projectors[1],
            'heat': heat, 'graph': graph,
            'gc': np.array([sum(b in 'GC' for b in key) for key in keys], dtype=float),
            'dimensions': [1, 7, 15, 9], 'ranks': {'additive': 8, 'pairwise_model': 23, 'saturated': 32},
            'traces_before_anchor': {'pairwise_half': 15.5, 'heat': float(np.trace(heat))},
            'traces_after_anchor': {'pairwise_half': 14.5, 'heat': float(np.trace(heat))-1}}


def _cell_text(cell, shared):
    """Interpret one explicitly selected cell; callers enforce column access."""
    if cell.find('x:f', NS) is not None:
        raise ValueError('Formula cells are not permitted count or identity inputs')
    typ = cell.attrib.get('t')
    if typ == 'inlineStr':
        node = cell.find('x:is', NS)
        return None if node is None else ''.join(node.itertext())
    node = cell.find('x:v', NS)
    value = None if node is None else node.text
    if typ == 's' and value is not None:
        try:
            return ''.join(shared[int(value)].itertext())
        except (ValueError, IndexError) as error:
            raise ValueError('Invalid shared-string reference') from error
    if typ in ('b', 'e'):
        raise ValueError('Boolean/error cells are not valid counts or keys')
    return value


def count_value(value):
    try:
        number = Decimal(value) if value is not None else Decimal('NaN')
    except (InvalidOperation, TypeError, ValueError) as error:
        raise ValueError('Missing or invalid count') from error
    if not number.is_finite() or number < 0 or number != number.to_integral_value():
        raise ValueError('Counts must be finite nonnegative integers')
    answer = int(number)
    try:
        finite = math.isfinite(float(answer))
    except OverflowError:
        finite = False
    if not finite:
        raise ValueError('Count exceeds the finite float64 logarithm range')
    return answer


def selected_counts(workbook, keys, phase):
    """Decode metadata A and selected count fields only, never other outcomes.

    Train permits B,C,E,F,G. Validate permits D,G. Unselected rows have only A
    interpreted. Shared strings are resolved only when a selected cell needs one.
    """
    if phase not in ('train', 'validate'):
        raise ValueError('Explicit train or validate phase required')
    columns = TRAIN_COLUMNS if phase == 'train' else TEST_COLUMNS
    permitted = {'A', *columns}
    wanted, found = set(keys), {}
    header_seen = False
    with ZipFile(workbook) as archive:
        shared_root = ET.fromstring(archive.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in archive.namelist() else []
        shared = list(shared_root)
        sheets = sorted(name for name in archive.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml', name))
        if len(sheets) != 1:
            raise ValueError('Exactly one worksheet is required by the frozen source mapping')
        with archive.open(sheets[0]) as stream:
            for _, row in ET.iterparse(stream, events=('end',)):
                if row.tag != TAG+'row':
                    continue
                cells = {}
                for cell in row:
                    column = re.sub(r'\d', '', cell.attrib.get('r', ''))
                    if column in permitted:
                        if column in cells:
                            raise ValueError('Duplicate permitted column within a row')
                        cells[column] = cell
                if row.attrib.get('r') == '1':
                    if header_seen:
                        raise ValueError('Repeated header row')
                    header_seen = True
                    observed = {col: _cell_text(cells[col], shared) if col in cells else None for col in permitted}
                    if observed != {col: HEADERS[col] for col in permitted}:
                        raise ValueError('Frozen library header mismatch')
                else:
                    if not header_seen:
                        raise ValueError('Header must precede data rows')
                    key = _cell_text(cells['A'], shared) if 'A' in cells else None
                    if key in wanted:
                        if key in found:
                            raise ValueError('Duplicate selected key')
                        found[key] = {col: count_value(_cell_text(cells[col], shared) if col in cells else None)
                                      for col in columns}
                row.clear()
    if not header_seen or set(found) != wanted:
        raise ValueError('Missing selected keys or header; no vertex exclusions permitted')
    return {col: [found[key][col] for key in keys] for col in columns}


def verify_protocol(protocol_path, eligibility_path, workbook):
    protocol = json.loads(Path(protocol_path).read_text(encoding='utf-8'))
    if protocol.get('status') != 'frozen':
        raise ValueError('Protocol must have status frozen before outcome access')
    expected = {'protocol_version': 'gill_reference_pilot_v1', 'pseudocounts': list(PSEUDOCOUNTS),
                'primary_pseudocount': 0.5, 'training_RNA_columns': ['B', 'C'],
                'training_DNA_columns': ['E', 'F'], 'validation_RNA_column': 'D',
                'validation_DNA_column': 'G', 'target': 'pairwise_half',
                'baseline_list': list(BASELINES), 'pairwise_shrinkage': 0.5,
                'heat_time': 'log(2)/4 with unnormalized Laplacian'}
    for key, value in expected.items():
        if protocol.get(key) != value:
            raise ValueError('Protocol setting mismatch: '+key)
    spec = design(protocol.get('keys', []))
    if protocol.get('reference_key') not in spec['keys']:
        raise ValueError('Deposited reference missing from frozen keys')
    if digest(eligibility_path) != protocol.get('eligibility_sha256'):
        raise ValueError('Eligibility source hash mismatch')
    metadata = json.loads(Path(eligibility_path).read_text(encoding='utf-8'))
    if metadata.get('status') != 'deposited_pipeline_forward_orientation_and_reference_fibre_design_eligible':
        raise ValueError('Eligibility status mismatch')
    pilot = metadata['reference_pilot']
    required = {'vertices': 32, 'edges': 112, 'components': [32], 'factor_graph': 'K2 Cartesian K4 Cartesian K4',
                'regular_degree': 7, 'diameter': 3, 'eightmer_nt_position_cardinalities': [2, 1, 1, 4, 1, 1, 4, 1],
                'variable_eightmer_positions_one_based': [1, 4, 7],
                'nucleotide_additive_design_rank_including_intercept': 8,
                'codon_additive_design_rank_including_intercept': 8,
                'two_factor_interaction_dimension': 15, 'three_factor_interaction_dimension': 9,
                'all_32_pass_strict_DNA_gt50': True, 'all_32_positive_in_each_DNA_library': True}
    if any(pilot.get(key) != value for key, value in required.items()):
        raise ValueError('Eligibility design/coverage mismatch')
    if metadata['inputs']['pMG23.dna']['sha256'] != protocol.get('construct_sha256'):
        raise ValueError('Construct provenance mismatch')
    workbook_hash = digest(workbook)
    if workbook_hash != protocol.get('workbook_sha256') or workbook_hash != metadata['inputs']['cds_processed.xlsx']['sha256']:
        raise ValueError('Workbook source hash mismatch')
    return protocol, spec


def matvec(matrix, vector):
    return np.array([math.fsum(float(a)*float(b) for a, b in zip(row, vector)) for row in matrix])


def anchor(field, reference):
    vector = np.asarray(field, dtype=float)
    if vector.shape != (32,) or not np.isfinite(vector).all():
        raise ValueError('Complete finite 32-vector required')
    return vector-vector[reference]


def isotonic_fit(inputs, response):
    """Exact-float tie groups and deterministic weighted PAVA; no tuning tolerance."""
    inputs = np.asarray(inputs, dtype=float)
    response = np.asarray(response, dtype=float)
    if inputs.shape != response.shape or inputs.ndim != 1 or not np.isfinite(inputs).all() or not np.isfinite(response).all():
        raise ValueError('Finite equally sized isotonic vectors required')
    knots = sorted(set(map(float, inputs)))
    blocks = []
    weights = []
    for index, knot in enumerate(knots):
        values = [float(y) for x, y in zip(inputs, response) if x == knot]
        weights.append(len(values))
        blocks.append([index, index+1, len(values), math.fsum(values)])
        while len(blocks) > 1 and blocks[-2][3]/blocks[-2][2] > blocks[-1][3]/blocks[-1][2]:
            right, left = blocks.pop(), blocks.pop()
            blocks.append([left[0], right[1], left[2]+right[2], math.fsum((left[3], right[3]))])
    fitted = [0.0]*len(knots)
    for begin, end, weight, total in blocks:
        fitted[begin:end] = [total/weight]*(end-begin)
    canonical = 0
    previous = None
    for value in fitted:
        if previous is None or value != previous:
            canonical += 1
        previous = value
    return {'input_knots': knots, 'fitted_levels': fitted, 'group_weights': weights,
            'distinct_inputs': len(knots), 'constant_blocks': canonical}


def isotonic_predict(fit, inputs):
    return np.interp(inputs, fit['input_knots'], fit['fitted_levels'])


def fit_models(training_field, spec, reference):
    z = np.asarray(training_field, dtype=float)
    if z.shape != (32,) or not np.isfinite(z).all():
        raise ValueError('Finite complete training field required')
    mean = math.fsum(map(float, z))/32
    gc = spec['gc']
    centered_gc = gc-math.fsum(map(float, gc))/32
    slope = math.fsum(float(x)*float(y-mean) for x, y in zip(centered_gc, z))/math.fsum(float(x*x) for x in centered_gc)
    additive = matvec(spec['additive'], z)
    isotonic = isotonic_fit(additive, z)
    raw = {'constant': np.full(32, mean), 'gc': mean+slope*centered_gc,
           'additive': additive, 'isotonic_additive': isotonic_predict(isotonic, additive),
           'pairwise_half': additive+0.5*matvec(spec['projectors'][2], z),
           'heat': matvec(spec['heat'], z), 'training_mean': z.copy()}
    return {name: anchor(raw[name], reference).tolist() for name in MODELS}, isotonic


def score_counts(rna, dna, pseudocount):
    rna = np.asarray(rna, dtype=float)
    dna = np.asarray(dna, dtype=float)
    if rna.ndim == 1:
        rna = rna[None, :]
    if dna.ndim == 1:
        dna = dna[None, :]
    return np.array([math.fsum(math.log2(float(rna[k, i])+pseudocount) for k in range(len(rna)))/len(rna)
                     - math.fsum(math.log2(float(dna[k, i])+pseudocount) for k in range(len(dna)))/len(dna)
                     for i in range(32)])


def anova_energy(field, spec):
    return {str(order): math.fsum(float(v*v) for v in matvec(p, field))/32
            for order, p in enumerate(spec['projectors'])}


def floyd_capacities(field, graph):
    z = np.asarray(field, dtype=float)
    capacity = np.full((len(z), len(z)), -np.inf)
    np.fill_diagonal(capacity, z)
    for u, neighbors in enumerate(graph):
        for v in neighbors:
            capacity[u, v] = min(z[u], z[v])
    for via in range(len(z)):
        capacity = np.maximum(capacity, np.minimum(capacity[:, via, None], capacity[None, via, :]))
    return capacity


def widest_capacities(field, graph):
    """Separate heap implementation used to challenge each Floyd result."""
    result = []
    for source in range(len(field)):
        best = [-math.inf]*len(field)
        best[source] = float(field[source])
        heap = [(-best[source], source)]
        while heap:
            neg_value, u = heappop(heap)
            value = -neg_value
            if value != best[u]:
                continue
            for v in graph[u]:
                candidate = min(value, float(field[v]))
                if candidate > best[v]:
                    best[v] = candidate
                    heappush(heap, (-candidate, v))
        result.append(best)
    return np.asarray(result)


def geometry(field, spec):
    z = np.asarray(field, dtype=float)
    if z.shape != (32,) or not np.isfinite(z).all():
        raise ValueError('Geometry requires a finite complete 32-vector')
    capacity = floyd_capacities(z, spec['graph'])
    if not np.array_equal(capacity, widest_capacities(z, spec['graph'])):
        raise AssertionError('Independent widest and Floyd capacities disagree')
    pairs = [{'i': i, 'j': j, 'capacity': float(capacity[i, j]),
              'barrier': float(min(z[i], z[j])-capacity[i, j])}
             for i in range(32) for j in range(i+1, 32)]
    if any(p['barrier'] < 0 for p in pairs):
        raise AssertionError('Negative endpoint-relative barrier')
    thresholds = []
    for threshold in sorted(set(map(float, z)), reverse=True):
        remaining = {i for i in range(32) if z[i] >= threshold}
        retained = len(remaining)
        sizes = []
        while remaining:
            seed = min(remaining)
            remaining.remove(seed)
            queue, size = deque([seed]), 0
            while queue:
                u = queue.popleft()
                size += 1
                for v in spec['graph'][u]:
                    if v in remaining:
                        remaining.remove(v)
                        queue.append(v)
            sizes.append(size)
        thresholds.append({'threshold': threshold, 'retained': retained,
                           'components': len(sizes), 'component_sizes': sorted(sizes, reverse=True)})
    return {'maximum_barrier': max(p['barrier'] for p in pairs), 'pair_count': len(pairs),
            'pairs': pairs, 'critical_thresholds': thresholds, 'anova_mean_square_by_order': anova_energy(z, spec)}


def metrics(prediction, outcome, reference):
    error = np.asarray(prediction, dtype=float)-np.asarray(outcome, dtype=float)
    if error.shape != (32,) or not np.isfinite(error).all() or error[reference] != 0:
        raise ValueError('Finite reference-anchored error vector required')
    mse = math.fsum(float(v*v) for v in error)/32
    signed = math.fsum(map(float, error))/32
    return {'rmse32': math.sqrt(mse),
            'rmse31': math.sqrt(math.fsum(float(error[i]**2) for i in range(32) if i != reference)/31),
            'centered_rmse': math.sqrt(math.fsum(float((v-signed)**2) for v in error)/32),
            'mean_signed_error': signed, 'errors': error.tolist()}


def decision(losses):
    best = min(losses[name]['rmse32'] for name in BASELINES)
    target = losses['pairwise_half']['rmse32']
    return {'best_baseline_rmse32': best, 'target_rmse32': target,
            'best_baseline_names': [name for name in BASELINES if losses[name]['rmse32'] == best],
            'passes_10_percent_hurdle': best > 0 and target <= 0.9*best,
            'strict_improvement': best > 0 and target < best,
            'relative_improvement': None if best == 0 else 1-target/best}


def _new_output(path):
    if Path(path).exists():
        raise FileExistsError('A new output directory is required')


def _identity(protocol_path, eligibility_path, workbook):
    return {'protocol_sha256': digest(protocol_path), 'eligibility_sha256': digest(eligibility_path),
            'workbook_sha256': digest(workbook), 'runner_sha256': digest(__file__)}


def _finish(outdir, payloads, started, identity, mode):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=False)
    for name, payload in payloads.items():
        write_json(outdir/name, payload)
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
        dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True))
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    receipt = {'mode': mode, 'status': 'completed', 'inputs': identity,
               'outputs': {name: digest(outdir/name) for name in payloads},
               'command': sys.argv, 'runtime_seconds': time.perf_counter()-started,
               'repository': {'commit': commit, 'dirty': dirty},
               'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'platform': platform.platform()},
               'scientific_scope': 'fixed 32-genotype measurement prediction, no independent biological inference'}
    write_json(outdir/'receipt.json', receipt)
    return receipt


def train(protocol_path, eligibility_path, workbook, outdir):
    _new_output(outdir)
    started = time.perf_counter()
    protocol, spec = verify_protocol(protocol_path, eligibility_path, workbook)
    before = _identity(protocol_path, eligibility_path, workbook)
    reference = spec['keys'].index(protocol['reference_key'])
    counts = selected_counts(workbook, spec['keys'], 'train')
    if counts['B'] == counts['C']:
        raise ValueError('Identical selected RNA vectors: duplicate-library identity unresolved')
    if any(min(counts[col]) <= 0 for col in ('E', 'F', 'G')):
        raise ValueError('Frozen positive DNA coverage no longer holds')
    if any(sum(counts[col][i] for col in ('E', 'F', 'G')) <= 50 for i in range(32)):
        raise ValueError('Frozen strict combined DNA coverage no longer holds')
    predictions, training_geometry, table = {}, {}, []
    for c in PSEUDOCOUNTS:
        z = score_counts([counts['B'], counts['C']], [counts['E'], counts['F']], c)
        model_predictions, isotonic = fit_models(z, spec, reference)
        predictions[str(c)] = {'training_unanchored': z.tolist(), 'training_anchored': anchor(z, reference).tolist(),
                               'models': model_predictions, 'isotonic': isotonic}
        training_geometry[str(c)] = {name: geometry(field, spec) for name, field in model_predictions.items()}
        table.extend({'pseudocount': c, 'model': name, 'key': key, 'prediction': value}
                     for name, field in model_predictions.items() for key, value in zip(spec['keys'], field))
    if before != _identity(protocol_path, eligibility_path, workbook):
        raise RuntimeError('Input changed during training')
    payload = {'schema': 'gill_reference_predictions_v1', 'identity': before,
               'settings': numerical_settings(), 'keys': spec['keys'], 'reference_key': protocol['reference_key'],
               'reference_index': reference, 'training_counts': counts, 'predictions': predictions,
               'factorial_dimensions': spec['dimensions'], 'ranks': spec['ranks'],
               'traces_before_anchor': spec['traces_before_anchor'], 'traces_after_anchor': spec['traces_after_anchor'],
               'exposure': {'D_cell_values_accessed': False, 'permitted_count_columns': list(TRAIN_COLUMNS),
                            'selected_rows': 32, 'RNA_libraries_used_for_fit': ['B', 'C']}}
    return _finish(outdir, {'predictions.json': payload, 'prediction_table.json': table,
                           'training_geometry.json': training_geometry}, started, before, 'train')


def validate(protocol_path, eligibility_path, workbook, prediction_path, prediction_sha256, outdir):
    _new_output(outdir)
    started = time.perf_counter()
    # This check precedes every workbook hash/open operation in validation.
    if not re.fullmatch(r'[0-9a-f]{64}', prediction_sha256 or '') or digest(prediction_path) != prediction_sha256:
        raise ValueError('Saved prediction SHA-256 mismatch before validation access')
    saved = json.loads(Path(prediction_path).read_text(encoding='utf-8'))
    protocol, spec = verify_protocol(protocol_path, eligibility_path, workbook)
    before = _identity(protocol_path, eligibility_path, workbook)
    if saved.get('schema') != 'gill_reference_predictions_v1' or saved.get('identity') != before:
        raise ValueError('Frozen prediction source identity mismatch')
    if saved.get('settings') != numerical_settings() or saved.get('keys') != spec['keys'] or saved.get('reference_key') != protocol['reference_key']:
        raise ValueError('Frozen prediction settings or reference mismatch')
    reference = spec['keys'].index(protocol['reference_key'])
    if saved.get('reference_index') != reference:
        raise ValueError('Frozen reference index mismatch')
    for c in PSEUDOCOUNTS:
        block = saved['predictions'][str(c)]
        if list(block['models']) != list(MODELS):
            raise ValueError('Incomplete frozen model set')
        for field in block['models'].values():
            vector = np.asarray(field, dtype=float)
            if vector.shape != (32,) or not np.isfinite(vector).all() or vector[reference] != 0:
                raise ValueError('Invalid frozen prediction vector')
    observed = selected_counts(workbook, spec['keys'], 'validate')
    counts = saved['training_counts']
    if observed['G'] != counts['G']:
        raise ValueError('Previously inspected DNA denominator changed')
    if observed['D'] in (counts['B'], counts['C']):
        raise ValueError('Identical selected RNA vectors: duplicate-library identity unresolved')
    counts = {**counts, **observed}
    results, decisions, scenarios, table = {}, {}, {}, []
    for c in PSEUDOCOUNTS:
        z = score_counts(counts['D'], counts['G'], c)
        outcome = anchor(z, reference)
        losses = {name: metrics(saved['predictions'][str(c)]['models'][name], outcome, reference) for name in MODELS}
        decisions[str(c)] = decision(losses)
        results[str(c)] = {'test_unanchored': z.tolist(), 'test_anchored': outcome.tolist(),
                           'models': losses, 'observed_geometry': geometry(outcome, spec)}
        table.extend({'pseudocount': c, 'model': name, 'key': key,
                      'prediction': saved['predictions'][str(c)]['models'][name][i],
                      'outcome': float(outcome[i]), 'error': losses[name]['errors'][i]}
                     for name in MODELS for i, key in enumerate(spec['keys']))
        scenarios[str(c)] = {}
        for rna in ('B', 'C', 'D'):
            for dna in ('E', 'F', 'G'):
                field = anchor(score_counts(counts[rna], counts[dna], c), reference)
                scenarios[str(c)][rna+'/'+dna] = {'field': field.tolist(), 'geometry': geometry(field, spec)}
    if before != _identity(protocol_path, eligibility_path, workbook) or digest(prediction_path) != prediction_sha256:
        raise RuntimeError('Input changed during validation')
    summary = {'primary_pseudocount': 0.5, 'primary_decision': decisions['0.5'],
               'all_pseudocount_decisions': decisions,
               'same_improvement_direction_at_both_sensitivities': all(decisions[str(c)]['strict_improvement'] for c in (0.25, 1.0)),
               'primary_success_and_sensitivity_direction': decisions['0.5']['passes_10_percent_hurdle'] and all(decisions[str(c)]['strict_improvement'] for c in (0.25, 1.0)),
               'joint_scenarios': 'Nine overlapping library contrasts; not independent validation replicates or uncertainty regions',
               'scenario_maximum_barrier_ranges': {c: {'minimum': min(v['geometry']['maximum_barrier'] for v in group.values()),
                                                     'maximum': max(v['geometry']['maximum_barrier'] for v in group.values())}
                                                  for c, group in scenarios.items()},
               'scope': protocol.get('interpretation'), 'uncertainty': 'No CI, p-value or biological interaction inference'}
    payload = {'schema': 'gill_reference_validation_v1', 'identity': before, 'prediction_sha256': prediction_sha256,
               'keys': spec['keys'], 'reference_key': protocol['reference_key'], 'settings': numerical_settings(),
               'results': results, 'summary': summary,
               'validation_counts': observed, 'count_columns_accessed': list(TEST_COLUMNS)}
    return _finish(outdir, {'validation.json': payload, 'validation_table.json': table,
                           'joint_scenarios.json': scenarios}, started,
                   {**before, 'prediction_sha256': prediction_sha256}, 'validate')


def synthetic_keys():
    return [''.join(p) for p in product('CT', 'G', 'T', 'ACGT', 'C', 'C', 'ACGT', 'T')]


def synthetic(outdir):
    _new_output(outdir)
    started = time.perf_counter()
    spec = design(synthetic_keys())
    vertices = list(product(range(2), range(4), range(4)))
    levels = ([-1, 1], [-3, -1, 1, 3], [-3, -1, 1, 3])
    additive = np.array([sum(levels[i][v[i]] for i in range(3)) for v in vertices], dtype=float)
    pair = np.array([levels[0][v[0]]*levels[1][v[1]] for v in vertices], dtype=float)
    triple = np.array([math.prod(levels[i][v[i]] for i in range(3)) for v in vertices], dtype=float)
    if not np.array_equal(matvec(spec['additive'], additive), additive):
        raise AssertionError('Additive control not recovered')
    if not np.array_equal(matvec(spec['projectors'][2], pair), pair) or not np.array_equal(matvec(spec['projectors'][3], triple), triple):
        raise AssertionError('Interaction controls not identified')
    reference = 0
    positive, _ = fit_models(additive+2*pair+triple, spec, reference)
    positive_losses = {name: metrics(field, anchor(additive+pair, reference), reference) for name, field in positive.items()}
    if positive_losses['pairwise_half']['rmse32'] != 0 or not decision(positive_losses)['passes_10_percent_hurdle']:
        raise AssertionError('Known pairwise-shrinkage positive control failed')
    negative, _ = fit_models(additive, spec, reference)
    negative_losses = {name: metrics(field, anchor(additive, reference), reference) for name, field in negative.items()}
    if negative_losses['additive']['rmse32'] != 0 or decision(negative_losses)['passes_10_percent_hurdle']:
        raise AssertionError('Perfect additive baseline must prevent a success claim')
    valley = np.full(32, -1.0)
    valley[vertices.index((0, 0, 0))] = valley[vertices.index((1, 1, 0))] = 1.0
    fields = {'additive': additive, 'monotone_cubic': additive+additive**3,
              'monotone_tanh': np.tanh(additive/3), 'valley': valley}
    controls = {name: geometry(field, spec) for name, field in fields.items()}
    if any(controls[name]['maximum_barrier'] != 0 for name in ('additive', 'monotone_cubic', 'monotone_tanh')) or controls['valley']['maximum_barrier'] != 2:
        raise AssertionError('Known all-threshold barrier controls failed')
    if controls['monotone_cubic']['anova_mean_square_by_order']['2'] <= 0 or controls['monotone_tanh']['anova_mean_square_by_order']['2'] <= 0:
        raise AssertionError('Monotone transformations should expose raw interaction ambiguity')
    output = {'status': 'synthetic controls passed', 'workbook_opened': False, 'observed_counts_read': False,
              'settings': numerical_settings(), 'dimensions': spec['dimensions'], 'ranks': spec['ranks'],
              'traces_before_anchor': spec['traces_before_anchor'], 'traces_after_anchor': spec['traces_after_anchor'],
              'positive_model_errors': positive_losses, 'negative_model_errors': negative_losses,
              'positive_decision': decision(positive_losses), 'negative_decision': decision(negative_losses),
              'controls': controls, 'noise_alias_exact': bool(np.array_equal(np.ones(32)+(valley-1), valley))}
    return _finish(outdir, {'synthetic.json': output}, started, {'runner_sha256': digest(__file__)}, 'synthetic')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    for mode in ('synthetic', 'train', 'validate'):
        command = sub.add_parser(mode)
        command.add_argument('--outdir', type=Path, required=True)
        if mode != 'synthetic':
            command.add_argument('--protocol', type=Path, required=True)
            command.add_argument('--eligibility', type=Path, required=True)
            command.add_argument('--workbook', type=Path, required=True)
        if mode == 'validate':
            command.add_argument('--predictions', type=Path, required=True)
            command.add_argument('--prediction-sha256', required=True)
    args = parser.parse_args()
    if args.mode == 'synthetic':
        receipt = synthetic(args.outdir)
    elif args.mode == 'train':
        receipt = train(args.protocol, args.eligibility, args.workbook, args.outdir)
    else:
        receipt = validate(args.protocol, args.eligibility, args.workbook, args.predictions, args.prediction_sha256, args.outdir)
    print(json.dumps({'mode': receipt['mode'], 'status': receipt['status'], 'outputs': receipt['outputs']}, indent=2))


if __name__ == '__main__':
    main()
