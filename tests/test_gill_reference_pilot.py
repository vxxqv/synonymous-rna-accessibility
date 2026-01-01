"""Synthetic-only tests. No deposited workbook or observed count is opened."""
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
from xml.sax.saxutils import escape
from zipfile import ZipFile

import numpy as np
import pytest


SOURCE = Path(__file__).resolve().parents[1]/'analysis/programme/retrospective/gill_reference_pilot.py'
SPEC = importlib.util.spec_from_file_location('gill_reference_pilot_test_target', SOURCE)
pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)


def toy_rows():
    return [{'A': key, 'B': 0 if i == 0 else 2*(i+7), 'C': 3*(i+11),
             'D': 5*i+37, 'E': 20+i, 'F': 40+2*i, 'G': 60+3*i}
            for i, key in enumerate(pilot.synthetic_keys())]


def toy_workbook(path, rows=None, poison_d=False, shared_keys=False):
    """Minimal XLSX XML with optionally uninterpretable held-out cells."""
    rows = toy_rows() if rows is None else rows
    shared = [str(row['A']) for row in rows]
    xml_rows = []
    for index, row in enumerate([pilot.HEADERS, *rows], 1):
        cells = []
        for col, value in row.items():
            ref = col+str(index)
            if col == 'D' and poison_d:
                cell = f'<c r="{ref}" t="s"><v>9999999</v></c>'
            elif col == 'A' and index > 1 and shared_keys:
                cell = f'<c r="{ref}" t="s"><v>{index-2}</v></c>'
            elif index == 1 or col == 'A':
                cell = f'<c r="{ref}" t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'
            elif isinstance(value, tuple) and value[0] == 'formula':
                cell = f'<c r="{ref}"><f>{escape(value[1])}</f><v>0</v></c>'
            else:
                cell = f'<c r="{ref}"><v>{escape(str(value))}</v></c>'
            cells.append(cell)
        xml_rows.append(f'<row r="{index}">'+''.join(cells)+'</row>')
    with ZipFile(path, 'w') as archive:
        archive.writestr('xl/worksheets/sheet1.xml', f'<worksheet xmlns="{pilot.NS["x"]}"><sheetData>'+''.join(xml_rows)+'</sheetData></worksheet>')
        if shared_keys:
            archive.writestr('xl/sharedStrings.xml', f'<sst xmlns="{pilot.NS["x"]}">'+''.join('<si><t>'+escape(s)+'</t></si>' for s in shared)+'</sst>')
    return path


def toy_frozen_inputs(root, rows=None, poison_d=False):
    root.mkdir(parents=True)
    workbook = toy_workbook(root/'toy.xlsx', rows, poison_d)
    workbook_hash = pilot.digest(workbook)
    eligibility = root/'eligibility.json'
    pilot.write_json(eligibility, {
        'status': 'deposited_pipeline_forward_orientation_and_reference_fibre_design_eligible',
        'inputs': {'pMG23.dna': {'sha256': 'a'*64}, 'cds_processed.xlsx': {'sha256': workbook_hash}},
        'reference_pilot': {
            'vertices': 32, 'edges': 112, 'components': [32],
            'factor_graph': 'K2 Cartesian K4 Cartesian K4', 'regular_degree': 7, 'diameter': 3,
            'eightmer_nt_position_cardinalities': [2, 1, 1, 4, 1, 1, 4, 1],
            'variable_eightmer_positions_one_based': [1, 4, 7],
            'nucleotide_additive_design_rank_including_intercept': 8,
            'codon_additive_design_rank_including_intercept': 8,
            'two_factor_interaction_dimension': 15, 'three_factor_interaction_dimension': 9,
            'all_32_pass_strict_DNA_gt50': True, 'all_32_positive_in_each_DNA_library': True}})
    protocol = root/'protocol.json'
    pilot.write_json(protocol, {
        'status': 'frozen', 'protocol_version': 'gill_reference_pilot_v1',
        'pseudocounts': [0.25, 0.5, 1.0], 'primary_pseudocount': 0.5,
        'training_RNA_columns': ['B', 'C'], 'training_DNA_columns': ['E', 'F'],
        'validation_RNA_column': 'D', 'validation_DNA_column': 'G',
        'target': 'pairwise_half', 'baseline_list': list(pilot.BASELINES),
        'pairwise_shrinkage': 0.5, 'heat_time': 'log(2)/4 with unnormalized Laplacian',
        'keys': pilot.synthetic_keys(), 'reference_key': pilot.synthetic_keys()[3],
        'eligibility_sha256': pilot.digest(eligibility), 'construct_sha256': 'a'*64,
        'workbook_sha256': workbook_hash, 'interpretation': 'Synthetic test only'})
    return protocol, eligibility, workbook


def test_projectors_heat_and_graph():
    design = pilot.design(pilot.synthetic_keys())
    assert design['dimensions'] == [1, 7, 15, 9]
    assert design['ranks'] == {'additive': 8, 'pairwise_model': 23, 'saturated': 32}
    assert sum(map(len, design['graph']))//2 == 112
    assert design['traces_before_anchor']['pairwise_half'] == 15.5
    assert design['traces_after_anchor']['pairwise_half'] == 14.5
    adjacency = np.zeros((32, 32))
    for i, neighbors in enumerate(design['graph']):
        adjacency[i, list(neighbors)] = 1
    values, vectors = np.linalg.eigh(7*np.eye(32)-adjacency)
    independent_heat = (vectors*np.exp(-math.log(2)/4*values))@vectors.T
    np.testing.assert_allclose(design['heat'], independent_heat, atol=2e-15, rtol=2e-14)
    assert np.linalg.matrix_rank(design['projectors'][2]) == 15
    assert np.linalg.matrix_rank(design['projectors'][3]) == 9


def test_weighted_pava_and_exact_ties():
    fit = pilot.isotonic_fit([0, 0, 1, 2], [4, 2, 1, 5])
    assert fit['group_weights'] == [2, 1, 1]
    assert fit['fitted_levels'] == [7/3, 7/3, 5]
    assert fit['constant_blocks'] == 2
    np.testing.assert_allclose(pilot.isotonic_predict(fit, [-1, .5, 1.5, 3]), [7/3, 7/3, 11/3, 5])
    tied = pilot.isotonic_fit([0, np.nextafter(0.0, 1.0), 1], [2, 2, 2])
    assert tied['distinct_inputs'] == 3
    assert tied['constant_blocks'] == 1


@pytest.mark.parametrize('shared_keys', [False, True])
def test_training_never_decodes_d_or_unselected_counts(tmp_path, monkeypatch, shared_keys):
    rows = toy_rows()+[{'A': 'AAAAAAAA', 'B': ('formula', '1/0'), 'C': 'NaN', 'D': 'bad'}]
    workbook = toy_workbook(tmp_path/'poison.xlsx', rows, poison_d=True, shared_keys=shared_keys)
    original = pilot._cell_text
    decoded = []
    def spy(cell, shared):
        ref = cell.attrib['r']
        decoded.append(ref)
        assert not ref.startswith('D')
        assert ref not in ('B34', 'C34')
        return original(cell, shared)
    monkeypatch.setattr(pilot, '_cell_text', spy)
    counts = pilot.selected_counts(workbook, pilot.synthetic_keys(), 'train')
    assert counts['B'][0] == 0
    assert set(counts) == {'B', 'C', 'E', 'F', 'G'}
    assert len(decoded) == 6+32*6+1


@pytest.mark.parametrize('value', ['-1', '0.5', 'NaN', 'Infinity', '-Infinity', '', None, '1e10000'])
def test_invalid_counts_fail(value):
    with pytest.raises(ValueError):
        pilot.count_value(value)


def test_zero_and_exact_integral_counts_retained():
    assert pilot.count_value('0') == 0
    assert pilot.count_value('123.000') == 123
    assert pilot.count_value('1e3') == 1000


@pytest.mark.parametrize('case', ['duplicate', 'missing_key', 'missing_count', 'formula'])
def test_invalid_coverage_fails(tmp_path, case):
    rows = toy_rows()
    if case == 'duplicate':
        rows.append(rows[0])
    elif case == 'missing_key':
        rows.pop()
    elif case == 'missing_count':
        del rows[0]['B']
    else:
        rows[0]['B'] = ('formula', '1+1')
    workbook = toy_workbook(tmp_path/'bad.xlsx', rows)
    with pytest.raises(ValueError):
        pilot.selected_counts(workbook, pilot.synthetic_keys(), 'train')


def test_protocol_requires_frozen_before_workbook_access(tmp_path, monkeypatch):
    protocol, eligibility, workbook = toy_frozen_inputs(tmp_path/'inputs')
    value = json.loads(protocol.read_text())
    value['status'] = 'draft'
    pilot.write_json(protocol, value)
    def forbidden(*args):
        raise AssertionError('No data hashing or reading is allowed with a draft protocol')
    monkeypatch.setattr(pilot, 'digest', forbidden)
    with pytest.raises(ValueError, match='status frozen'):
        pilot.verify_protocol(protocol, eligibility, workbook)


def test_protocol_hash_and_reference_mismatches_fail(tmp_path):
    protocol, eligibility, workbook = toy_frozen_inputs(tmp_path/'inputs')
    pilot.verify_protocol(protocol, eligibility, workbook)
    value = json.loads(protocol.read_text())
    value['reference_key'] = 'AAAAAAAA'
    pilot.write_json(protocol, value)
    with pytest.raises(ValueError, match='reference missing'):
        pilot.verify_protocol(protocol, eligibility, workbook)
    value['reference_key'] = pilot.synthetic_keys()[3]
    value['workbook_sha256'] = '0'*64
    pilot.write_json(protocol, value)
    with pytest.raises(ValueError, match='Workbook source hash'):
        pilot.verify_protocol(protocol, eligibility, workbook)


def test_widest_geometry_against_threshold_connectivity():
    spec = pilot.design(pilot.synthetic_keys())
    field = np.array([((i*13)%19)-9 for i in range(32)], dtype=float)
    geometry = pilot.geometry(field, spec)
    assert geometry['pair_count'] == 496
    for pair in geometry['pairs']:
        expected = -math.inf
        for threshold in sorted(set(field)):
            if min(field[pair['i']], field[pair['j']]) < threshold:
                continue
            seen, todo = {pair['i']}, [pair['i']]
            while todo:
                i = todo.pop()
                for j in spec['graph'][i]:
                    if j not in seen and field[j] >= threshold:
                        seen.add(j)
                        todo.append(j)
            if pair['j'] in seen:
                expected = threshold
        assert pair['capacity'] == expected


def test_validation_hash_guard_precedes_any_workbook_access(tmp_path, monkeypatch):
    prediction = tmp_path/'predictions.json'
    prediction.write_text('{}')
    def forbidden(*args):
        raise AssertionError('Protocol/data must not be accessed before prediction SHA verification')
    monkeypatch.setattr(pilot, 'verify_protocol', forbidden)
    with pytest.raises(ValueError, match='SHA-256 mismatch before validation'):
        pilot.validate('absent_protocol', 'absent_eligibility', 'absent_workbook', prediction, '0'*64, tmp_path/'out')


def test_complete_toy_training_validation_and_heldout_denominator_exclusion(tmp_path, monkeypatch):
    protocol, eligibility, workbook = toy_frozen_inputs(tmp_path/'input')
    train_dir, validate_dir = tmp_path/'train', tmp_path/'validate'
    training_receipt = pilot.train(protocol, eligibility, workbook, train_dir)
    predictions = train_dir/'predictions.json'
    saved = json.loads(predictions.read_text())
    assert saved['exposure']['D_cell_values_accessed'] is False
    assert list(saved['training_counts']) == list(pilot.TRAIN_COLUMNS)
    assert len(json.loads((train_dir/'prediction_table.json').read_text())) == 3*7*32
    assert training_receipt['outputs']['predictions.json'] == pilot.digest(predictions)
    original = pilot._cell_text
    def validate_spy(cell, shared):
        assert cell.attrib['r'][0] in 'ADG'
        return original(cell, shared)
    with monkeypatch.context() as patch:
        patch.setattr(pilot, '_cell_text', validate_spy)
        pilot.validate(protocol, eligibility, workbook, predictions, pilot.digest(predictions), validate_dir)
    result = json.loads((validate_dir/'validation.json').read_text())
    assert set(result['results']) == {'0.25', '0.5', '1.0'}
    assert len(json.loads((validate_dir/'validation_table.json').read_text())) == 3*7*32
    assert len(json.loads((validate_dir/'joint_scenarios.json').read_text())['0.5']) == 9
    for block in result['results'].values():
        assert list(block['models']) == list(pilot.MODELS)
        for metric in block['models'].values():
            assert len(metric['errors']) == 32
            assert metric['errors'][saved['reference_index']] == 0
            assert metric['rmse31'] == pytest.approx(metric['rmse32']*math.sqrt(32/31))
        assert block['observed_geometry']['pair_count'] == 496
    changed = toy_rows()
    for row in changed:
        row['G'] *= 19
    p2, e2, w2 = toy_frozen_inputs(tmp_path/'changed_G', changed)
    pilot.train(p2, e2, w2, tmp_path/'train_changed_G')
    saved2 = json.loads((tmp_path/'train_changed_G/predictions.json').read_text())
    assert saved['predictions'] == saved2['predictions']
    assert saved['identity'] != saved2['identity']


def test_duplicate_rna_and_nonpositive_dna_stop(tmp_path):
    for case in ('duplicate', 'zero_DNA'):
        rows = toy_rows()
        for row in rows:
            if case == 'duplicate':
                row['C'] = row['B']
            else:
                row['G'] = 0
        protocol, eligibility, workbook = toy_frozen_inputs(tmp_path/case, rows)
        with pytest.raises(ValueError, match='duplicate-library|positive DNA'):
            pilot.train(protocol, eligibility, workbook, tmp_path/(case+'_output'))
        assert not (tmp_path/(case+'_output')).exists()


def test_saved_settings_guard_before_d_access(tmp_path, monkeypatch):
    protocol, eligibility, workbook = toy_frozen_inputs(tmp_path/'input', poison_d=True)
    output = tmp_path/'train'
    pilot.train(protocol, eligibility, workbook, output)
    predictions = output/'predictions.json'
    saved = json.loads(predictions.read_text())
    saved['settings']['heat_time'] = 0
    pilot.write_json(predictions, saved)
    def forbidden(*args):
        raise AssertionError('Changed settings must be rejected before selected cell access')
    monkeypatch.setattr(pilot, 'selected_counts', forbidden)
    with pytest.raises(ValueError, match='settings or reference mismatch'):
        pilot.validate(protocol, eligibility, workbook, predictions, pilot.digest(predictions), tmp_path/'validate')


def test_synthetic_controls_open_no_workbook_and_output_guard(tmp_path, monkeypatch):
    def forbidden(*args):
        raise AssertionError('Synthetic mode must never open a workbook')
    monkeypatch.setattr(pilot, 'ZipFile', forbidden)
    output = tmp_path/'synthetic'
    receipt = pilot.synthetic(output)
    result = json.loads((output/'synthetic.json').read_text())
    assert receipt['status'] == 'completed'
    assert result['workbook_opened'] is False
    assert result['positive_decision']['passes_10_percent_hurdle']
    assert not result['negative_decision']['passes_10_percent_hurdle']
    assert result['controls']['valley']['maximum_barrier'] == 2
    assert result['controls']['monotone_cubic']['maximum_barrier'] == 0
    with pytest.raises(FileExistsError):
        pilot.synthetic(output)


def test_synthetic_cli(tmp_path):
    run = subprocess.run([sys.executable, str(SOURCE), 'synthetic', '--outdir', str(tmp_path/'cli')],
                         capture_output=True, text=True, check=True)
    assert json.loads(run.stdout)['mode'] == 'synthetic'


def test_invalid_design_fails():
    keys = pilot.synthetic_keys()
    with pytest.raises(ValueError):
        pilot.design(keys[:-1])
    with pytest.raises(ValueError):
        pilot.design([keys[0], *keys[:-1]])
    with pytest.raises(ValueError):
        pilot.design(list(reversed(keys)))
