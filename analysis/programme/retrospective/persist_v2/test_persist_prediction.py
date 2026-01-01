"""Synthetic tests only. Never opens the deposited combined outcome table."""
from pathlib import Path
import copy
import json
import math
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import persist_prediction as p


@pytest.fixture
def synthetic():
    rng = np.random.default_rng(9021)
    train = [f'D{i:02d}' for i in range(22)]
    hold = [f'H{i:02d}' for i in range(15)]
    ids = train+hold
    width = len(p.COMPOSITION_NAMES)+len(p.CODON_NAMES)+1
    x = rng.normal(size=(37, width))
    x[:, 0] = 1
    x[:, 2] = x[:, 1]
    y = -2+1.3*x[:22, 3]-0.7*x[:22, -1]+rng.normal(scale=.3, size=22)
    groups = [[key] for key in train[:12]]+[train[12:14], train[14:16], train[16:22]]
    protocol = {'development_ids':train,'holdout_ids':hold,'development_groups':groups}
    return ({key:row.tolist() for key,row in zip(ids,x)},dict(zip(train,y)),protocol,
            {key:format(i,'016b') for i,key in enumerate(ids)})


def test_hand_matrix_and_two_formulas():
    matrix = np.zeros((5,5))
    matrix[1,4], matrix[2,3] = .25, .5
    report = p.mup_from_upper_bpp(matrix,4)
    assert report['value'] == report['row_formula_value'] == .625
    assert not report['clipped']


@pytest.mark.parametrize('kind', ['zero_index','lower','diagonal','over_one','negative','row_sum','nan','shape'])
def test_probability_input_rejections(kind):
    matrix = np.zeros((5,5))
    if kind == 'zero_index': matrix[0,3] = .1
    if kind == 'lower': matrix[3,1] = .1
    if kind == 'diagonal': matrix[2,2] = .1
    if kind == 'over_one': matrix[1,4] = 1.1
    if kind == 'negative': matrix[1,4] = -.1
    if kind == 'row_sum': matrix[1,3],matrix[1,4] = .7,.6
    if kind == 'nan': matrix[1,4] = np.nan
    if kind == 'shape': matrix = matrix[1:,1:]
    with pytest.raises(ValueError):
        p.mup_from_upper_bpp(matrix,4)


def test_probability_boundary_and_no_clipping():
    assert p.mup_from_upper_bpp(np.zeros((2,2)),1)['value'] == 1
    m = np.zeros((3,3)); m[1,2] = 1
    assert p.mup_from_upper_bpp(m,2)['value'] == 0
    m[1,2] = -1e-13
    report = p.mup_from_upper_bpp(m,2)
    assert report['value'] > 1 and not report['clipped']


def test_toy_installed_interface_only():
    settings = json.loads((Path(__file__).parent/'metadata_settings.json').read_text())['vienna']
    report = p.fold_one('GCGCAAAAGCGC',settings)
    assert report['matrix_shape'] == [13,13]
    assert 0 <= report['value'] <= 1 and report['rescaling_calls'] == 1
    assert report['effective_pf_scale'] > 0
    RNA, md = p.configured_md(settings)
    fc = RNA.fold_compound('GCGCAAAAGCGC',md,RNA.OPTION_MFE|RNA.OPTION_PF)
    _, mfe = fc.mfe(); fc.exp_params_rescale(mfe); fc.pf()
    upper = np.asarray(fc.bpp(symmetric=False)); symmetric = np.asarray(fc.bpp(symmetric=True))
    assert np.allclose(symmetric,upper+upper.T,rtol=0,atol=1e-15)
    assert math.isclose(report['value'],1-float(symmetric[1:,1:].sum())/12,rel_tol=0,abs_tol=1e-12)


def test_version_mismatch_stops_before_fold(monkeypatch):
    settings = json.loads((Path(__file__).parent/'metadata_settings.json').read_text())['vienna']
    settings = copy.deepcopy(settings); settings['version'] = 'not-installed'
    with pytest.raises(ValueError,match='version'):
        p.configured_md(settings)


def test_declared_feature_vector_on_toy_string():
    row = {'RNA.sequence':'ATGGCCTAA','CDS.sequence.longest':'ATGGCCTAA','Barcode':'GCA'}
    values = p.sequence_features(row,.75)
    assert len(values) == 83
    assert sum(values[:4]) == pytest.approx(1)
    assert sum(values[4:20]) == pytest.approx(1)
    assert values[20] == pytest.approx(2/3)
    assert sum(values[21:82]) == pytest.approx(1)
    assert values[-1] == .75
    assert values[21+p.CODONS.index('AUG')] == .5
    assert values[21+p.CODONS.index('GCC')] == .5
    assert not any(c in p.CODONS for c in p.STOPS)


def test_frozen_feature_artifact_consistency():
    row = {'RNA.sequence':'ATGGCCTAA','CDS.sequence.longest':'ATGGCCTAA','Barcode':'GCA'}
    report = {'value':.75,'length':9,'normalized_sequence_sha256':p.text_sha('AUGGCCUAA'),
              'rescaling_calls':1,'clipped':False}
    data = {'ids':['D'],'names':list(p.COMPOSITION_NAMES+p.CODON_NAMES+(p.MUP_NAME,)),
            'values':{'D':p.sequence_features(row,.75)},'diagnostics':{'D':report}}
    p.validate_feature_artifact(data,{'all_ids':['D']},{'D':row})
    data['values']['D'][0] += .01
    with pytest.raises(ValueError,match='differs'):
        p.validate_feature_artifact(data,{'all_ids':['D']},{'D':row})


def test_feature_artifact_nonfinite_gate():
    row = {'RNA.sequence':'ATGGCCTAA','CDS.sequence.longest':'ATGGCCTAA','Barcode':'GCA'}
    data = {'ids':['D'],'names':list(p.COMPOSITION_NAMES+p.CODON_NAMES+(p.MUP_NAME,)),
            'values':{'D':[math.nan]*83},'diagnostics':{}}
    with pytest.raises(ValueError,match='Nonfinite'):
        p.validate_feature_artifact(data,{'all_ids':['D']},{'D':row})


@pytest.mark.parametrize('text',['','acgu','ACGN','ACG U','ACG&U'])
def test_alphabet_not_repaired(text):
    with pytest.raises(ValueError): p.normalize_rna(text)


def test_ridge_matches_augmented_least_squares():
    x = np.array([[1.,0.,2.],[1.,1.,4.],[1.,2.,6.],[1.,3.,8.]])
    y = np.array([-2.,0.,4.,7.])
    fitted = p.fit_ridge(x,y,1.)
    active = fitted['active_columns']
    assert active == [1,2]
    z = (x[:,active]-np.mean(x[:,active],axis=0))/np.std(x[:,active],axis=0)
    direct = np.linalg.lstsq(np.vstack([z,np.eye(2)]),np.r_[y-y.mean(),np.zeros(2)],rcond=None)[0]
    assert np.allclose(fitted['coefficients'],direct,rtol=0,atol=1e-12)


def test_homogeneous_tie_rule_and_zero_range():
    losses = {.1:1.,1.:1.+5e-13,10.:2.,100.:3.}
    for a in [1e-6,.001,1.,10.,1e6]:
        selected,_ = p.select_penalty({k:a*v for k,v in losses.items()}, np.array([2.,3.])*a+7.)
        assert selected == 1.
    exact = {.1:0.,1.:0.,10.:0.,100.:0.}
    selected,detail = p.select_penalty(exact,np.array([5.,5.]))
    assert selected == 100. and detail['tolerance'] == 0.
    tiny = {.1:0.,1.:1e-30,10.:1e-20,100.:1e-10}
    assert p.select_penalty(tiny,[5.,5.])[0] == .1


@pytest.mark.parametrize('a,c',[(1e-6,-3.),(.001,7.),(.25,-3.),(1.,7.),(10.,-3.),(1000.,7.),(1e6,-3.)])
def test_complete_positive_affine_equivariance(synthetic,a,c):
    features,y,protocol,cds = synthetic
    original = p.fit_six(features,y,protocol,cds)
    transformed = p.fit_six(features,{i:a*v+c for i,v in y.items()},protocol,cds)
    for model in p.MODELS:
        want = a*np.asarray(original['predictions'][model])+c
        assert np.allclose(transformed['predictions'][model],want,rtol=2e-10,atol=2e-10*max(a,1e-6))
    for model in original['tuning']:
        assert transformed['tuning'][model]['selected_penalty'] == original['tuning'][model]['selected_penalty']


def test_constant_targets_no_positive_gain(synthetic):
    features,y,protocol,cds = synthetic
    result = p.fit_six(features,{i:0. for i in y},protocol,cds)
    for model in p.MODELS: assert result['predictions'][model] == [0.]*15
    for details in result['tuning'].values(): assert details['selected_penalty'] == 100.
    evaluation = p.evaluate(result,{i:0. for i in protocol['holdout_ids']})
    assert not evaluation['at_least_ten_percent_better_than_every_baseline']
    assert len(evaluation['zero_error_baselines']) == 5
    assert all(v['spearman_rank_agreement'] is None for v in evaluation['models'].values())


def test_exact_ten_percent_rule_and_zero_comparator():
    maes = {name:1. for name in p.MODELS}; maes[p.PRIMARY] = .9
    assert p.gain_decision(maes)['at_least_ten_percent_better_than_every_baseline']
    maes[p.MODELS[0]] = 0.; maes[p.PRIMARY] = 0.
    assert not p.gain_decision(maes)['at_least_ten_percent_better_than_every_baseline']
    assert p.gain_decision(maes)['fractional_gains'][p.MODELS[0]] is None


def write_toy_table(path, train, hold, y, hold_payload):
    data = [b'SequenceID\tdegradation_coef_in_cell\tAUP.Vienna.FULL\tother']
    for key in train:
        data.append(key.encode()+b'\t'+str(y[key]).encode()+b'\t\xff\t\xff')
    for key in hold:
        data.append(key.encode()+b'\t'+hold_payload+b'\t\xff\t\xff')
    path.write_bytes(b'\n'.join(data)+b'\n')


def test_selective_training_access_and_holdout_poison(synthetic,tmp_path):
    features,y,protocol,cds = synthetic
    path = tmp_path/'toy.tsv'
    runs = []
    for payload in [b'\xff',b'not-a-number',b'1e300',b'-999']:
        write_toy_table(path,protocol['development_ids'],protocol['holdout_ids'],y,payload)
        observed = []
        selected = p.read_selected_targets(path,protocol['development_ids'],'degradation_coef_in_cell',lambda *v:observed.append(v))
        assert [k for k,_ in observed] == protocol['development_ids']
        runs.append(p.fit_six(features,selected,protocol,cds))
    assert all(run == runs[0] for run in runs[1:])


def test_holdout_features_do_not_set_training_preparation(synthetic):
    features,y,protocol,cds = synthetic
    first = p.fit_six(features,y,protocol,cds)
    altered = copy.deepcopy(features)
    for key in protocol['holdout_ids']: altered[key] = (np.asarray(altered[key])*1e5+3).tolist()
    second = p.fit_six(altered,y,protocol,cds)
    assert first['tuning'] == second['tuning'] and first['fitted_ridges'] == second['fitted_ridges']


def test_fold_preparation_excludes_held_group(synthetic):
    features,y,protocol,_ = synthetic
    ids,group = protocol['development_ids'],protocol['development_groups'][0]
    x = np.asarray([features[i] for i in ids]); target = np.asarray([y[i] for i in ids])
    _,a = p.tune_ridge(x,target,ids,protocol['development_groups'])
    x[0] = x[0]*1e6+100; target[0] += 1e6
    _,b = p.tune_ridge(x,target,ids,protocol['development_groups'])
    for penalty in ['0.1','1.0','10.0','100.0']:
        assert a['folds'][penalty][0]['preparation'] == b['folds'][penalty][0]['preparation']
        assert a['folds'][penalty][0]['held_ids'] == group


@pytest.mark.parametrize('bad',['nan','inf','-inf','missing','duplicate'])
def test_target_coverage_failures(tmp_path,bad):
    path = tmp_path/'toy.tsv'
    if bad == 'missing': data = 'SequenceID\ty\nX\t0\n'
    elif bad == 'duplicate': data = 'SequenceID\ty\nD\t0\nD\t1\n'
    else: data = f'SequenceID\ty\nD\t{bad}\n'
    path.write_text(data)
    with pytest.raises(ValueError): p.read_selected_targets(path,['D'],'y')


def test_negative_and_zero_targets_preserved(tmp_path):
    path = tmp_path/'toy.tsv'; path.write_text('SequenceID\ty\nD\t-3\nE\t0\n')
    assert p.read_selected_targets(path,['D','E'],'y') == {'D':-3.,'E':0.}


def test_extra_holdout_target_in_training_rejected(synthetic):
    features,y,protocol,cds = synthetic
    y = dict(y); y[protocol['holdout_ids'][0]] = 999
    with pytest.raises(ValueError): p.fit_six(features,y,protocol,cds)


def test_prediction_hash_guard_precedes_git_access(tmp_path,monkeypatch):
    pred = tmp_path/'prediction.json'; pred.write_text('{}')
    frozen = tmp_path/'freeze.json'; frozen.write_text(json.dumps({'prediction_sha256':'wrong'}))
    monkeypatch.setattr(p.subprocess,'check_output',lambda *a,**k:pytest.fail('Git should not run after a bad local hash'))
    with pytest.raises(ValueError,match='hash mismatch'): p.verify_prediction_commit(pred,frozen)


def test_prediction_commit_mismatch_rejected(tmp_path,monkeypatch):
    pred = tmp_path/'prediction.json'; pred.write_text('{}')
    frozen = tmp_path/'freeze.json'; frozen.write_text(json.dumps({'prediction_sha256':p.sha256(pred),'commit':'abc','git_path':'predictions.json','repository':'toy'}))
    monkeypatch.setattr(p.subprocess,'check_output',lambda *a,**k:b'different')
    with pytest.raises(ValueError,match='Committed'): p.verify_prediction_commit(pred,frozen)
    monkeypatch.setattr(p.subprocess,'check_output',lambda *a,**k:pred.read_bytes())
    assert p.verify_prediction_commit(pred,frozen)['commit'] == 'abc'


def test_unapproved_execution_freeze_rejected_before_source_access(tmp_path):
    frozen = tmp_path/'freeze.json'; frozen.write_text(json.dumps({'status':'proposed'}))
    with pytest.raises(ValueError,match='does not authorize'):
        p.verify_execution_freeze(frozen,tmp_path,tmp_path/'nonexistent.json','features')


def test_runtime_version_guard_before_source_target_access(tmp_path,monkeypatch):
    source = tmp_path/'toy.py'; source.write_text('pass\n')
    monkeypatch.setattr(p,'__file__',str(source))
    protocol = tmp_path/'protocol.json'
    protocol.write_text(json.dumps({'software':{'Python':'wrong','NumPy':np.__version__},
                                   'source_sha256':{'unread_target.tsv':'irrelevant'}}))
    frozen = tmp_path/'freeze.json'
    frozen.write_text(json.dumps({'status':'authorized_by_root','authorized_stages':['train'],
                                 'input_sha256':{'toy.py':p.sha256(source),'protocol.json':p.sha256(protocol)}}))
    with pytest.raises(ValueError,match='Python or NumPy'):
        p.verify_execution_freeze(frozen,tmp_path,protocol,'train')


def test_near_neighbor_group_transitivity():
    ids = ['a','b','c','d']; sequences = {'a':'0000','b':'0001','c':'0011','d':'1110'}
    assert p.groups_at_distance(ids,sequences,1) == [['a','b','c'],['d']]


def test_ranks_ties_and_influence(synthetic):
    assert p.ranks([3,1,1,2]).tolist() == [4,1.5,1.5,3]
    assert p.rank_agreement([1,2,3],[3,2,1]) == -1
    features,y,protocol,cds = synthetic
    fitted = p.fit_six(features,y,protocol,cds)
    outcome = {i:float(j-5) for j,i in enumerate(protocol['holdout_ids'])}
    result = p.evaluate(fitted,outcome)
    assert len(result['rows']) == len(result['leave_one_holdout_row_out']) == 15
    assert len(result['models']) == 6 and len(result['influence_gain_ranges']) == 5

