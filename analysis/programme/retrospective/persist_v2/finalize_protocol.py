"""Package the metadata-only proposal and synthetic receipt. Never run real stages."""
from pathlib import Path
import json
import subprocess
import sys
import time
from persist_prediction import (sha256, write_new, COMPOSITION_NAMES, CODON_NAMES,
                                MUP_NAME, MODELS, PRIMARY, LAMBDAS)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    targets = ['synthetic_test_receipt.json', 'protocol.json', 'proposed_freeze.json']
    if any((HERE/name).exists() for name in targets):
        raise FileExistsError('Do not overwrite the proposal or test receipt')
    input_names = ['persist_prediction.py', 'test_persist_prediction.py', 'prepare_metadata.py',
                   'metadata_settings.json', 'metadata_settings_initial.json', 'Turner2004.par',
                   'preflight_incidents.json', 'protocol.md', 'finalize_protocol.py']
    before = {str((HERE/name).relative_to(ROOT)):sha256(HERE/name) for name in input_names}
    command = [sys.executable, '-m', 'pytest', str(HERE/'test_persist_prediction.py'), '-q']
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    elapsed = time.perf_counter()-started
    help_result = subprocess.run([sys.executable, str(HERE/'persist_prediction.py'), '--help'],
                                 cwd=ROOT, capture_output=True, text=True, check=False)
    after = {str((HERE/name).relative_to(ROOT)):sha256(HERE/name) for name in input_names}
    receipt = {'status':'pass' if result.returncode == 0 and help_result.returncode == 0 and before == after else 'fail',
               'command':command, 'wall_seconds':elapsed, 'exit_code':result.returncode,
               'stdout':result.stdout, 'stderr':result.stderr, 'input_sha256':before,
               'inputs_unchanged':before == after,
               'help_exit_code':help_result.returncode, 'help_stdout':help_result.stdout,
               'scope':'Synthetic data and one 12-position synthetic fold only. No real features or cellular targets.'}
    write_new(HERE/'synthetic_test_receipt.json', receipt)
    if receipt['status'] != 'pass':
        raise RuntimeError('Synthetic receipt failed; proposal not finalized')
    metadata = json.loads((HERE/'metadata_settings.json').read_text())
    protocol = dict(metadata)
    protocol.update({'protocol_id':'PERSIST_prediction_v2', 'status':'proposed_for_root_review',
                     'endpoint':'degradation_coef_in_cell', 'endpoint_units':'original deposited numeric units; absolute calibration unresolved',
                     'models':list(MODELS), 'primary_model':PRIMARY, 'ridge_lambdas':list(LAMBDAS),
                     'feature_names':list(COMPOSITION_NAMES+CODON_NAMES+(MUP_NAME,)),
                     'composition_columns':len(COMPOSITION_NAMES), 'codon_composition_columns':len(COMPOSITION_NAMES)+len(CODON_NAMES),
                     'primary_columns':len(COMPOSITION_NAMES)+len(CODON_NAMES)+1,
                     'primary_endpoint':'MAE on all 15 fixed holdout rows',
                     'success_rule':'primary MAE <= 0.9 * each baseline MAE, with every baseline MAE > 0',
                     'tuning':'leave one frozen development group out; total absolute held-fold error / 22',
                     'tie_rule':'largest lambda within 1e-12 times full development response range of minimum loss; exact equality if range zero',
                     'standardization':'training-fold mean and population standard deviation; drop exact constant training columns; unpenalized intercept',
                     'ridge_objective':'sum squared residuals + lambda times squared coefficient norm',
                     'secondary':['RMSE','Spearman rank agreement with average ties; undefined if constant'],
                     'influence':'all 15 leave-one-holdout-row-out gain records and ranges; dependent sensitivity, not CI',
                     'coverage':'all37 features, all22 development responses, all15 holdout responses; stop on a required gap, no subset',
                     'excluded_fields':['AUP.Vienna.FULL','all solution stability outcomes','all response-derived predictors'],
                     'authorization':'no real stages authorized by this proposal; root review and execution freeze required',
                     'prediction_commit_gate':True,
                     'reference_sources':['https://viennarna.readthedocs.io/en/latest/api_python.html',
                                          'https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/src/ViennaRNA/model.h']})
    protocol['source_sha256'].update(before)
    protocol['source_sha256'][str((HERE/'synthetic_test_receipt.json').relative_to(ROOT))] = sha256(HERE/'synthetic_test_receipt.json')
    for relative, expected in protocol['source_sha256'].items():
        if sha256(ROOT/relative) != expected:
            raise ValueError(f'Source changed before finalization: {relative}')
    write_new(HERE/'protocol.json', protocol)
    manifest = {'status':'proposed_for_root_review', 'authorized_stages':[],
                'input_sha256':dict(protocol['source_sha256']),
                'real_sequences_folded':0, 'cellular_targets_decoded':0,
                'required_next_action':'Independent source/protocol review, then root commits the exact freeze before authorizing features only.'}
    manifest['input_sha256'][str((HERE/'protocol.json').relative_to(ROOT))] = sha256(HERE/'protocol.json')
    write_new(HERE/'proposed_freeze.json', manifest)
    print(json.dumps({'tests':receipt['status'], 'seconds':elapsed,
                      'protocol_sha256':sha256(HERE/'protocol.json'),
                      'implementation_sha256':sha256(HERE/'persist_prediction.py'),
                      'proposed_freeze_sha256':sha256(HERE/'proposed_freeze.json'),
                      'real_stages_executed':False}))


if __name__ == '__main__':
    main()
