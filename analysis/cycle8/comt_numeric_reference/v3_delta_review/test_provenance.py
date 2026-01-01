"""Synthetic stage lineage only. Scientific calls are replaced in memory by sentinels."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2) + '\n', encoding='utf8')


def main(candidate, pilot, out):
    candidate, pilot, out = Path(candidate).resolve(), Path(pilot).resolve(), Path(out).resolve()
    if out.exists():
        raise FileExistsError(out)
    before = {p.name: sha(p) for p in candidate.iterdir() if p.is_file()}
    sys.path.insert(0, str(candidate))
    driver = importlib.import_module('driver')
    calls = []
    driver.norm.normalized = lambda rows, status, contract: calls.append('synthetic_normalized_sentinel') or {'invented': True}
    driver.norm.contrast_summaries = lambda result, contract: calls.append('synthetic_summary_sentinel') or {'invented': True}
    checks = []
    out.mkdir(parents=True)
    with tempfile.TemporaryDirectory(dir=out) as td:
        root = Path(td)
        source_paths = {'protocol': pilot / 'protocol_v2.md',
            'amendment': pilot / 'prospective_reconstruction_amendment.md',
            'integration_addendum': candidate / 'integration_addendum.md',
            'contract': candidate / 'selective_contract.json',
            'normalization_amendment': pilot / 'normalization_contract_amendment_v3.md'}
        legacy_hashes = {k: sha(v) for k, v in source_paths.items() if k != 'normalization_amendment'}
        predecessor = root / 'predecessor.json'
        dump(predecessor, {'schema': 'COMT_A1_code_freeze_v1', 'named_input_hashes': legacy_hashes, 'files': {}})
        source_paths['predecessor_code_manifest'] = predecessor
        new_code = root / 'successor.json'
        code_files = {p.name: sha(p) for p in candidate.glob('*.py')}
        for name in ('apply_source_rf.R', '../proxy_models.py', '../source_helpers.R', '../modeling_utils.r', '../summarization_utils.r'):
            code_files[name] = sha(candidate / name)
        dump(new_code, {'schema': 'COMT_A1_code_freeze_v1', 'files': code_files,
                        'named_input_hashes': {k: sha(v) for k, v in source_paths.items()}})
        source_paths['code_manifest'] = new_code
        old_pins = {**legacy_hashes, 'code_manifest': sha(predecessor)}
        ext, rf = root / 'extract', root / 'rf'
        ext.mkdir(); rf.mkdir()
        for name in driver.STAGE_FILES['extract']:
            (ext / name).write_text('record_id\n' if name.startswith('rows_') else '{}\n')
        for name in driver.STAGE_FILES['rf']:
            (rf / name).write_text('record_id\tVAR_ID\tfrequency_pass\tglobally_rejected\n' if name.startswith('row_status') else 'invented\n')
        em = driver.seal(ext, 'extract', 'invented', old_pins)
        ep = ext / 'sealed_manifest.json'
        rm = driver.seal(rf, 'rf', 'invented', {**old_pins, 'extract_manifest': sha(ep)})
        rp = rf / 'sealed_manifest.json'
        def auth():
            paths = {**source_paths, 'extract_manifest': ep, 'rf_manifest': rp}
            return {'schema': 'COMT_A1_authorization_v1', 'phase': 'normalize',
                    'status': 'approved_frozen_normalization_successor',
                    'historical_reproduction_claim': False, 'outcome_access_authorized': True,
                    'synthetic_fixture_only': True,
                    'inputs': {k: {'path': str(v), 'sha256': sha(v)} for k, v in paths.items()}}
        ap = root / 'authorization.json'
        dump(ap, auth())
        receipt = driver.execute('normalize', ap, root / 'valid')
        assert calls == ['synthetic_normalized_sentinel', 'synthetic_summary_sentinel']
        assert set(driver.COMMON_PINS) <= set(receipt['input_hashes'])
        calls.clear(); checks.append('valid_legacy_extract_RF_successor_seals_reach_only_invented_sentinels')
        for phase in ('extract', 'rf'):
            a = auth(); a['phase'] = phase; dump(ap, a)
            try: driver.execute(phase, ap, root / ('forbidden_' + phase))
            except ValueError: pass
            else: raise AssertionError('Successor allowed repeated extraction/RF')
            assert not calls
            checks.append('refuse_repeated_' + phase)
        a = auth(); a['status'] = 'approved_frozen_preoutcome'; dump(ap, a)
        try: driver.execute('normalize', ap, root / 'old_status')
        except ValueError: pass
        else: raise AssertionError('Legacy status authorized successor')
        checks.append('refuse_legacy_authorization_status')
        for phase, original, target in [('extract', em, ep), ('rf', rm, rp)]:
            for pin in driver.LEGACY_PINS:
                changed = copy.deepcopy(original)
                changed['input_hashes'][pin] = '0' * 64
                dump(target, changed)
                dump(ap, auth())  # re-pin the modified test manifest to isolate semantic pin checking
                destination = root / ('changed_' + phase + '_' + pin)
                try: driver.execute('normalize', ap, destination)
                except ValueError as e: assert 'provenance mismatch' in str(e)
                else: raise AssertionError('Changed predecessor pin accepted')
                assert not destination.exists() and not calls
                checks.append('refuse_' + phase + '_' + pin + '_semantic_mismatch')
            dump(target, original)
        # A changed RF->extraction link is rejected even if the files themselves hash correctly.
        changed = copy.deepcopy(rm); changed['input_hashes']['extract_manifest'] = '0' * 64
        dump(rp, changed); dump(ap, auth())
        try: driver.execute('normalize', ap, root / 'wrong_extraction_link')
        except ValueError as e: assert 'another extraction' in str(e)
        else: raise AssertionError('Cross-extraction RF accepted')
        dump(rp, rm); checks.append('refuse_cross_extraction_RF_link')
        # Every successor stage must carry both added pins and the successor code hash.
        for pin in ('normalization_amendment', 'predecessor_code_manifest'):
            a = auth(); a['inputs'].pop(pin); dump(ap, a)
            try: driver.execute('normalize', ap, root / ('missing_' + pin))
            except ValueError: pass
            else: raise AssertionError('Missing successor pin accepted')
            assert not calls
            checks.append('refuse_missing_' + pin)
        expected = {k: sha(v) for k, v in source_paths.items()}
        oldnorm = root / 'oldnorm'; oldnorm.mkdir()
        for n in driver.STAGE_FILES['normalize']:(oldnorm / n).write_text('{}\n')
        driver.seal(oldnorm, 'normalize', 'invented', old_pins)
        try: driver.stage_manifest(oldnorm / 'sealed_manifest.json', 'normalize', expected)
        except ValueError: pass
        else: raise AssertionError('Legacy normalization accepted as successor')
        checks.append('refuse_old_normalization_as_successor_stage')
        try: driver.execute('normalize', root / 'NONEXISTENT_AUTH', root / 'valid')
        except FileExistsError: pass
        else: raise AssertionError('Overwrite guard lost')
        checks.append('overwrite_guard_before_authorization_read')
        for flag in ('-O', '-OO'):
            run = subprocess.run([sys.executable, flag, str(candidate / 'driver.py'), '--help'], capture_output=True)
            assert run.returncode != 0 and b'Optimized Python' in run.stderr
            checks.append('optimization_refused_' + flag)
    assert before == {p.name: sha(p) for p in candidate.iterdir() if p.is_file()}
    result = {'status': 'PASS_INDEPENDENT_SYNTHETIC_SUCCESSOR_PROVENANCE_GUARDS',
        'checks': checks, 'check_count': len(checks), 'candidate_files_unchanged': True,
        'actual_source_fields_or_outcomes_decoded': False, 'scientific_calls_mocked_in_memory': True,
        'driver_sha256': sha(candidate / 'driver.py'), 'test_sha256': sha(__file__)}
    dump(out / 'receipt.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--pilot',required=True);p.add_argument('--out',required=True)
    a = p.parse_args();main(a.candidate,a.pilot,a.out)
