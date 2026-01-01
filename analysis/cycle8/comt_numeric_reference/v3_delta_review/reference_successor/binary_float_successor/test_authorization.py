"""Metadata-only provenance tests. No actual assay numeric body is opened."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import copy
import json
import os
from pathlib import Path
import tempfile
import reference_binary_float as ref


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n', encoding='utf8')


def main(pilot, predecessor, old_seal, old_code, new_code, out):
    pilot, predecessor, old_seal, old_code, new_code, out = map(Path, (pilot, predecessor, old_seal, old_code, new_code, out))
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    checks = []
    original_science = ref.scientific_comparison
    calls = []
    def sentinel(*args):
        calls.append('invented_sentinel')
        return {'status': 'INVENTED_SENTINEL_ONLY'}
    ref.scientific_comparison = sentinel
    try:
        temp_parent = Path('\\\\?\\' + str(out.resolve())) if os.name == 'nt' else out.resolve()
        with tempfile.TemporaryDirectory(dir=temp_parent) as temp:
            root = Path(temp)
            assert root.resolve().is_relative_to(temp_parent.resolve())
            paths = {'protocol': pilot / 'protocol_v2.md', 'integration_addendum': pilot / 'implementation_v3/integration_addendum.md',
                'normalization_amendment': pilot / 'normalization_contract_amendment_v3.md',
                'contract': pilot / 'implementation_v3/selective_contract.json', 'predecessor_reference': predecessor,
                'source_first_reference': Path(__file__).resolve().parents[1] / 'reference_source_first.py',
                'source_first_preaccess_seal': Path(__file__).resolve().parents[1] / 'successor_preaccess_seal.json',
                'original_replay_failure_seal': Path(__file__).resolve().parents[1] / 'actual_v3_failure_seal.json',
                'reference_contract_amendment': Path(__file__).with_name('reference_contract_amendment.md'),
                'predecessor_reference_seal': old_seal, 'predecessor_code_manifest': old_code, 'successor_code_manifest': new_code}
            for role in ('contrasts', 'predictions', 'evaluation'):
                paths[role] = root / (role + '.json');dump(paths[role], {'invented': True})
            legacy = {'protocol': ref.PROTOCOL_SHA256, 'amendment': ref.SOURCE_AMENDMENT_SHA256,
                'integration_addendum': ref.INTEGRATION_SHA256, 'contract': ref.sha(paths['contract']),
                'code_manifest': ref.PREDECESSOR_CODE_SHA256}
            current = {**legacy, 'code_manifest': ref.sha(new_code),
                'normalization_amendment': ref.NORMALIZATION_AMENDMENT_SHA256,
                'predecessor_code_manifest': ref.PREDECESSOR_CODE_SHA256}
            specs = [('extract_manifest', 'extract', {}), ('rf_manifest', 'rf', {}),
                ('normalization_manifest', 'normalize', {'contrasts.json': ref.sha(paths['contrasts'])}),
                ('prediction_manifest', 'predict', {'predictions.json': ref.sha(paths['predictions'])}),
                ('evaluation_manifest', 'evaluate', {'evaluation.json': ref.sha(paths['evaluation'])})]
            originals = {}
            for role, phase, files in specs:
                pins = dict(legacy if phase in ('extract', 'rf') else current)
                if phase == 'rf':pins['extract_manifest'] = ref.sha(paths['extract_manifest'])
                if phase == 'normalize':
                    pins['extract_manifest'] = ref.sha(paths['extract_manifest']);pins['rf_manifest'] = ref.sha(paths['rf_manifest'])
                if phase in ('predict', 'evaluate'):pins['normalization_manifest'] = ref.sha(paths['normalization_manifest'])
                if phase == 'evaluate':pins['prediction_manifest'] = ref.sha(paths['prediction_manifest'])
                payload = {'schema':'COMT_A1_stage_seal_v1','phase':phase,'input_hashes':pins,'files':files}
                paths[role] = root / (role + '.json');dump(paths[role],payload);originals[role] = payload
            auth_path = root / 'authorization.json'
            def auth():
                return {'schema':'COMT_independent_binary_float_authorization_v3','status':'explicitly_authorized',
                    'synthetic_fixture_only':True,'reference_sha256':ref.sha(ref.__file__),
                    'inputs':{k:{'path':str(v.resolve()),'sha256':ref.sha(v)} for k,v in paths.items()}}
            dump(auth_path,auth());assert ref.run(auth_path,root/'valid') == {'status':'INVENTED_SENTINEL_ONLY'}
            assert calls == ['invented_sentinel'];calls.clear();checks.append('correct_lineage_reaches_only_mock_science')
            for role, payload in originals.items():
                for key in ('protocol', 'integration_addendum', 'contract', 'code_manifest'):
                    altered = copy.deepcopy(payload);altered['input_hashes'][key] = '0'*64
                    dump(paths[role],altered);dump(auth_path,auth())
                    try:ref.run(auth_path,root/('bad_'+role+'_'+key))
                    except ValueError:pass
                    else:raise AssertionError('Changed stage pin accepted')
                    assert not calls
                    checks.append('reject_'+role+'_'+key)
                dump(paths[role],payload)
            for key in ('normalization_amendment','predecessor_reference','predecessor_reference_seal','predecessor_code_manifest', 'source_first_reference', 'source_first_preaccess_seal', 'original_replay_failure_seal', 'reference_contract_amendment'):
                a = auth();a['inputs'][key]['sha256'] = '0'*64;dump(auth_path,a)
                try:ref.run(auth_path,root/('bad_'+key))
                except ValueError:pass
                else:raise AssertionError('Changed fixed predecessor pin accepted')
                assert not calls;checks.append('reject_changed_'+key)
            # Non-synthetic reconstruction cannot omit any raw role.
            a = auth();a.pop('synthetic_fixture_only');dump(auth_path,a)
            try:ref.run(auth_path,root/'observed_without_raw')
            except ValueError:pass
            else:raise AssertionError('Observed full reconstruction silently downgraded')
            assert not calls;checks.append('observed_reconstruction_requires_all_raw_roles')
            # Failure recording uses invented invalid counts only and never repairs them.
            def invalid(*args):
                raise ref.InvalidDenominatorError({'Sample':'INVENTED','POS':1,'selected_DP':'80','sum_retained_CAO':'106','WT':'-26'})
            ref.scientific_comparison = invalid;dump(auth_path,auth())
            try:ref.run(auth_path,root/'invalid')
            except ref.InvalidDenominatorError:pass
            else:raise AssertionError('Invalid denominator did not stop')
            failure=json.loads((root/'invalid/failure.json').read_text())
            assert failure['status']=='STOP_INVALID_DENOMINATOR_NO_REPAIR' and failure['invalid_denominator']['WT']=='-26'
            assert failure['repair_or_threshold_change_attempted'] is False
            checks.append('invalid_denominator_failure_recorded_no_repair')
            a=auth();dump(auth_path,a);existing=root/'existing';existing.mkdir()
            try:ref.run(root/'NONEXISTENT_AUTH',existing)
            except FileExistsError:pass
            else:raise AssertionError('Overwrite allowed')
            checks.append('existing_output_refused_before_auth_read')
    finally:
        ref.scientific_comparison = original_science
    receipt={'status':'PASS_METADATA_ONLY_SUCCESSOR_AUTHORIZATION_TESTS','checks':checks,'check_count':len(checks),
        'reference_sha256':ref.sha(ref.__file__),'test_sha256':ref.sha(__file__),
        'metadata_code_manifest_hashes':{'predecessor':ref.sha(old_code),'successor':ref.sha(new_code)},
        'actual_scientific_payloads_read':False,'scientific_calls_mocked_in_memory':True}
    dump(out/'receipt.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',required=True);p.add_argument('--predecessor',required=True);p.add_argument('--old-seal',required=True);p.add_argument('--old-code',required=True);p.add_argument('--new-code',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.pilot,a.predecessor,a.old_seal,a.old_code,a.new_code,a.out)
