"""Replay the frozen independent FAS audits offline from explicit archive paths."""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python is not allowed')
if sys.version_info < (3, 10):
    raise RuntimeError('Python 3.10 or newer is required; tested with CPython 3.12.14')
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

MANIFEST_SHA = '2aaf5ceee510b1fb51458e6c733152d1864232ef8ba5062c3791071a1a1982bc'
REVIEW_RELATIVE = Path('analysis/cycle6/result_reviews')

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def native(path):
    value = str(path.resolve())
    if os.name == 'nt' and not value.startswith('\\\\?\\'):
        value = '\\\\?\\UNC\\' + value[2:] if value.startswith('\\\\') else '\\\\?\\' + value
    return Path(value)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def checked(path, item):
    raw = path.read_bytes()
    require(len(raw) == item['bytes'] and sha(raw) == item['sha256'], 'Pinned input mismatch: '+item['path'])
    return raw

def execute(command, cwd):
    process = subprocess.run([sys.executable, '-B', *map(str, command)], cwd=cwd,
                             capture_output=True, text=True, check=False)
    require(process.returncode == 0, 'Independent verifier failed: '+process.stdout+process.stderr)
    return process.stdout.strip()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    repo, out = native(args.repo), native(args.out)
    require(not out.exists(), 'Refusing existing output directory')
    require(not out.is_relative_to(repo), 'Output must be outside the computational archive')
    started = time.perf_counter()
    reviews = repo / REVIEW_RELATIVE
    manifest_raw = (reviews/'replay_manifest.json').read_bytes()
    require(sha(manifest_raw) == MANIFEST_SHA, 'Pinned replay manifest mismatch')
    manifest = json.loads(manifest_raw)
    require(manifest['schema'] == 'fas_independent_offline_replay_v1', 'Wrong replay schema')
    seen = set()
    for item in manifest['files']:
        relative = Path(item['path'])
        require(not relative.is_absolute() and '..' not in relative.parts, 'Unsafe manifest path')
        require(item['path'] not in seen, 'Duplicate manifest input')
        seen.add(item['path'])
        checked(repo/relative, item)
    # Both compressed bytes and decompressed exact result are pinned before any write.
    compressed = (repo/'results/cycle6/full_result_v1.json.gz').read_bytes()
    result_raw = gzip.decompress(compressed)
    require(len(result_raw) == manifest['decompressed_result_bytes'] and
            sha(result_raw) == manifest['decompressed_result_sha256'], 'Decompressed result mismatch')
    out.mkdir(parents=True, exist_ok=False)
    result_path = out/'full_result_v1.json'
    result_path.write_bytes(result_raw)
    numeric_script = reviews/'numerical/verify_counts_offline.py'
    synthetic_stdout = execute([numeric_script, '--mode', 'synthetic', '--out', out/'synthetic.json'], out)
    numerical_stdout = execute([numeric_script, '--mode', 'observed', '--authorized-observed',
                                '--repo', repo, '--input-directory', repo/'data/cycle6/fas_counts',
                                '--result', result_path, '--out', out/'numerical.json'], out)
    contract = repo/'analysis/cycle6/fas_pilot/source_contract.json'
    contract_hash = sha(contract.read_bytes())
    mathematical_stdout = execute([reviews/'mathematical/verify_result.py', '--result', result_path,
                                   '--result-sha256', manifest['decompressed_result_sha256'],
                                   '--contract', contract, '--contract-sha256', contract_hash,
                                   '--out', out/'mathematical'], out)
    expected_numeric = json.loads((reviews/'numerical/observed_receipt.json').read_bytes())
    actual_numeric = json.loads((out/'numerical.json').read_bytes())
    require(actual_numeric.pop('git_history_checked_in_this_replay') is False, 'Offline provenance flag')
    for key in ('runtime_seconds', 'verifier_sha256'):
        expected_numeric.pop(key)
        actual_numeric.pop(key)
    require(actual_numeric == expected_numeric, 'Numerical audit differs from original saved evidence')
    expected_synthetic = json.loads((reviews/'numerical/synthetic_receipt.json').read_bytes())
    actual_synthetic = json.loads((out/'synthetic.json').read_bytes())
    for key in ('runtime_seconds', 'verifier_sha256'):
        expected_synthetic.pop(key)
        actual_synthetic.pop(key)
    require(actual_synthetic == expected_synthetic, 'Synthetic safeguards differ from original evidence')
    expected_math = json.loads((reviews/'mathematical/verification.json').read_bytes())
    actual_math = json.loads((out/'mathematical/verification.json').read_bytes())
    require(actual_math == expected_math, 'Mathematical audit differs from original saved evidence')
    for item in manifest['files']:
        checked(repo/item['path'], item)
    require((reviews/'replay_manifest.json').read_bytes() == manifest_raw, 'Manifest changed during replay')
    require(result_path.read_bytes() == result_raw, 'Decompressed result changed during replay')
    receipt = {
        'status':'PASS', 'schema':'fas_independent_offline_replay_receipt_v1',
        'manifest_sha256':MANIFEST_SHA, 'wrapper_sha256':sha(Path(__file__).read_bytes()),
        'verified_archived_inputs':len(manifest['files']),
        'decompressed_result_sha256':sha(result_raw),
        'original_numerical_evidence_equal_except_runtime_and_adapter_hash':True,
        'original_mathematical_evidence_equal_as_complete_json':True,
        'selected_field_scalar_comparisons':actual_numeric['selected_field_scalar_comparisons'],
        'exact_numerical_discrepancy':actual_numeric['maximum_exact_numerical_discrepancy'],
        'certificates_checked':actual_math['certificate_count'],
        'pair_records_checked':actual_math['pair_records_checked'],
        'synthetic_checks':len(actual_synthetic['checks']),
        'git_history_checked':False, 'network_used':False,
        'runtime':sys.version, 'runtime_seconds':time.perf_counter()-started,
        'stdout':{'synthetic':synthetic_stdout,'numerical':numerical_stdout,'mathematical':mathematical_stdout},
        'outputs':{str(path.relative_to(out)).replace('\\','/'):sha(path.read_bytes()) for path in
                   [out/'synthetic.json', out/'numerical.json', out/'mathematical/verification.json']},
        'scope':'Exact offline replay of the fixed saved-data audits, not a new validation experiment or Git-history authentication.'}
    (out/'replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:receipt[key] for key in ['status','verified_archived_inputs','certificates_checked','pair_records_checked','exact_numerical_discrepancy']}))

if __name__ == '__main__':
    main()
