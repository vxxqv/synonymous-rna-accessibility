"""Replay exact B, C, D, E and G checks after byte-verified materialization.

The checker implementations are preserved scientific artifacts. This adapter
only supplies relocated paths and new output locations; it does not fit models.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def long(path):
    value = str(path.resolve())
    return Path('\\\\?\\' + value) if os.name == 'nt' and not value.startswith('\\\\?\\') else Path(value)


def main():
    if not __debug__:
        raise RuntimeError('Optimized execution is refused')
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    workspace, out = long(args.workspace), long(args.out)
    if out.exists():
        raise FileExistsError('Fresh output directory required')
    materialization = json.loads((workspace / 'materialization_receipt.json').read_text())
    if materialization['status'] not in {'MATERIALIZED_BYTE_VERIFIED', 'MATERIALIZED_AND_COMT_VERIFIED'}:
        raise ValueError('Successful materialization required')
    base = workspace / 'work/research_cycle8'
    # Recheck original logical bytes before any checker is run.
    for record in materialization['restored']:
        if sha(base / record['path']) != record['sha256']:
            raise ValueError('Materialized scientific bytes changed: ' + record['path'])
    out.mkdir(parents=True)
    b1 = base / 'cprime_results/primary_v1/result.json'
    b = base / 'cprime_followup'
    d = base / 'context_design_audit'
    e = base / 'fas_certificate_audit'
    g = base / 'observation_identification_extension'
    # D2's frozen checker writes beside itself. Copy identical code into a fresh
    # sibling directory, preserving its relative-source resolution convention.
    d2 = base / 'context_transport_sharp_boundary'
    fresh_d2 = d2 / 'portable_boundary_review'
    fresh_d2.mkdir()
    shutil.copyfile(d2 / 'root_review/check.py', fresh_d2 / 'check.py')
    commands = [
        ('B', b / 'independent_review/verify.py', ['--b1', b1, '--result', b / 'results_v1/result.json', '--protocol', b / 'protocol.md', '--derivation', b / 'derivation.md', '--analysis', b / 'analysis.py', '--out', out / 'B'], out / 'B/receipt.json'),
        ('D1', d / 'root_review/check.py', ['--design', d, '--source', b1, '--out', out / 'D1'], out / 'D1/verification.json'),
        ('D2', fresh_d2 / 'check.py', [], fresh_d2 / 'verification.json'),
        ('E', e / 'independent_review/verify.py', ['--source', workspace / 'work/research_cycle6/actual_v1/result.json', '--author-result', e / 'results_v1/result.json', '--robustness', workspace / 'work/research_cycle7/results/robustness.json', '--output', out / 'E.json'], out / 'E.json'),
        ('CG', g / 'independent_review/verify_symbolic.py', ['--source-dir', g, '--dao-gate', base / 'dao_gate', '--out', out / 'CG'], out / 'CG/receipt.json'),
    ]
    results = []
    for label, script, flags, receipt in commands:
        log = out / (label + '.log')
        with log.open('w', encoding='utf-8') as stream:
            process = subprocess.run([sys.executable, '-B', str(script), *map(str, flags)], stdout=stream, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(label + ' failed; preserved log: ' + str(log))
        record = json.loads(receipt.read_text())
        expected = {'B': 'PASS', 'D1': 'PASS_SEPARATE_EXACT_CERTIFICATE_AND_DESIGN_RECONSTRUCTION', 'D2': 'PASS_INDEPENDENT_EXACT_BOUNDARY_CERTIFICATES', 'E': 'passed', 'CG': 'PASS_WITH_CONDITIONING_WORDING_CLARIFICATION'}
        if record['status'] != expected[label]:
            raise ValueError(label + ' did not report scientific PASS')
        results.append({'check': label, 'status': record['status'], 'checker_sha256': sha(script), 'receipt_sha256': sha(receipt)})
        print(json.dumps(results[-1]), flush=True)
    summary = {'status': 'PASS_ALL_FIVE_EXACT_REPLAYS', 'materialization_receipt_sha256': sha(workspace / 'materialization_receipt.json'), 'checks': results, 'adapter_sha256': sha(Path(__file__)), 'scope': 'Independent numerical and symbolic reconstruction on frozen source bytes. No new biological replication, calibration or outcome selection.'}
    (out / 'receipt.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
