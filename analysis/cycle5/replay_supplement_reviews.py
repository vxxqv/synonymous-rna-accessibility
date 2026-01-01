"""Portable location adapter for the independently authored, byte-pinned oracles.

No mathematical function, input, weighting rule or result is changed. The adapter
sets only repository/output locations, executes the historical checker and
requires equality with its authenticated result, except for Python version.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path

PINS = {
    'development': {
        'directory': 'results/cycle5/development/supplement_review',
        'source': 'verify_development.py',
        'source_sha256': '4db97846dba4e165a6cb63861ef369492b74dd531c726f8ef4cc65eee72cfd30',
        'result_sha256': '872b0ba211762aa1753117569126861fc4f32b19c1b377b554b3449d8a489cc4',
    },
    'reserved_validation': {
        'directory': 'results/cycle5/reserved_validation/supplement_review/correction_v2',
        'source': 'verify_reserved.py',
        'source_sha256': '43cf8b1412bc571877dc8d63aa199b948e4df458f4a171977f4f56d80ac28311',
        'result_sha256': '7c231a7682ddc04c185bdaf3f9633cb0d6de83606b93092cdf4792c66d1a8bd3',
    },
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is not permitted')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--arm', choices=PINS, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    repo, out = args.repo.resolve(strict=True), args.out.resolve()
    if out.exists() or out == repo or repo in out.parents:
        raise RuntimeError('Use a new output directory outside the repository')
    pin = PINS[args.arm]
    directory = repo / pin['directory']
    source = directory / pin['source']
    original = directory / 'independent_results.json'
    if sha(source) != pin['source_sha256'] or sha(original) != pin['result_sha256']:
        raise RuntimeError('Authenticated oracle or result bytes changed')
    out.mkdir()
    spec = importlib.util.spec_from_file_location('independent_cycle5_oracle', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.REPO = repo
    module.HERE = out
    with contextlib.redirect_stdout(io.StringIO()):
        module.main()
    produced = json.loads((out/'independent_results.json').read_text(encoding='utf-8'))
    expected = json.loads(original.read_text(encoding='utf-8'))
    produced_runtime = produced.pop('runtime')
    expected_runtime = expected.pop('runtime')
    if produced != expected:
        raise RuntimeError('Replayed result differs beyond Python runtime version')
    if sha(source) != pin['source_sha256'] or sha(original) != pin['result_sha256']:
        raise RuntimeError('Authenticated files changed during replay')
    receipt = {
        'status': 'verified_exact_result_except_runtime', 'arm': args.arm,
        'oracle_sha256': pin['source_sha256'], 'reference_result_sha256': pin['result_sha256'],
        'adapter_sha256': sha(Path(__file__)), 'runtime': produced_runtime,
        'reference_runtime': expected_runtime,
        'result_sha256': sha(out/'independent_results.json'),
        'scope': 'Selected saved-output audit; no raw assay reconstruction or new outcome model.',
    }
    (out/'replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))

if __name__ == '__main__':
    main()
