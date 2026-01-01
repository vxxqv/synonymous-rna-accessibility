"""Bounded adapter guards only. No historical scientific main is run."""
if not __debug__:
    raise RuntimeError('Run audit checks without optimization')
import ast
import contextlib
import copy
import hashlib
import importlib.util
import io
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[2]
REPO = WORK / 'outputs/synonymous-rna-accessibility'
SOURCE = REPO / 'analysis/cycle5/replay_supplement_reviews.py'
PIN = '999f4b4ea8fc453daecfd7aa130412c1d321277cc298a85722e8cf4b6bed7aab'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok, text):
    if not ok: raise RuntimeError(text)
require(sha(SOURCE) == PIN, 'Adapter hash mismatch')
spec = importlib.util.spec_from_file_location('adapter_under_review', SOURCE)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
checks = []
def record(name): checks.append(name)
def run_cli(repo, arm, out, flags=()):
    return subprocess.run([sys.executable, *flags, str(SOURCE), '--repo', str(repo), '--arm', arm, '--out', str(out)], capture_output=True, text=True)
def reject(name, repo, arm, out, flags=(), text=None):
    existed = out.exists()
    result = run_cli(repo, arm, out, flags)
    require(result.returncode != 0, name + ' incorrectly passed')
    if text: require(text in result.stderr, name + ' wrong failure: ' + result.stderr)
    require(existed or not out.exists(), name + ' created output')
    record(name)

source_hashes = {str(SOURCE.relative_to(WORK)): PIN}
with tempfile.TemporaryDirectory(prefix='guards_', dir=HERE) as temp:
    # Extended-length path only for this deeply nested Windows test staging.
    base = Path('\\\\?\\' + str(Path(temp).resolve())) if os.name == 'nt' else Path(temp)
    for arm, pin in adapter.PINS.items():
        repo = base / ('archive_' + arm)
        directory = repo / pin['directory']
        directory.mkdir(parents=True)
        original_source = REPO / pin['directory'] / pin['source']
        original_result = REPO / pin['directory'] / 'independent_results.json'
        source_hashes[str(original_source.relative_to(WORK))] = sha(original_source)
        source_hashes[str(original_result.relative_to(WORK))] = sha(original_result)
        require(sha(original_source) == pin['source_sha256'], arm+' source pin')
        require(sha(original_result) == pin['result_sha256'], arm+' result pin')
        local_source = directory / pin['source']
        local_result = directory / 'independent_results.json'
        local_source.write_bytes(original_source.read_bytes())
        local_result.write_bytes(original_result.read_bytes())
        # An empty archive-like root is deliberately enough for guard tests.
        # It has neither Git nor the scientific inputs, so scientific replay
        # cannot accidentally run during these tests.
        reject(arm+': inside checkout', repo, arm, repo/'new_output', text='outside the repository')
        reject(arm+': checkout itself', repo, arm, repo, text='outside the repository')
        existing = base / ('existing_'+arm); existing.mkdir()
        sentinel = existing/'sentinel'; sentinel.write_bytes(b'unchanged')
        reject(arm+': existing output', repo, arm, existing, text='outside the repository')
        require(sentinel.read_bytes() == b'unchanged', 'Existing output touched')
        local_source.write_bytes(original_source.read_bytes()+b'\n# tampered\n')
        reject(arm+': changed oracle', repo, arm, base/('source_fail_'+arm), text='Authenticated oracle or result bytes changed')
        local_source.write_bytes(original_source.read_bytes())
        local_result.write_bytes(original_result.read_bytes()+b' ')
        reject(arm+': changed result', repo, arm, base/('result_fail_'+arm), text='Authenticated oracle or result bytes changed')
        local_result.write_bytes(original_result.read_bytes())
        expected = json.loads(local_result.read_bytes())
        # Exercise the adapter's relocation and equality machinery with a
        # stubbed loader. The pinned oracle bytes are never modified or run.
        for mode in ('exact', 'runtime_only', 'scientific_drift', 'source_drift_during', 'result_drift_during'):
            out = base / (arm+'_'+mode)
            fake = SimpleNamespace()
            def exec_stub(module, mode=mode):
                def fake_main():
                    require(module.REPO == repo.resolve(), 'REPO relocation wrong')
                    require(module.HERE == out.resolve(), 'HERE relocation wrong')
                    value = copy.deepcopy(expected)
                    if mode == 'runtime_only': value['runtime'] = 'audit-stub-different-version'
                    if mode == 'scientific_drift': value['verdict'] = 'INTENTIONAL_MISMATCH'
                    (module.HERE/'independent_results.json').write_text(json.dumps(value),encoding='utf-8')
                    if mode == 'source_drift_during': local_source.write_bytes(original_source.read_bytes()+b'\n')
                    if mode == 'result_drift_during': local_result.write_bytes(original_result.read_bytes()+b' ')
                module.main = fake_main
            fake_spec = SimpleNamespace(loader=SimpleNamespace(exec_module=exec_stub))
            argv = [str(SOURCE),'--repo',str(repo),'--arm',arm,'--out',str(out)]
            error = None
            with patch.object(sys,'argv',argv), patch.object(adapter.importlib.util,'spec_from_file_location',return_value=fake_spec), patch.object(adapter.importlib.util,'module_from_spec',return_value=fake), contextlib.redirect_stdout(io.StringIO()):
                try: adapter.main()
                except RuntimeError as exc: error=str(exc)
            should_pass = mode in ('exact','runtime_only')
            require((error is None)==should_pass, arm+': '+mode+' unexpected result '+str(error))
            require((out/'replay_receipt.json').exists()==should_pass, 'False or missing pass receipt')
            if mode == 'scientific_drift': require('differs beyond' in error, 'Wrong scientific-drift failure')
            if mode.endswith('_during'): require('changed during replay' in error, 'Wrong post-run integrity failure')
            if should_pass:
                receipt=json.loads((out/'replay_receipt.json').read_text())
                require(receipt['adapter_sha256']==PIN and receipt['oracle_sha256']==pin['source_sha256'], 'Receipt provenance mismatch')
            local_source.write_bytes(original_source.read_bytes())
            local_result.write_bytes(original_result.read_bytes())
            record(arm+': stub '+mode)
        # Inspect independently authored imports and fixed saved input pins.
        tree=ast.parse(original_source.read_text(encoding='utf-8'))
        imports=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import): imports += [n.name.split('.')[0] for n in node.names]
            elif isinstance(node,ast.ImportFrom): imports.append((node.module or '').split('.')[0])
        require(set(imports) <= {'collections','fractions','hashlib','json','math','pathlib','platform'}, 'Non-stdlib/author oracle import')
        record(arm+': oracle imports standard library only')
        archived_receipt=REPO/'results/cycle5/portable_review_checks'/arm/'replay_receipt.json'
        r=json.loads(archived_receipt.read_bytes())
        require(r['adapter_sha256']==PIN and r['oracle_sha256']==pin['source_sha256'] and r['reference_result_sha256']==pin['result_sha256'] and r['result_sha256']==pin['result_sha256'], 'Root replay receipt pin mismatch')
        source_hashes[str(archived_receipt.relative_to(WORK))]=sha(archived_receipt)
        record(arm+': archived real replay receipt exact bytes')
    reject('missing root',base/'missing','development',base/'missingroot_out')
    root_file=base/'not_directory';root_file.write_text('not a root')
    reject('file root',root_file,'development',base/'fileroot_out')
    malformed=base/'empty_root';malformed.mkdir()
    reject('empty root',malformed,'development',base/'emptyroot_out')
    for flag in ('-O','-OO'):
        reject(flag+' before path access',base/'missing','development',base/('optimized_'+flag[1:]),flags=(flag,),text='Optimized Python is not permitted')
    reject('unrecognized arm',malformed,'not_an_arm',base/'badarm',text='invalid choice')
    require(str(base.resolve()).removeprefix('\\\\?\\').startswith(str(HERE.resolve())+os.sep), 'Temporary cleanup escaped owned directory')
    shutil.rmtree(base)
for path,pin in source_hashes.items(): require(sha(WORK/path)==pin,'Original changed: '+path)
record('all original reviewed hashes unchanged')
receipt={'verdict':'PASS_BOUNDED_ADAPTER_REVIEW','checks':checks,'count':len(checks),'source_hashes':source_hashes,'test_script_sha256':sha(Path(__file__)),'scientific_oracle_mains_run_in_this_review':False,'new_outcomes_accessed':False,'scope':'Guard CLI failures, loader-stub relocation and equality checks, AST imports and archived real-replay receipts. Scientific suites not rerun.'}
output=HERE/'checks.json'
if output.exists(): raise RuntimeError('Do not overwrite checks.json')
output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'verdict':receipt['verdict'],'count':len(checks),'receipt_sha256':sha(output)},indent=2))
