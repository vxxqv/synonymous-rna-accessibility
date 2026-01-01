"""Frozen metadata/configuration preflight only. No sequence or folding calls."""
import sys
if not __debug__:
    raise RuntimeError('Disabled assertions are refused before input access')
sys.dont_write_bytecode = True
import argparse
import hashlib
import importlib.util
import json
import platform
import traceback
import warnings
from pathlib import Path

PREFIX = 'analysis/programme/retrospective/persist_v2/'
PINNED = {
    PREFIX + 'protocol.json': 'a3dc686a23acaa12af2f2f44154b8baac819cb6df3a0a6d81328dd3ebdc5fd24',
    PREFIX + 'persist_prediction.py': '962f953ed17574db96f5ba0cde4827486830244941ea745bb34fb697a5514aaa',
    PREFIX + 'Turner2004.par': 'f4a71a542533ce5420515004a0b51f6949e93d4087204185a35eb38fff30d466',
}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def parameter_text_bytes(raw):
    return raw.replace(b'\r\n', b'\n')

def preflight(repo, output):
    repo, output = Path(repo).resolve(), Path(output).resolve()
    if output.exists() or output.is_relative_to(repo):
        raise ValueError('Output must be fresh and outside the input repository')
    hashes = {name: digest(repo / name) for name in PINNED}
    if hashes != PINNED:
        raise ValueError('Frozen configuration/source hash mismatch')
    output.mkdir(parents=True)
    report = {'status': 'started', 'input_sha256_before': hashes,
              'python': platform.python_version(), 'platform': platform.platform(),
              'folds': 0, 'sequence_rows_read': 0, 'outcome_reads': 0,
              'stage': 'version_import'}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            if platform.python_implementation() != 'CPython' or sys.version_info[:2] != (3, 12):
                raise RuntimeError('CPython3.12 is required')
            import numpy as np
            import RNA
            report.update({'numpy': np.__version__, 'ViennaRNA': RNA.__version__,
                           'pf_float_precision_code': RNA.pf_float_precision()})
            if np.__version__ != '2.5.4' or RNA.__version__ != '2.7.2' or RNA.pf_float_precision() != 0:
                raise RuntimeError('Frozen version or PF precision differs')
            report['stage'] = 'author_import'
            spec = importlib.util.spec_from_file_location('frozen_persist_config', repo / (PREFIX + 'persist_prediction.py'))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            protocol = json.loads((repo / (PREFIX + 'protocol.json')).read_text(encoding='utf-8'))
            report['stage'] = 'configured_md'
            loaded, md = module.configured_md(protocol['vienna'])
            report['stage'] = 'snapshot_md'
            snapshot = module.snapshot_md(md, loaded)
            report['effective_md'] = snapshot
            if snapshot != protocol['vienna']['effective_md']:
                raise RuntimeError('Frozen effective MD mismatch')
            report['stage'] = 'parameter_export'
            exported = output / 'Turner2004.exported.par'
            if not bool(loaded.params_save(str(exported))):
                raise RuntimeError('Parameter export failed')
            report['exported_parameter_sha256'] = digest(exported)
            report['parameter_comparison'] = 'Exact after CRLF to LF only; raw input and export hashes retained'
            report['exported_parameter_LF_sha256'] = hashlib.sha256(parameter_text_bytes(exported.read_bytes())).hexdigest()
            if parameter_text_bytes(exported.read_bytes()) != parameter_text_bytes((repo / (PREFIX + 'Turner2004.par')).read_bytes()):
                raise RuntimeError('Frozen parameter byte mismatch')
            report['stage'] = 'complete_without_folding'
            report['status'] = 'configuration_only_passed'
    except BaseException as exc:
        report['status'] = 'failed'
        report['exception'] = repr(exc)
        report['traceback'] = traceback.format_exc()
        raise
    finally:
        report['input_sha256_after'] = {name: digest(repo / name) for name in PINNED}
        report['inputs_immutable'] = report['input_sha256_after'] == hashes
        (output / 'preflight.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    preflight(args.repo_root, args.output)

if __name__ == '__main__':
    main()
