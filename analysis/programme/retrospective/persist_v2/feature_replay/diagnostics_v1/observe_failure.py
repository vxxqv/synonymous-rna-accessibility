"""Read only the feature replay's generated diagnostics. Never import science code."""
import sys
if not __debug__:
    raise RuntimeError('Disabled assertions are refused before input access')
import argparse
import hashlib
import json
from pathlib import Path

ROOT_FILES = ('STARTED.json', 'environment.json', 'FAILURE.json')
WORKER_FILES = ('stdout.txt', 'stderr.txt', 'row.json', 'Turner2004.exported.par')

def observe(root):
    root = Path(root).resolve()
    if not root.exists():
        return {'status': 'diagnostic_directory_absent', 'root': str(root), 'files': []}
    if not root.is_dir():
        raise ValueError('Diagnostic root must be a directory')
    candidates = [(root / name, name) for name in ROOT_FILES]
    for child in sorted(root.iterdir()):
        if child.is_symlink():
            raise ValueError('Symlink in diagnostic root is refused')
        if child.is_dir():
            candidates.extend((child / name, child.name + '/' + name) for name in WORKER_FILES)
    result = []
    for path, relative in candidates:
        if path.is_symlink():
            raise ValueError('Symlink diagnostic file is refused')
        if not path.exists():
            continue
        if not path.is_file() or not path.resolve().is_relative_to(root):
            raise ValueError('Diagnostic path is not a contained regular file')
        raw = path.read_bytes()
        item = {'file': relative, 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
        # Logs are JSON escaped so embedded GitHub workflow commands stay inert.
        # Full original stream bytes remain in the uploaded replay directory.
        if path.name != 'Turner2004.exported.par':
            item['text_utf8_with_replacement'] = raw.decode('utf-8', errors='replace')
        result.append(item)
    return {'status': 'observation_only', 'root': str(root), 'files': result,
            'scientific_modules_imported': False, 'worker_reruns': 0,
            'full_stream_bytes_require_artifact_upload': True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    # Each line starts with JSON syntax, never with an executable workflow command.
    print(json.dumps(observe(args.root), ensure_ascii=True, allow_nan=False))

if __name__ == '__main__':
    main()
