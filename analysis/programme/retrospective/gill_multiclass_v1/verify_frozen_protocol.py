"""Replay the unchanged, outcome-blind root check in a temporary directory."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
args = parser.parse_args()
if args.out.exists():
    raise FileExistsError('Use a fresh result path')
manifest = json.loads((HERE/'freeze_manifest.json').read_text(encoding='utf-8'))
for name, expected in manifest['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, name
with tempfile.TemporaryDirectory(prefix='gill_protocol_check_') as folder:
    root = Path(folder)
    (root/'inventory.json').write_bytes(gzip.decompress((HERE/'inventory.json.gz').read_bytes()))
    shutil.copyfile(HERE/'freeze_proposal.json',root/'freeze_proposal.json')
    shutil.copyfile(HERE/'root_protocol_review_original.py',root/'root_protocol_review.py')
    result = subprocess.run([sys.executable,str(root/'root_protocol_review.py')],check=True,capture_output=True,text=True)
    actual = json.loads((root/'root_protocol_checks.json').read_text(encoding='utf-8'))
    saved = json.loads((HERE/'root_protocol_checks.json').read_text(encoding='utf-8'))
    assert actual == saved, 'Independent sequence-only check changed'
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(actual,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'exact_match','outcome_values_opened':False,'classes':actual['classes']}))
