"""Record inspection/proposal hashes. Does not execute scientific code."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2] / 'outputs' / 'synonymous-rna-accessibility'
# HERE is work/persist_cycle3/linux_feature_diagnostic; workspace is parents[2].
FILES = (
    '.github/workflows/verify.yml',
    'requirements-tested.txt',
    'analysis/programme/retrospective/persist_v2/feature_replay/replay.py',
    'analysis/programme/retrospective/persist_v2/feature_replay/manifest.json',
    'analysis/programme/retrospective/persist_v2/feature_replay/README.md',
    'analysis/programme/retrospective/persist_v2/persist_prediction.py',
    'analysis/programme/retrospective/persist_v2/protocol.json',
    'analysis/programme/retrospective/persist_v2/Turner2004.par',
    'results/programme/persist_v2/feature_replay/verification_receipt.json',
)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

manifest = {
    'scope': 'Read-only diagnostic investigation and synthetic tests; no actual folding',
    'repository_review_revision': '602ac0559efcc9c5e10410d85f6a164761dc2b3f',
    'github_run_id': 40178313369, 'github_job_id': 121123714821,
    'retained_artifact_count': 0,
    'actual_child_stdout_stderr_available': False,
    'actual_child_exit_code_available': False,
    'cause': 'unresolved; no diagnosis inferred from generic parent exception',
    'real_folds_executed': 0, 'outcome_reads': 0,
    'shared_files_modified': [],
    'inspected_repository_sha256': {name: sha(REPO / name) for name in FILES},
    'delivered_sha256': {str(path.relative_to(HERE)).replace('\\', '/'): sha(path)
                         for path in sorted(HERE.iterdir())
                         if path.is_file() and path.name != 'manifest.json'},
}
with (HERE / 'manifest.json').open('x', encoding='utf-8') as stream:
    json.dump(manifest, stream, indent=2); stream.write('\n')
print(json.dumps({'manifest_sha256': sha(HERE / 'manifest.json'),
                  'inspected_files': len(FILES), 'delivered_files': len(manifest['delivered_sha256'])}))
