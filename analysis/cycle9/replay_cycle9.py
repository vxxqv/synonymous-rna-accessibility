"""Materialize only sealed numerical inputs and run unchanged independent checks."""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python is refused before input/output access')
import argparse, csv, gzip, hashlib, json, os, platform, subprocess
from pathlib import Path, PurePosixPath

def need(value, message):
    if not value:
        raise ValueError(message)

def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()

def long_path(path):
    s = str(Path(path).resolve())
    if os.name == 'nt' and not s.startswith('\\\\?\\'):
        s = '\\\\?\\UNC\\' + s[2:] if s.startswith('\\\\') else '\\\\?\\' + s
    return Path(s)

def safe(root, name):
    need(isinstance(name, str) and name and '\\' not in name and ':' not in name, 'Invalid relative path')
    p = PurePosixPath(name)
    need(not p.is_absolute() and all(x not in ('', '.', '..') for x in name.split('/')), 'Unsafe relative path')
    q = root.joinpath(*p.parts)
    for ancestor in [q, *q.parents]:
        need(not ancestor.is_symlink(), 'Symlink archive/input path refused')
        if ancestor == root:
            break
    return q

def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')

def normalize_pin_paths(pins):
    result = {}
    for path, h in pins.items():
        s = path.replace('\\', '/')
        if not s.startswith('work/'):
            need('/work/' in s, 'Unrecognized review input path')
            s = 'work/' + s.rsplit('/work/', 1)[1]
        need(s.startswith(('work/research_cycle8/', 'work/research_cycle9/')), 'Unexpected logical input namespace')
        need(s not in result or result[s] == h, 'Conflicting canonical input hash')
        result[s] = h
    return result

def compare(received, expected, name):
    # These unchanged verifiers include absolute path labels in their receipts.
    # Only those labels are canonicalized; every hash and scientific field stays.
    a, b = dict(received), dict(expected)
    for key in ('input_hashes_before_after',):
        if key in a or key in b:
            need(key in a and key in b, name + ': missing hash map')
            need(normalize_pin_paths(a.pop(key)) == normalize_pin_paths(b.pop(key)), name + ': input hash mismatch')
    need(a == b, name + ': numerical/scientific receipt differs')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, required=True, help='Archive-like root containing analysis/cycle9')
    ap.add_argument('--out', type=Path, required=True, help='New directory outside the supplied archive root')
    ap.add_argument('--manifest-sha256', required=True, help='Trusted FILES.json SHA256 from the reviewed handoff')
    ap.add_argument('--verify-only', action='store_true', help='Check all archive bytes without creating output')
    a = ap.parse_args()
    need(platform.python_implementation() == 'CPython' and sys.version_info[:2] == (3,12), 'CPython3.12 required')
    repo, out = long_path(a.repo), long_path(a.out)
    need(not out.exists(), 'Fresh output directory required')
    need(repo != out and repo not in out.parents and out not in repo.parents, 'Output must be outside archive and its ancestors')
    mp = safe(repo, 'analysis/cycle9/FILES.json')
    mb = mp.read_bytes()
    need(sha_bytes(mb) == a.manifest_sha256.lower(), 'Manifest hash mismatch')
    man = json.loads(mb)
    need(man['schema'] == 'Cycle9_offline_numerical_replay_v1', 'Manifest schema')
    seen, destinations, restored, logical = set(), set(), [], {}
    for row in man['files']:
        name = row['archive_path']
        need(name.startswith('analysis/cycle9/'), 'Unexpected archive namespace')
        need(name.casefold() not in seen, 'Duplicate archive path')
        seen.add(name.casefold())
        p = safe(repo, name)
        b = p.read_bytes()
        need(len(b) == row['stored_bytes'] and sha_bytes(b) == row['stored_sha256'], 'Stored hash/size mismatch: ' + name)
        need(row['compression'] in ('none', 'gzip'), 'Unsupported compression')
        z = gzip.decompress(b) if row['compression'] == 'gzip' else b
        need(len(z) == row['logical_bytes'] and sha_bytes(z) == row['logical_sha256'], 'Logical hash/size mismatch: ' + name)
        logical[name] = z
        if row['restore_path'] is not None:
            q = safe(out, row['restore_path'])
            need(row['restore_path'].casefold() not in destinations, 'Duplicate restoration path')
            destinations.add(row['restore_path'].casefold())
            restored.append((q, z, row['restore_path'], row['logical_sha256']))
    actual = {p.relative_to(repo).as_posix().casefold() for p in safe(repo, 'analysis/cycle9').rglob('*') if p.is_file()}
    need(actual == seen | {'analysis/cycle9/files.json'}, 'Uninventoried or missing Cycle9 archive file')
    for e in man['checks']:
        need(e['script'].casefold() in destinations, 'Checker not restored')
        safe(out, e['output'])
        need(e['reference_archive_path'] in logical, 'Reference receipt absent')
    if a.verify_only:
        print(json.dumps({'status':'PASS_ARCHIVE_MANIFEST','files':len(seen),'manifest_sha256':a.manifest_sha256}))
        return
    out.mkdir(parents=True)
    for dest, b, name, h in restored:
        dest.parent.mkdir(parents=True, exist_ok=True)
        need(not dest.exists(), 'Restoration collision')
        dest.write_bytes(b)
    logdir = out/'replay_logs'
    logdir.mkdir()
    results = []
    try:
        for e in man['checks']:
            log = logdir/(e['name']+'.log')
            script = safe(out, e['script'])
            with log.open('w', encoding='utf-8') as stream:
                run = subprocess.run([sys.executable, '-B', str(script)], cwd=out, stdout=stream, stderr=subprocess.STDOUT)
            need(run.returncode == 0, e['name'] + ' failed; retained log: ' + str(log))
            path = safe(out, e['output'])
            got = json.loads(path.read_text(encoding='utf-8'))
            expected = json.loads(logical[e['reference_archive_path']])
            compare(got, expected, e['name'])
            for key, value in e['required'].items():
                need(got[key] == value, e['name'] + ': required coverage mismatch')
            results.append({'name':e['name'],'status':'PASS_EXACT_SCIENTIFIC_RECEIPT','coverage':e['required'],
                            'receipt_sha256':sha_bytes(path.read_bytes()),'unchanged_checker_sha256':sha_bytes(script.read_bytes())})
        # An additional exported table must agree including every profile row.
        extra = man['profile_csv']
        with safe(out, extra['output']).open(newline='', encoding='utf-8') as f:
            got = list(csv.DictReader(f))
        import io
        expected = list(csv.DictReader(io.StringIO(logical[extra['reference_archive_path']].decode('utf-8'))))
        need(got == expected and len(got) == 60, 'Complete profile CSV mismatch')
        for p, b, name, h in restored:
            need(sha_bytes(p.read_bytes()) == h, 'Input changed during replay: ' + name)
        # Recheck original archive bytes, including code, after all children.
        for row in man['files']:
            need(sha_bytes(safe(repo,row['archive_path']).read_bytes()) == row['stored_sha256'], 'Archive mutated during replay')
        result = {'status':'PASS_ALL_THREE_CYCLE9_NUMERICAL_REPLAYS','manifest_sha256':a.manifest_sha256,
                  'archive_files':len(seen),'restored_files':len(restored),'checks':results,'profile_csv_rows':60,
                  'runtime':{'implementation':platform.python_implementation(),'python':platform.python_version(),'platform':platform.platform()},
                  'network_used':False,'git_used':False,'author_models_or_RF_rerun':False,'frozen_files_modified':False,
                  'receipt_comparison':'All parsed fields equal; absolute input path labels canonicalized to work-relative suffixes, all hash values retained.',
                  'scope':'Replay of existing exposed numerical verifications only, not new outcomes or historical access-order authentication.'}
        dump(out/'replay_receipt.json',result)
        print(json.dumps(result,indent=2))
    except BaseException as exc:
        dump(out/'replay_failure.json',{'status':'FAILED_PRESERVED','completed_checks':results,'error':repr(exc)})
        raise

if __name__ == '__main__':
    main()
