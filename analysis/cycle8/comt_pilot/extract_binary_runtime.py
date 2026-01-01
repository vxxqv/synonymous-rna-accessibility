from pathlib import Path,PurePosixPath
import json,hashlib,zipfile,os
P=Path(__file__).resolve().parent
m=json.loads((P/'runtime/download_manifest.json').read_text());lib=P/'runtime/library'
if lib.exists():raise RuntimeError('fresh library required')
lib.mkdir()
n=0
for r in m['packages']:
 p=P/r['file']
 if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError('changed package')
 with zipfile.ZipFile(p) as z:
  for info in z.infolist():
   rel=PurePosixPath(info.filename)
   if rel.is_absolute() or '..' in rel.parts or rel.parts[0]!=r['package']:raise RuntimeError('unsafe archive member')
   target=lib.joinpath(*rel.parts)
   ext=Path('\\\\?\\'+str(target.resolve())) if os.name=='nt' else target
   if info.is_dir():ext.mkdir(parents=True,exist_ok=True)
   else:
    ext.parent.mkdir(parents=True,exist_ok=True);ext.write_bytes(z.read(info));n+=1
(P/'runtime/extraction_receipt.json').write_text(json.dumps({'method':'Verified official Windows R binary ZIP extraction into task-local library','package_archives':len(m['packages']),'files':n,'download_manifest_sha256':hashlib.sha256((P/'runtime/download_manifest.json').read_bytes()).hexdigest(),'historical_installer_failure_preserved':'library_failed_longpath','no_global_library_modifications':True},indent=2)+'\n')
print('Extracted',n)
