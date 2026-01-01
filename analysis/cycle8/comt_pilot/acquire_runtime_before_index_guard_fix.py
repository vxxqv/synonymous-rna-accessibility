from pathlib import Path
import json,urllib.request,concurrent.futures,hashlib,zipfile
P=Path(__file__).resolve().parent;c=json.loads((P/'dependency_closure_proposal.json').read_text());sizes=json.loads((P/'dependency_size_check.json').read_text())
if not sizes['all_sizes_known'] or sizes['total_bytes']>100_000_000:raise RuntimeError('acquisition cost gate')
d=P/'runtime'/'downloads';d.mkdir(parents=True,exist_ok=True)
def get(item):
 n,p=item;url=f'https://cran.r-project.org/bin/windows/contrib/4.6/{n}_{p["Version"]}.zip'; dest=d/url.rsplit('/',1)[1]
 if dest.exists():raise RuntimeError('fresh download guard '+n)
 with urllib.request.urlopen(url,timeout=45) as r:b=r.read(40_000_001)
 if len(b)>40_000_000:raise RuntimeError('package size cap')
 if hashlib.md5(b).hexdigest()!=p['MD5sum']:raise RuntimeError('package index digest '+n)
 dest.write_bytes(b)
 with zipfile.ZipFile(dest) as z:
  if sum(x.file_size for x in z.infolist())>200_000_000:raise RuntimeError('expanded package size cap')
 return {'package':n,'version':p['Version'],'file':str(dest.relative_to(P)).replace('\\','/'),'url':url,'bytes':len(b),'md5':p['MD5sum'],'sha256':hashlib.sha256(b).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(get,c['binary_packages'].items()))
(P/'runtime'/'download_manifest.json').write_text(json.dumps({'packages':r,'total_bytes':sum(x['bytes'] for x in r)},indent=2)+'\n')
print('Downloaded and verified',len(r),sum(x['bytes'] for x in r))
