from pathlib import Path
import json,urllib.request,concurrent.futures,hashlib,zipfile
P=Path(__file__).resolve().parent;c=json.loads((P/'dependency_closure_proposal.json').read_text());sizes=json.loads((P/'dependency_size_check.json').read_text())
if not sizes['all_sizes_known'] or sizes['total_bytes']>100_000_000:raise RuntimeError('acquisition cost gate')
d=P/'runtime'/'downloads';d.mkdir(parents=True,exist_ok=True)
def get(item):
 n,p=item;url=f'https://cran.r-project.org/bin/windows/contrib/4.6/{n}_{p["Version"]}.zip'; dest=d/url.rsplit('/',1)[1]
 if dest.exists():b=dest.read_bytes()
 else:
  with urllib.request.urlopen(url,timeout=45) as r:b=r.read(40_000_001)
 if len(b)>40_000_000:raise RuntimeError('package size cap')
 if not dest.exists():dest.write_bytes(b)
 with zipfile.ZipFile(dest) as z:
  if sum(x.file_size for x in z.infolist())>200_000_000:raise RuntimeError('expanded package size cap')
  desc=z.read(n+'/DESCRIPTION').decode('utf8').replace('\r\n','\n')
  if '\nVersion: '+p['Version']+'\n' not in desc:raise RuntimeError('package version mismatch '+n)
  if '\nPackage: '+n+'\n' not in '\n'+desc:raise RuntimeError('package name mismatch '+n)
  md5checks=0
  for line in z.read(n+'/MD5').decode().splitlines():
   digest,filename=line.split(None,1);filename=filename.strip().lstrip('*').strip('"')
   if hashlib.md5(z.read(n+'/'+filename)).hexdigest()!=digest:raise RuntimeError('package internal checksum '+n+'/'+filename)
   md5checks+=1
 return {'package':n,'version':p['Version'],'file':str(dest.relative_to(P)).replace('\\','/'),'url':url,'bytes':len(b),'internal_md5_checks':md5checks,'sha256':hashlib.sha256(b).hexdigest(),'external_digest_available':False}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(get,c['binary_packages'].items()))
(P/'runtime'/'download_manifest.json').write_text(json.dumps({'packages':r,'total_bytes':sum(x['bytes'] for x in r)},indent=2)+'\n')
print('Downloaded and verified',len(r),sum(x['bytes'] for x in r))
