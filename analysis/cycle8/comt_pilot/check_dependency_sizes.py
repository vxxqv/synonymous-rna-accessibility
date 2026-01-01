from pathlib import Path
import json,urllib.request,concurrent.futures
P=Path(__file__).resolve().parent
c=json.loads((P/'dependency_closure_proposal.json').read_text())
def head(item):
 n,d=item; url=f'https://cran.r-project.org/bin/windows/contrib/4.6/{n}_{d["Version"]}.zip'
 try:
  with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=20) as r:return {'name':n,'url':url,'bytes':int(r.headers.get('Content-Length',-1)),'status':r.status}
 except Exception as e:return {'name':n,'url':url,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(head,c['binary_packages'].items()))
out={'files':r,'total_bytes':sum(x.get('bytes',0) for x in r),'all_sizes_known':all(x.get('bytes',-1)>=0 for x in r)}
(P/'dependency_size_check.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k!='files'});print([x for x in r if 'error' in x])
