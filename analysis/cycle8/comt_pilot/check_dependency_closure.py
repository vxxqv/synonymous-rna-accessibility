from inspect_dependency_sources import get,P
import gzip,json,re,csv
b=get('https://cran.r-project.org/bin/windows/contrib/4.6/PACKAGES.gz','CRAN_WINDOWS_PACKAGES.gz')
def parse(b):
 out={}
 for raw in gzip.decompress(b).decode().split('\n\n'):
  d={}; last=None
  for l in raw.splitlines():
   if l[:1].isspace() and last:d[last]+=' '+l.strip()
   elif ': ' in l:last,x=l.split(': ',1);d[last]=x
  if 'Package' in d:out[d['Package']]=d
 return out
pk=parse(b);inst={x['Package']:x for x in csv.DictReader((P/'installed_R_packages.tsv').open(),delimiter='\t')}
roots=['data.table','caret','randomForest','e1071']; todo=roots[:]; selected={}; installed={}
while todo:
 name=todo.pop()
 if name in selected or name in installed or name=='R':continue
 if name in inst:
  installed[name]=inst[name];continue
 d=pk[name];selected[name]=d
 for key in ['Depends','Imports','LinkingTo']:
  for x in d.get(key,'').split(','):
   dep=re.split(r'\s|\(',x.strip())[0]
   if dep:todo.append(dep)
result={'roots':roots,'not_a_final_runtime_lock':True,'binary_packages':dict(sorted(selected.items())),'system_R_packages':dict(sorted(installed.items())),'count_packages_to_acquire':len(selected)}
(P/'dependency_closure_proposal.json').write_text(json.dumps(result,indent=2)+'\n')
print('Count',len(selected));print([(x,d['Version']) for x,d in sorted(selected.items())]);print('Installed',[(x,d['Version']) for x,d in sorted(installed.items())])
