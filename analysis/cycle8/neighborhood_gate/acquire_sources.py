from pathlib import Path
import requests,json,hashlib,xml.etree.ElementTree as E

if not __debug__:raise RuntimeError('Optimized Python prohibited')
p=Path('work/research_cycle8/neighborhood_gate');p.mkdir(parents=True,exist_ok=True)
receipts=[]
for name,url in [
 ('tp53.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5626615/fullTextXML'),
 ('schmitz_crossref.json','https://api.crossref.org/works/10.1186/s12864-021-07462-z'),
 ('tp53_crossref.json','https://api.crossref.org/works/10.1158/1541-7786.MCR-17-0245')]:
 dst=p/name
 if dst.exists():raise FileExistsError(dst)
 try:
  r=requests.get(url,timeout=25);dst.write_bytes(r.content);receipts.append({'file':name,'url':url,'final_url':r.url,'status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()})
 except requests.RequestException as e:receipts.append({'file':name,'url':url,'error':str(e)})
src=Path('work/joint_function_cycle5/splicing_sources/schmitz.xml');(p/'schmitz.xml').write_bytes(src.read_bytes());receipts.append({'file':'schmitz.xml','reused_from':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
(p/'acquisition_receipt.json').write_text(json.dumps(receipts,indent=2));print(json.dumps(receipts,indent=2))
for name in ['schmitz','tp53']:
 try:r=E.fromstring((p/(name+'.xml')).read_bytes())
 except E.ParseError:continue
 print('TITLE',name,' '.join(r.find('.//article-title').itertext()))
 for sec in r.findall('.//sec'):
  ti=sec.find('title');title=' '.join(ti.itertext()) if ti is not None else ''
  print('SECTION',name,title)
 selected=[]
 for sec in r.findall('.//sec'):
  ti=sec.find('title');title=' '.join(ti.itertext()) if ti is not None else ''
  if 'method' in title.lower() or 'availability' in title.lower() or 'design of a synthetic' in title.lower():selected.append(' '.join(sec.itertext()))
 (p/(name+'_methods.txt')).write_text('\n\n'.join(selected),encoding='utf-8')
 for supp in r.findall('.//supplementary-material'):print('SUPPLEMENT',name,E.tostring(supp,encoding='unicode')[:1100])
