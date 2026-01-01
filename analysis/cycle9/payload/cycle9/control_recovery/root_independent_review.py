"""Fresh independent reconstruction from the 18 approved XLSX numeric cells."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib, json, posixpath, zipfile, xml.etree.ElementTree as E
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RNS='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
def require(x,s):
 if not x:raise ValueError(s)
def main():
 if not __debug__:raise RuntimeError('Optimized Python refused')
 out=HERE/'root_independent_review';require(not out.exists(),'Refuse overwrite');out.mkdir()
 source=ROOT/'work/research_cycle8/cprime_gate/supp_Supplementary Data 6.xlsx'
 h=hashlib.sha256(source.read_bytes()).hexdigest()
 require(h=='a3d18aff825b756711de622533cb730c1202e571f86a3b9b2e00edea6b073dd3','workbook hash')
 d=json.loads((HERE/'results_v1/result.json').read_text())
 allowed={f'{c}{r}' for c in 'FGHIJK' for r in range(15,18)}
 expected={'C15':'R1207R_ATT(G/A)','C16':'R1207R_(A/G)TTG','C17':'R1207R_(A/G)TT(G/A)'}
 with zipfile.ZipFile(source) as z:
  wb=E.fromstring(z.read('xl/workbook.xml'));sheet=next(s for s in wb.findall('m:sheets/m:sheet',NS) if s.attrib['name']=='validation K562')
  rels=E.fromstring(z.read('xl/_rels/workbook.xml.rels'));rid=sheet.attrib['{'+RNS+'}id']
  target=next(r.attrib['Target'] for r in rels if r.attrib['Id']==rid)
  member=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
  strings=[''.join(t.itertext()) for s in E.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS) for t in [s]]
  # Identity pass never converts any numeric cells.
  cells=E.fromstring(z.read(member));seen={}
  for c in cells.findall('.//m:c',NS):
   a=c.attrib.get('r')
   if a in expected:
    require(c.attrib.get('t')=='s','label representation');seen[a]=strings[int(c.find('m:v',NS).text)]
  require(seen==expected,'label contract')
  values={}
  for c in cells.findall('.//m:c',NS):
   a=c.attrib.get('r')
   if a not in allowed:continue
   require(c.attrib.get('t','n')=='n','numeric type')
   v=c.find('m:v',NS);require(v is not None and v.text is not None,'stored cache')
   values[a]=F(v.text);require(values[a]>0,'nonpositive admitted source')
 require(set(values)==allowed,'cell coverage')
 require({r['cell']:F(r['exact']) for r in d['source_cells']}==values,'lexical extraction agreement')
 state={}
 for name,row in [('M',15),('E',16),('C',17)]:
  state[name]={(i+1,j+1):values[f'{rc}{row}']/values[f'{dc}{row}'] for i,rc in enumerate('IJK') for j,dc in enumerate('FGH')}
  records=next(r for r in d['state_results'] if r['state']==name)
  require({(r['RNA_index'],r['DNA_index']):F(r['value']['exact']) for r in records['ratios']}==state[name],'all state ratios')
  require(F(records['envelope']['lower']['exact'])==min(state[name].values()),'state lower')
  require(F(records['envelope']['upper']['exact'])==max(state[name].values()),'state upper')
 # A fresh Cartesian reconstruction of the joint space preserves the common C value.
 pair={'M':[],'E':[]};joint=[]
 for c,m,e in product(state['C'].values(),state['M'].values(),state['E'].values()):joint.append((c/m,c/e))
 require(len(joint)==729,'joint denominator')
 for name in ('M','E'):
  pair[name]=[c/x for c,x in product(state['C'].values(),state[name].values())]
  require(len(pair[name])==81,'contrast denominator')
  require(max(pair[name])<F(1,2),'uniform below half')
 require(all(a<F(1,2) and b<F(1,2) for a,b in joint),'joint secondary classification')
 # Compare every reported quotient by its numeric multiset, preserving multiplicity.
 actual=sorted(F(r['value']['exact']) for r in d['contrast_quotients'])
 require(actual==sorted(pair['M']+pair['E']),'all162 contrasts')
 # Synthetic algebraic guards: exact equality is not strict, shared C cannot split.
 require(not(F(1)<1) and F(1,2)<=F(1,2),'strict/weak boundary')
 toy=[(c/F(1),c/F(2)) for c in [F(1),F(3)]]
 require((F(1),F(3,2)) not in toy,'shared-C coupling')
 result={'status':'PASS','workbook_sha256':h,'result_sha256':hashlib.sha256((HERE/'results_v1/result.json').read_bytes()).hexdigest(),
 'source_cells':18,'state_ratios':27,'contrast_quotients':162,'joint_completions':729,
 'range':{s:{'lower':str(min(v)),'upper':str(max(v)),'lower_float':float(min(v)),'upper_float':float(max(v))} for s,v in pair.items()},
 'scope':'Fresh independent arithmetic from allowlisted source cells; descriptive source-state contrast, no biological independence or isolated-edit inference.',
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
