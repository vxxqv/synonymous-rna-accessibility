"""Selective, hash-authorized cPRIME finite-set analysis. No automatic acquisition."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
from decimal import Decimal,InvalidOperation
from statistics import median
import argparse,csv,hashlib,io,json,posixpath,zipfile,xml.etree.ElementTree as ET

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    from interval_core import jsonable
    p.write_text(json.dumps(jsonable(v),indent=2)+'\n',encoding='utf-8')

def read_rows(path,eligible,contexts):
    """Decode only row IDs/context and the three declared paired outcome columns."""
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)
        book=ET.fromstring(z.read('xl/workbook.xml'))
        sheet=next(s for s in book.find('s:sheets',NS) if s.attrib['name']=='Supplementary_Table_3_Multi_exo')
        rid=sheet.attrib['{'+REL+'}id']
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        target=next(n.attrib['Target'] for n in rels if n.attrib['Id']==rid)
        xmlpath=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        def decode(cell):
            if cell is None:return None
            if cell.find('s:f',NS) is not None:raise ValueError('Outcome formulas require a separate contract')
            typ=cell.attrib.get('t')
            v=cell.find('s:v',NS)
            if typ=='inlineStr':return ''.join(cell.find('s:is',NS).itertext())
            if v is None:return None
            if typ=='s':return ''.join(strings[int(v.text)].itertext())
            if typ=='e':return 'EXCEL_ERROR'
            return v.text
        def value(cell):
            s=decode(cell)
            if s is None or str(s).strip().lower() in {'','na','nan','none','n/a','excel_error'}:return None
            try:d=Decimal(s)
            except InvalidOperation:raise ValueError('Unexpected nonnumeric outcome encoding')
            return F(d) if d.is_finite() else None
        rows=[];headers=None
        for ev,row in ET.iterparse(io.BytesIO(z.read(xmlpath)),events=('end',)):
            if row.tag!='{'+NS['s']+'}row':continue
            cells={''.join(c for c in x.attrib['r'] if c.isalpha()):x for x in row}
            if row.attrib.get('r')=='1':
                expected={'A':'combined_var_id','B':'var_id','E':'exon_number','F':'syn_var','G':'cell_line','H':'RPM_DNA','I':'RPM_RNA','J':'cprime_score'}
                headers={k:decode(cells.get(k)) for k in expected}
                if any(headers.get(k)!=v for k,v in expected.items()):raise ValueError('Source schema mismatch')
                row.clear();continue
            key=decode(cells.get('A'))
            if key not in eligible:row.clear();continue
            cell=decode(cells.get('G'));record={'source_row':int(row.attrib['r']),'state':key,'cell':cell,
             'intended':decode(cells.get('B')),'marker':decode(cells.get('F')),'exon':int(Decimal(decode(cells.get('E'))))}
            info=eligible[key]
            if cell not in contexts or record['intended']!=info['var_id'] or record['marker']!=info['syn_var'] or record['exon']!=int(info['exon_number']):
                raise ValueError('Source identity mismatch before outcome decoding')
            for field,col in [('DNA','H'),('RNA','I'),('source_score','J')]:record[field]=value(cells.get(col))
            ratio_available=all(record[f] is not None and record[f]>0 for f in ['DNA','RNA'])
            record['ratio']=record['RNA']/record['DNA'] if ratio_available else None
            compare_available=ratio_available and record['source_score'] is not None
            record['score_difference']=record['source_score']-record['ratio'] if compare_available else None
            record['source_score_consistent']=abs(record['score_difference'])<=F(1,10**8)*max(F(1),abs(record['ratio'])) if compare_available else None
            usable=ratio_available and record['source_score'] is not None and record['source_score']>0
            record['usable']=usable
            record['missing_reason']=None if usable else 'missing_nonfinite_or_nonpositive_required_value'
            rows.append(record);row.clear()
    return rows

def analyze(rows,eligibility):
    from interval_core import marker_bounds,common_order,widen
    states=eligibility['states'];contexts=eligibility['cell_contexts'];groups=eligibility['groups']
    by=defaultdict(list)
    for row in rows:
        key=row['state'];info=states[key]
        if row['cell'] not in contexts or row['intended']!=info['var_id'] or row['marker']!=info['syn_var'] or row['exon']!=int(info['exon_number']):
            raise ValueError('Source identity mismatch')
        by[key,row['cell']].append(row)
    if any(len(by[k,c])!=3 for k in states for c in contexts):raise ValueError('Expected three metadata rows per state/cell')
    if len(rows)!=len(states)*len(contexts)*3:raise ValueError('Extra source rows')
    bad=[r for r in rows if r['source_score_consistent'] is False]
    if bad:return {'status':'SOURCE_REPRODUCTION_FAILED','source_rows':rows,'inconsistent_rows':[r['source_row'] for r in bad]}
    observations=[];fields={c:{} for c in contexts}
    for k in sorted(states):
        for c in contexts:
            rs=by[k,c];qs=[r['ratio'] for r in rs if r['usable']]
            complete=len(qs)==3;pair=(min(qs),max(qs)) if complete else (F(0),None)
            fields[c][k]=pair
            observations.append({'state':k,'cell':c,'usable_rows':len(qs),'complete':complete,'lower':pair[0],'upper':pair[1],
                                 'median':median(qs) if complete else None,'source_rows':[r['source_row'] for r in rs]})
    obs={(r['state'],r['cell']):r for r in observations}
    sensitivity=[]
    for gamma in [F(1),F(5,4),F(3,2),F(2)]:
        fs={c:{v:widen(p,gamma) for v,p in d.items()} for c,d in fields.items()}
        comparisons=[]
        for e,ids in sorted(groups.items()):
            for c in contexts:
                b=marker_bounds([fs[c][v] for v in ids]);points=[obs[v,c]['median'] for v in ids]
                comparisons.append({'intended':e,'cell':c,'markers':ids,'marker_count':len(ids),
                  'complete':all(obs[v,c]['complete'] for v in ids),'minimum_fold':b['minimum'],'maximum_fold':b['maximum'],
                  'maximum_attained':b['maximum_attained'],'point_fold':max(points)/min(points) if all(x is not None for x in points) else None,
                  'common_response_feasible':b['minimum']==1,
                  'thresholds':{str(t):{'minimum_at_least':b['minimum']>=t,'maximum_at_most':b['maximum'] is not None and b['maximum']<=t} for t in [F(5,4),F(3,2),F(2),F(4)]}})
        order=[]
        for exon in eligibility['exons']:
            ids=sorted(v for v in states if int(states[v]['exon_number'])==exon)
            test=common_order({c:{v:fs[c][v] for v in ids} for c in contexts})
            order.append({'exon':exon,'state_count':len(ids),**test})
        sensitivity.append({'gamma':gamma,'group_cell_results':comparisons,'exon_shared_order_results':order,
         'summary':{'group_cells':len(comparisons),'complete_group_cells':sum(r['complete'] for r in comparisons),
          'strictly_positive_lower_dispersion':sum(r['minimum_fold']>1 for r in comparisons),
          'minimum_at_least_twofold':sum(r['minimum_fold']>=2 for r in comparisons),
          'maximum_at_most_twofold':sum(r['maximum_fold'] is not None and r['maximum_fold']<=2 for r in comparisons),
          'shared_order_infeasible_exons':[r['exon'] for r in order if r['status']=='infeasible']}})
    return {'status':'COMPLETE_DESCRIPTIVE_AUDIT','source_rows':rows,'state_cells':observations,'sensitivities':sensitivity,
     'scope':'Finite local-state proxy ranges, not calibrated confidence, isolated intended-edit effects, full haplotypes, additive validation or independent-study replication.'}

def main():
    if not __debug__:raise RuntimeError('Optimized execution is not permitted')
    p=argparse.ArgumentParser();p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();here=Path(__file__).resolve().parent;gate=here.parent/'cprime_gate'
    auth=json.loads(a.authorization.read_text());freeze=here/'freeze.json'
    if auth.get('permit_outcomes') is not True or auth.get('freeze_sha256')!=sha(freeze):raise ValueError('Missing matching outcome authorization')
    frozen=json.loads(freeze.read_text())
    for name,pin in frozen['files'].items():
        if sha(here/name)!=pin:raise ValueError('Frozen input changed: '+name)
    eligibility=json.loads((here/'eligibility.json').read_text())
    for name,pin in eligibility['input_sha256'].items():
        if sha(gate/name)!=pin:raise ValueError('Source input changed: '+name)
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    from interval_core import jsonable
    rows=read_rows(gate/'supp_Supplementary Data 8.xlsx',eligibility['states'],eligibility['cell_contexts'])
    result=analyze(rows,eligibility)
    result['freeze_sha256']=sha(freeze);result['source_workbook_sha256']=eligibility['input_sha256']['supp_Supplementary Data 8.xlsx']
    write(a.out/'result.json',result)
    write(a.out/'receipt.json',{'status':result['status'],'source_rows':len(rows),'result_sha256':sha(a.out/'result.json'),'freeze_sha256':sha(freeze)})
    print(json.dumps({'status':result['status'],'source_rows':len(rows),
                      'summaries':jsonable([s['summary'] for s in result.get('sensitivities',[])])}))
if __name__=='__main__':main()
