"""Known-answer pre-outcome checks; all examples are invented."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,permutations
from tempfile import TemporaryDirectory
import random,json,zipfile,hashlib
from interval_core import marker_bounds,common_order,jsonable
from run_pilot import read_rows,analyze

def check(condition,message):
    if not condition:raise AssertionError(message)

def main():
    here=Path(__file__).resolve().parent
    elementary=[(F(i),F(j)) for i in range(1,5) for j in range(i,5)]
    count=0
    for n in [2,3]:
        for ints in product(elementary,repeat=n):
            candidates=sorted(set(x for pair in ints for x in pair))
            values=[max(ys)/min(ys) for ys in product(*[[x for x in candidates if l<=x<=u] for l,u in ints])]
            b=marker_bounds(ints)
            check(b['minimum']==min(values) and b['maximum']==max(values),'Extrema mismatch')
            count+=1
    b=marker_bounds([(F(1),F(10)),(F(4),F(6))])
    check(b['maximum']==6 and b['minimum']==1,'Same-marker extrema bug')
    check(marker_bounds([(0,None)])['maximum']==1,'Singleton bug')
    check(marker_bounds([(0,None),(2,3)])['maximum'] is None,'Unknown upper')
    check(marker_bounds([(0,None),(1,2),(4,5)])['minimum']==2,'Unknown hid existing separation')
    rng=random.Random(271828);oracle_checks=0;cyclic=0
    for _ in range(500):
        fields={c:{v:rng.choice(elementary) for v in ['a','b','c']} for c in ['x','y']}
        result=common_order(fields)
        possible=[]
        for order in permutations(['a','b','c']):
            okay=True
            for d in fields.values():
                val=F(0)
                for v in order:
                    val=max(val,d[v][0]);okay &= val<=d[v][1]
            if okay:possible.append(order)
        check((result['status']=='feasible')==bool(possible),'Shared order oracle mismatch')
        if result['status']=='feasible':
            for c,d in fields.items():
                prev=F(0)
                for v in result['latent_order']:
                    val=result['responses'][c][v]
                    check(d[v][0]<=val<=d[v][1] and val>=prev,'Invalid witness');prev=val
        else:
            cyclic+=1
            check(result['cycle'][0]==result['cycle'][-1],'Cycle not closed')
            check(all(e['origins'] and all(z['upper']<z['lower'] for z in e['origins']) for e in result['cycle_edges']),'False cycle edge')
        oracle_checks+=1
    check(common_order({'x':{'a':(0,None),'b':(0,None)},'y':{'a':(0,None),'b':(0,None)}})['status']=='feasible','All unknown')
    check(common_order({'x':{'a':(1,2),'b':(2,3)},'y':{'a':(2,3),'b':(1,2)}})['status']=='feasible','Endpoint tie')
    fields={'x':{'a':(1,1),'b':(2,2)},'y':{'a':(2,2),'b':(1,1)}}
    check(common_order(fields)['status']=='infeasible','Known reversal')
    eligibility={'states':{v:{'var_id':'edit','syn_var':v,'exon_number':1} for v in ['a','b']},
                 'groups':{'edit':['a','b']},'cell_contexts':['x','y'],'exons':[1]}
    rows=[]
    for v,c in product(['a','b'],['x','y']):
        for r in range(3):
            q=F(1 if v=='a' else 3)
            rows.append({'source_row':len(rows)+2,'state':v,'cell':c,'intended':'edit','marker':v,'exon':1,
                         'DNA':F(2),'RNA':2*q,'source_score':q,'ratio':q,'score_difference':F(0),
                         'source_score_consistent':True,'usable':True,'missing_reason':None})
    result=analyze(rows,eligibility)
    check(result['sensitivities'][0]['summary']['minimum_at_least_twofold']==2,'Positive group fixture')
    shuffled=list(reversed(rows));other=analyze(shuffled,eligibility)
    check(result['sensitivities']==other['sensitivities'],'Unordered source rows changed results')
    missing=[dict(r) for r in rows];missing[0].update({'usable':False,'ratio':None})
    check(analyze(missing,eligibility)['sensitivities'][0]['summary']['complete_group_cells']==1,'Missing group not propagated')
    bad=[dict(r) for r in rows];bad[0]['source_score_consistent']=False
    check(analyze(bad,eligibility)['status']=='SOURCE_REPRODUCTION_FAILED','Score failure hidden')
    # Minimal valid XLSX with deliberately poisoned unselected cells and rows.
    mainns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    def c(col,row,value,numeric=False):
        return f'<c r="{col}{row}"><v>{value}</v></c>' if numeric else f'<c r="{col}{row}" t="inlineStr"><is><t>{value}</t></is></c>'
    headers=['combined_var_id','var_id','unused','unused','exon_number','syn_var','cell_line','RPM_DNA','RPM_RNA','cprime_score','POISON']
    xmlrows=['<row r="1">'+''.join(c(chr(65+i),1,v) for i,v in enumerate(headers))+'</row>']
    for i,r in enumerate(rows,2):
        vals=[r['state'],'edit','unused','unused',1,r['marker'],r['cell'],2,int(r['RNA']),int(r['ratio']),'DO_NOT_DECODE']
        xmlrows.append(f'<row r="{i}">'+''.join(c(chr(65+j),i,v,j in [4,7,8,9]) for j,v in enumerate(vals))+'</row>')
    xmlrows.append('<row r="99">'+c('A',99,'NOT_ELIGIBLE')+c('H',99,'POISON',True)+'</row>')
    with TemporaryDirectory(dir=here,prefix='synthetic_') as tmp:
        p=Path(tmp)/'fixture.xlsx'
        with zipfile.ZipFile(p,'w') as z:
            z.writestr('xl/workbook.xml',f'<workbook xmlns="{mainns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Supplementary_Table_3_Multi_exo" r:id="rId1"/></sheets></workbook>')
            z.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
            z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{mainns}"><sheetData>'+''.join(xmlrows)+'</sheetData></worksheet>')
        parsed=read_rows(p,eligibility['states'],eligibility['cell_contexts'])
        check(len(parsed)==12 and all(r['ratio']==r['source_score'] for r in parsed),'Selective reader')
        check(analyze(parsed,eligibility)['sensitivities']==result['sensitivities'],'Reader result mismatch')
    receipt={'status':'PASS','finite_bound_oracles':count,'shared_order_permutation_oracles':oracle_checks,
             'cyclic_fixtures':cyclic,'known_positive_negative_tie_missing_cases':True,
             'selective_reader_poison_checks':True,'biological_outcomes_read':False,
             'code_sha256':{n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in ['interval_core.py','run_pilot.py','test_synthetic.py']}}
    (here/'synthetic_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':
    if not __debug__:raise RuntimeError('Optimized execution not allowed')
    main()
