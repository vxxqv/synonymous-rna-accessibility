"""Invented inputs only. Never opens the published workbook or result."""
if not __debug__:raise RuntimeError('Optimized Python prohibited')
from pathlib import Path
from fractions import Fraction as Q
from copy import deepcopy
import random,json,zipfile,xml.etree.ElementTree as E,hashlib
import openpyxl
import reference as r
from compare_and_run import compare

p=Path(__file__).parent/'synthetic_v1'
if p.exists():raise FileExistsError(p)
p.mkdir(parents=True)
eligible={'states':{s:{'var_id':'g'+str(1+i//2),'syn_var':'m'+str(i),'exon_number':17 if i<2 else 25} for i,s in enumerate(['s1','s2','s3','s4'])},
          'groups':{'g1':['s1','s2'],'g2':['s3','s4']},'cell_contexts':['K562','KBM7'],'exons':[17,25]}
w=openpyxl.Workbook();s=w.active;s.title='Supplementary_Table_3_Multi_exo'
s.append(['combined_var_id','var_id','hgvsc','hgvsp','exon_number','syn_var','cell_line','RPM_DNA','RPM_RNA','cprime_score'])
for state,meta in eligible['states'].items():
 for cell in eligible['cell_contexts']:
  base=(1 if state=='s1' else 2) if cell=='K562' else (2 if state=='s1' else 1)
  for j in range(3):
   q=base if state in ['s1','s2'] else j+1+(state=='s4')
   s.append([state,meta['var_id'],'n','p',meta['exon_number'],meta['syn_var'],cell,10,10*q,q])
s.append(['OUTSIDE','ignored',None,None,17,'marker','K562',1,2,2]);source=p/'invented.xlsx';w.save(source)
z=zipfile.ZipFile(source);data={n:z.read(n) for n in z.namelist()};t=E.fromstring(data['xl/worksheets/sheet1.xml'])
for c in t.findall('.//m:c',r.N):
 if c.attrib['r']=='H26':c.find('m:v',r.N).text='UNSELECTED_POISON_MUST_NOT_BE_NUMERICALLY_PARSED'
 if c.attrib['r']=='I14':c.find('m:v',r.N).text='1.0000000000000001'
 if c.attrib['r']=='J14':c.find('m:v',r.N).text='0.10000000000000001'
data['xl/worksheets/sheet1.xml']=E.tostring(t)
with zipfile.ZipFile(source,'w') as out:
 for n,b in data.items():out.writestr(n,b)
rows,roundtrip=r.read_selected(source,eligible);assert len(rows)==24 and len(roundtrip)>=2
out=r.reconstruct(rows,eligible);assert out['status']=='RECONSTRUCTED'
assert out['sensitivities'][0]['exon_shared_order_results'][0]['status']=='infeasible'
assert out['sensitivities'][0]['exon_shared_order_results'][1]['status']=='feasible'
assert compare(r.ready(out),r.ready(out),eligible)['source_rows']==24
shuffled=deepcopy(rows);random.Random(9).shuffle(shuffled)
again=r.reconstruct(shuffled,eligible)
for x,y in zip(out['sensitivities'],again['sensitivities']):assert x==y

boxes={('a','c'):{'lower':Q(1),'upper':Q(10),'median':Q(5),'complete':True},
       ('b','c'):{'lower':Q(4),'upper':Q(6),'median':Q(5),'complete':True}}
b=r.group_bound('g','c',['a','b'],boxes,Q(1));assert b['minimum_fold']==1 and b['maximum_fold']==6
boxes['a','c'].update(lower=Q(0),upper=None,median=None,complete=False)
b=r.group_bound('g','c',['a','b'],boxes,Q(1));assert b['minimum_fold']==1 and b['maximum_fold'] is None and not b['maximum_attained']
bad=deepcopy(rows);bad[0].update(usable=False,ratio=Q(1),source_score=Q(-1),score_difference=Q(-2),source_score_consistent=False)
assert r.reconstruct(bad,eligible)['status']=='SOURCE_REPRODUCTION_FAILED'
missing=deepcopy(rows);missing[0].update(usable=False,missing_reason=['RNA'],RNA=None,ratio=None,score_difference=None,source_score_consistent=None)
mo=r.reconstruct(missing,eligible);assert not next(s for s in mo['state_cells'] if s['state']==missing[0]['state'] and s['cell']==missing[0]['cell'])['complete']
assert compare(r.ready(mo),r.ready(mo),eligible)['state_cells']==8
counter=0
from itertools import product
for vals in product([1,2,3],repeat=4):
 states={('a','K562'):{'lower':Q(vals[0]),'upper':Q(vals[0])},('b','K562'):{'lower':Q(vals[1]),'upper':Q(vals[1])},
         ('a','KBM7'):{'lower':Q(vals[2]),'upper':Q(vals[2])},('b','KBM7'):{'lower':Q(vals[3]),'upper':Q(vals[3])}}
 e={'states':{'a':{'exon_number':17},'b':{'exon_number':17}},'cell_contexts':['K562','KBM7']}
 result=r.shared_order(17,e,states,Q(1));opposite=(vals[0]-vals[1])*(vals[2]-vals[3])<0
 assert (result['status']=='infeasible')==opposite;counter+=1
receipt={'status':'PASS_SYNTHETIC_ONLY','real_workbook_or_result_opened':False,'synthetic_paired_rows':len(rows),
         'raw_decimal_roundtrip_differences':len(roundtrip),'two_context_order_cases':counter,
         'checks':['unselected_poison_ignored','within_row_pairs_preserved','exact_decimal_source','all_gamma_outputs',
                   'row_permutation','distinct_marker_maximum','unknown_state','nonpositive_score_reproduction_failure',
                   'feasible_positive_response','opposite_order_cycle','complete_record_comparator'],
         'files':{f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in ['reference.py','compare_and_run.py','test_synthetic.py']}}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
