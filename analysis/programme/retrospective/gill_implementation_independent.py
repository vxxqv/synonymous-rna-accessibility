# Archived from the pre-outcome synthetic review; original SHA256 72eeb8aa83770da643badb13e0e3e052c862570a3977a7838ef2a0201e4a7e3a
"""Independent toy-only Gill runner audit. Never opens the observed workbook."""
from pathlib import Path
from fractions import Fraction
from collections import deque
from itertools import product, combinations
import hashlib
import importlib.util
import json
import math
import random
import tempfile
from unittest.mock import patch
from zipfile import ZipFile,ZIP_DEFLATED
import xml.etree.ElementTree as ET
import numpy as np

if not __debug__:
    raise RuntimeError('Audit checks require assertions enabled')
REPO=Path(__file__).resolve().parents[3]
SOURCE=REPO/'analysis'/'programme'/'retrospective'/'gill_reference_pilot.py'
specmod=importlib.util.spec_from_file_location('reviewed_gill_runner',SOURCE)
m=importlib.util.module_from_spec(specmod);specmod.loader.exec_module(m)
before=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
assert before=='f57f084c85fbfae3b27888e55284ff41e10fb178d719351fe94c40d17c9c2b83', 'Reviewed runner hash changed'
keys=[''.join(v) for v in product('CT','G','T','ACGT','C','C','ACGT','T')]
vertices=list(product(range(2),range(4),range(4)))
spec=m.design(keys)
graph=[[j for j,y in enumerate(vertices) if sum(a!=b for a,b in zip(x,y))==1] for x in vertices]
assert [list(a) for a in spec['graph']]==graph
rng=random.Random(2026102504)
counts={'factorial_fields':0,'isotonic_cases':0,'geometry_fields':0,'score_checks':0,
        'invalid_count_checks':0,'selective_reader_checks':0,'end_to_end_toy_runs':0,'freeze_failures':0}

def projection_reference(values,order):
    # Inclusion-exclusion of conditional arithmetic means, independent of tensor matrices.
    out=[]
    for vertex in vertices:
        total=Fraction(0)
        for selected in combinations(range(3),order):
            for size in range(order+1):
                for subset in combinations(selected,size):
                    matching=[Fraction(int(y)) for x,y in zip(vertices,values) if all(x[k]==vertex[k] for k in subset)]
                    total += (-1)**(order-size)*sum(matching)/len(matching)
        out.append(float(total))
    return np.array(out)

for trial in range(60):
    z=np.array([rng.randrange(-20,21) for _ in keys],dtype=float)
    parts=[projection_reference(z,o) for o in range(4)]
    for order in range(4):
        assert np.array_equal(m.matvec(spec['projectors'][order],z),parts[order])
    predictions,_=m.fit_models(z,spec,3)
    target=parts[0]+parts[1]+parts[2]/2
    assert np.array_equal(predictions['pairwise_half'],target-target[3])
    assert np.array_equal(predictions['additive'],parts[0]+parts[1]-(parts[0]+parts[1])[3])
    counts['factorial_fields']+=1
adj=np.zeros((32,32))
for u,neighbors in enumerate(graph):adj[u,neighbors]=1
L=7*np.eye(32)-adj
eigs,basis=np.linalg.eigh(L)
reference_heat=(basis*np.exp(-math.log(2)/4*eigs))@basis.T
assert np.max(abs(reference_heat-spec['heat']))<2e-15
assert abs(float(np.trace(spec['heat']))-25/4*(1+1/math.sqrt(2)))<1e-13

def pava_reference(x,y):
    knots=sorted(set(x));weights=[x.count(k) for k in knots]
    totals=[sum(Fraction(v) for xx,v in zip(x,y) if xx==k) for k in knots]
    return knots,[max(min(sum(totals[a:b+1])/sum(weights[a:b+1]) for b in range(i,len(knots))) for a in range(i+1)) for i in range(len(knots))]

for trial in range(150):
    n=rng.randrange(1,15)
    x=[rng.randrange(-3,4) for _ in range(n)]
    y=[rng.randrange(-10,11) for _ in range(n)]
    knots,fit=pava_reference(x,y);got=m.isotonic_fit(x,y)
    assert got['input_knots']==knots
    assert np.max(abs(np.array(got['fitted_levels'])-np.array(list(map(float,fit)))))<1e-12
    assert got['constant_blocks']==len(set(fit))
    assert list(m.isotonic_predict(got,[min(x)-1,max(x)+1]))==[float(fit[0]),float(fit[-1])]
    counts['isotonic_cases']+=1

def geometry_reference(z):
    capacities=np.full((32,32),np.nan);thresholds=[]
    for threshold in sorted(set(z),reverse=True):
        remaining={i for i,v in enumerate(z) if v>=threshold};retained=len(remaining);sizes=[]
        while remaining:
            u=min(remaining);remaining.remove(u);component={u};queue=deque([u])
            while queue:
                for v in graph[queue.popleft()]:
                    if v in remaining:remaining.remove(v);component.add(v);queue.append(v)
            for u in component:
                for v in component:
                    if np.isnan(capacities[u,v]):capacities[u,v]=threshold
            sizes.append(len(component))
        thresholds.append({'threshold':float(threshold),'retained':retained,'components':len(sizes),'component_sizes':sorted(sizes,reverse=True)})
    return capacities,thresholds

for trial in range(25):
    z=[rng.randrange(-5,6) for _ in keys]
    capacity,thresholds=geometry_reference(z);got=m.geometry(z,spec)
    assert got['critical_thresholds']==thresholds and got['pair_count']==496
    for row in got['pairs']:
        i,j=row['i'],row['j']
        assert row['capacity']==capacity[i,j] and row['barrier']==min(z[i],z[j])-capacity[i,j]
    counts['geometry_fields']+=1

for c in [.25,.5,1]:
    r1=[0]+list(range(1,32));r2=[i+3 for i in range(32)]
    d1=[100+i for i in range(32)];d2=[103+2*i for i in range(32)]
    got=m.score_counts([r1,r2],[d1,d2],c)
    ref=np.array([(math.log2(a+c)+math.log2(b+c)-math.log2(d+c)-math.log2(e+c))/2 for a,b,d,e in zip(r1,r2,d1,d2)])
    assert np.max(abs(got-ref))<2e-15
    counts['score_checks']+=32
for bad in [None,'nan','inf','-inf','-1','1.5','not-a-count','1e9999']:
    try:m.count_value(bad)
    except ValueError:counts['invalid_count_checks']+=1
    else:raise AssertionError(('invalid count accepted',bad))
assert m.count_value('0')==0 and m.count_value('1e2')==100

HERE=REPO/'results/programme/gill_reference_pilot/independent_audit/implementation_synthetic_replay'
HERE.mkdir(parents=True,exist_ok=False)
XMLNS='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
def toy_workbook(path,poison=False,missing=None,duplicate=False,bad_selected=None,bad_header=None):
    root=ET.Element('worksheet',xmlns=XMLNS);sheet=ET.SubElement(root,'sheetData')
    header=ET.SubElement(sheet,'row',r='1')
    def cell(row,address,text):
        c=ET.SubElement(row,'c',r=address,t='inlineStr');s=ET.SubElement(c,'is');ET.SubElement(s,'t').text=str(text)
    for col,text in m.HEADERS.items():cell(header,col+'1','POISON_D_HEADER' if poison and col=='D' else ('WRONG' if bad_header==col else text))
    data=[]
    for i,key in enumerate(keys):
        if key==missing:continue
        row=ET.SubElement(sheet,'row',r=str(i+2));cell(row,'A'+str(i+2),key)
        vals={'B':i,'C':i+3,'D':i+7,'E':100+i,'F':103+2*i,'G':110+3*i}
        for col,value in vals.items():
            cell(row,col+str(i+2),'UNOPENED_D_OUTCOME' if poison and col=='D' else (bad_selected[1] if bad_selected and col==bad_selected[0] and i==0 else value))
        data.append(vals)
    extra=ET.SubElement(sheet,'row',r='100');cell(extra,'A100','AAAAAAAA')
    for col in 'BCDEFG':
        node=ET.SubElement(extra,'c',r=col+'100');ET.SubElement(node,'f').text='UNRELATED_PROTEIN_OUTCOME';ET.SubElement(node,'v').text='INVALID'
    if duplicate:
        row=ET.SubElement(sheet,'row',r='101');cell(row,'A101',keys[0])
        for col in 'BCDEFG':cell(row,col+'101',100)
    with ZipFile(path,'w',ZIP_DEFLATED) as z:z.writestr('xl/worksheets/sheet1.xml',ET.tostring(root))
    return data

with tempfile.TemporaryDirectory(prefix='toy_only_',dir=HERE) as temporary:
    tmp=Path(temporary);book=tmp/'toy.xlsx';toy_workbook(book,poison=True)
    touched=[];original=m._cell_text
    def observe(cell,shared):
        touched.append(cell.attrib['r']);return original(cell,shared)
    with patch.object(m,'_cell_text',side_effect=observe):got=m.selected_counts(book,keys,'train')
    assert got['B'][0]==0 and len(got['B'])==32
    assert not any(a.startswith('D') for a in touched)
    assert not any(a.endswith('100') and not a.startswith('A') for a in touched)
    counts['selective_reader_checks']+=1
    for kwargs in [{'missing':keys[-1]},{'duplicate':True},{'bad_selected':('B','-1')},{'bad_selected':('C','0.5')},{'bad_header':'E'}]:
        toy_workbook(book,**kwargs)
        try:m.selected_counts(book,keys,'train')
        except ValueError:counts['selective_reader_checks']+=1
        else:raise AssertionError(('bad toy accepted',kwargs))
    toy_workbook(book)
    protocol=json.loads((REPO/'analysis/programme/retrospective/gill_reference_protocol_v1.json').read_text())
    eligibility=json.loads((REPO/'results/programme/retrospective/gill_amplicon_pilot_v2.json').read_text())
    eligibility['inputs']['cds_processed.xlsx']['sha256']=m.digest(book)
    eligpath=tmp/'eligibility.json';m.write_json(eligpath,eligibility)
    protocol.update({'status':'frozen','keys':keys,'workbook_sha256':m.digest(book),'eligibility_sha256':m.digest(eligpath)})
    protocolpath=tmp/'protocol.json';m.write_json(protocolpath,protocol)
    train_dir=tmp/'training';m.train(protocolpath,eligpath,book,train_dir)
    counts['end_to_end_toy_runs']+=1
    saved=json.loads((train_dir/'predictions.json').read_text())
    assert saved['exposure']['D_cell_values_accessed'] is False
    for c in m.PSEUDOCOUNTS:
        assert set(saved['predictions'][str(c)]['models'])==set(m.MODELS)
    # A wrong prediction hash must fail before any workbook open/hash.
    with patch.object(m,'verify_protocol',side_effect=AssertionError('workbook stage reached before bad-hash stop')):
        try:m.validate(protocolpath,eligpath,book,train_dir/'predictions.json','0'*64,tmp/'bad_hash')
        except ValueError:counts['freeze_failures']+=1
        else:raise AssertionError('wrong prediction hash accepted')
    m.validate(protocolpath,eligpath,book,train_dir/'predictions.json',m.digest(train_dir/'predictions.json'),tmp/'validation')
    counts['end_to_end_toy_runs']+=1
    validation=json.loads((tmp/'validation'/'validation.json').read_text())
    assert validation['count_columns_accessed']==['D','G']
    for block in validation['results'].values():
        for metric in block['models'].values():
            assert abs(metric['rmse31']-metric['rmse32']*math.sqrt(32/31))<1e-12
            assert metric['errors'][saved['reference_index']]==0
    # Changing only G cannot influence any fitted model, although coverage retains G.
    train_counts=saved['training_counts'];first=m.score_counts([train_counts['B'],train_counts['C']],[train_counts['E'],train_counts['F']],.5)
    changed=dict(train_counts);changed['G']=[v+10000 for v in changed['G']]
    second=m.score_counts([changed['B'],changed['C']],[changed['E'],changed['F']],.5)
    assert np.array_equal(first,second)
    all_zero={name:{'rmse32':0.0} for name in m.MODELS}
    assert m.decision(all_zero)['passes_10_percent_hurdle'] is False

after=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
assert before==after,'source changed during independent run'
receipt={'status':'passed','scope':'Independent synthetic-only implementation audit; actual workbook never opened',
         'seed':2026102504,'checks':counts,'runner_sha256_before':before,'runner_sha256_after':after,
         'independent_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'metadata_protocol_only_actual_inputs':['gill_reference_protocol_v1.json','gill_amplicon_pilot_v2.json'],
         'observed_workbook_opened':False,'observed_RNA_values_read':False}
(HERE/'checks.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
