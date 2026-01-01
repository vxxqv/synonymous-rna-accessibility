"""Synthetic-only tests. This file never opens the observed Gill workbook."""
import gzip, hashlib, importlib.util, json, math, shutil, subprocess, sys
from itertools import product
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('gill_multiclass',HERE/'gill_multiclass.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def item(last,cid,part,stratum):
    keys=sorted(''.join(v) for v in product('AC',['AAA','GAA'],['AAA','AAC'],[last]))
    edges=sum(sum(a!=b for a,b in zip(x,y))==1 for i,x in enumerate(keys) for y in keys[:i])
    return {'class_id':cid,'observed_keys':keys,'expected_vertices':8,'complete':True,'expected_component_sizes':[8],'expected_edges':edges,'codon_main_effect_rank':4,'codon_residual_dimensions':4,'codon_factors_all_cliques':True,'candidate_reproducibility_eligible':True,'candidate_partition':part,'reference_protein':False,'central_amino_acid_stratum_157_158':stratum}

def workbook(path,classes,poison_d=False,poison_unselected=True,allzero=False):
    keys=sorted(k for x in classes for k in x['observed_keys']);rows={k:i+2 for i,k in enumerate(keys)}
    for x in classes:x['observed_workbook_rows']=[rows[k] for k in x['observed_keys']]
    root=ET.Element(m.TAG+'worksheet');data=ET.SubElement(root,m.TAG+'sheetData')
    def text(row,col,r,value):
        c=ET.SubElement(row,m.TAG+'c',r=col+str(r),t='inlineStr');s=ET.SubElement(c,m.TAG+'is');ET.SubElement(s,m.TAG+'t').text=value
    row=ET.SubElement(data,m.TAG+'row',r='1')
    for col,h in m.HEADERS.items():
        if poison_d and col=='D':c=ET.SubElement(row,m.TAG+'c',r='D1',t='s');ET.SubElement(c,m.TAG+'v').text='99999999'
        else:text(row,col,1,h)
    for i,k in enumerate(keys):
        rn=i+2;row=ET.SubElement(data,m.TAG+'row',r=str(rn));text(row,'A',rn,k)
        for col in 'BCDEFG':
            c=ET.SubElement(row,m.TAG+'c',r=col+str(rn),t='n')
            if poison_d and col=='D':c.set('t','s');value='99999999'
            else:
                j='BCDEFG'.index(col);value=0 if allzero else (64+17*i+11*j)%137
            ET.SubElement(c,m.TAG+'v').text=str(value)
    if poison_unselected:
        rn=len(keys)+2;row=ET.SubElement(data,m.TAG+'row',r=str(rn));text(row,'A',rn,'TTTTTTTT')
        for col in 'BCDEFG':c=ET.SubElement(row,m.TAG+'c',r=col+str(rn),t='s');ET.SubElement(c,m.TAG+'v').text='99999999'
    with ZipFile(path,'w') as z:z.writestr('xl/worksheets/sheet1.xml',ET.tostring(root))
    return keys

def fixture(tmp_path,poison_d=False,allzero=False):
    classes=[item('T','synthetic_dev','development','AA'),item('G','synthetic_val','validation','CC')]
    source=tmp_path/'toy.xlsx';keys=workbook(source,classes,poison_d=poison_d,allzero=allzero)
    raw=(json.dumps({'classes':classes})+'\n').encode();inv=tmp_path/'inventory.json.gz';inv.write_bytes(gzip.compress(raw,mtime=0))
    p={'id':'GILL_MULTICLASS_REPRODUCIBILITY_V1','status':'frozen_before_nonreference_RNA_outcomes','inventory_sha256':hashlib.sha256(raw).hexdigest(),'source_workbook_sha256':m.digest(source),'class_ids_by_partition':{'development':['synthetic_dev'],'validation':['synthetic_val'],'reference_stratum_quarantine':[]},'primary_class_count':1,'primary_vertices':8,'development_class_count':1,'development_vertices':8,'pseudocounts':[.25,.5,1.0],'primary_pseudocount':.5,'training_RNA_columns':['B','C'],'training_DNA_columns':['E','F'],'validation_RNA_column':'D','validation_DNA_column':'G','primary_descriptive_reproducibility_target':.5,'uncertainty':'synthetic test, no CI','response_scale_limit':'synthetic test'}
    pp=tmp_path/'protocol.json';m.write_json(pp,p)
    args=SimpleNamespace(protocol=pp,inventory=inv,workbook=source,output=tmp_path/'training',expected_protocol_sha256=m.digest(pp))
    return args,classes,keys

def validated(tmp_path,allzero=False):
    args,classes,keys=fixture(tmp_path,allzero=allzero);m.train(args)
    args.training=args.output;args.expected_training_receipt_sha256=m.digest(args.training/'receipt.json');args.output=tmp_path/'validation'
    m.validate(args);return args,classes,keys

def test_train_poisoned_D_header_rows_and_unselected(tmp_path):
    a,classes,keys=fixture(tmp_path,poison_d=True);m.train(a)
    summary=m.read_json(a.output/'training_summary.json')
    assert 'D' not in summary['decoded_columns'];assert summary['vertices']==16
    counts=m.read_gz(a.output/'counts.json.gz');assert len(counts)==16 and set(counts[keys[0]])==set(m.TRAIN)

def test_selective_validate_skips_other_RNA(tmp_path):
    classes=[item('T','x','validation','AA')];p=tmp_path/'toy.xlsx';keys=workbook(p,classes)
    with ZipFile(p) as z:xml=z.read('xl/worksheets/sheet1.xml')
    x=ET.fromstring(xml)
    for cell in x.iter(m.TAG+'c'):
        if re_col(cell.attrib['r']) in ('B','C','E','F'):
            cell.clear();cell.attrib.update(r='B999',t='s');ET.SubElement(cell,m.TAG+'v').text='9999999'
    # Independent direct reader test with normal A/D/G and forbidden invalid cells.
    for row in x.iter(m.TAG+'row'):
        for i,c in enumerate(row):
            if c.attrib.get('r')=='B999':c.set('r',('B','C','E','F')[i%4]+row.attrib['r'])
    # Preserve unique column identities using a simpler hand-built three-column source.
    root=ET.Element(m.TAG+'worksheet');sd=ET.SubElement(root,m.TAG+'sheetData')
    for rn,key in [(1,'Sequence')]+[(i+2,k) for i,k in enumerate(keys)]:
        row=ET.SubElement(sd,m.TAG+'row',r=str(rn))
        for col,value in [('A',key),('D',m.HEADERS['D'] if rn==1 else '0'),('G',m.HEADERS['G'] if rn==1 else '5')]:
            c=ET.SubElement(row,m.TAG+'c',r=col+str(rn),t='inlineStr');ET.SubElement(ET.SubElement(c,m.TAG+'is'),m.TAG+'t').text=value
        c=ET.SubElement(row,m.TAG+'c',r='B'+str(rn),t='s');ET.SubElement(c,m.TAG+'v').text='99999999'
    with ZipFile(p,'w') as z:z.writestr('xl/worksheets/sheet1.xml',ET.tostring(root))
    values,audit=m.read_selected(p,keys,m.TEST);assert len(values)==8 and 'B' not in audit

def re_col(r):return ''.join(x for x in r if not x.isdigit())

@pytest.mark.parametrize('x',['-1','0.1','NaN','Infinity','bad',None])
def test_invalid_counts(x):
    with pytest.raises(ValueError):m.integer(x)

def test_count_zero_ties_undefined_and_complete_pipeline(tmp_path):
    a,_,_=validated(tmp_path,allzero=True);v=m.read_json(a.output/'validation.json')
    assert v['primary_result']['cosine'] is None and v['primary_result']['status']=='information_limited'
    assert v['primary_result']['zero_energy_classes']=={'training':1,'validation':1}
    metrics=m.read_gz(a.output/'class_metrics.json.gz')['0.5']['validation'][0]
    assert metrics['RNA_rank_correlations']==[None]*3
    assert len(metrics['decomposition_terms'])==45
    assert all(x['maximum']==0 for x in metrics['scenario_barriers'])

def test_projection_independent_lstsq_and_invariance():
    g=item('T','x','validation','AA');d=m.design(g);rng=np.random.default_rng(17)
    keys=d['keys'];cols=[[1.]*8]
    for a,b in m.SEGMENTS:
        levels=sorted({k[a:b] for k in keys})
        cols += [[float(k[a:b]==v) for k in keys] for v in levels[1:]]
    X=np.array(cols).T
    for _ in range(10):
        f=rng.normal(size=8);p,q=m.project(f,d);ols=X@np.linalg.lstsq(X,f,rcond=None)[0]
        assert np.allclose(p,ols,atol=1e-12);assert np.allclose(X.T@q,0,atol=1e-12)
        assert np.allclose(m.project(f+13.5,d)[1],q,atol=1e-12)
        assert np.allclose(m.project(q,d)[1],q,atol=1e-12)
    assert np.array_equal(m.project(np.ones(8)*.123456789,d)[1],np.zeros(8))

def test_maximum_barrier_independent_floyd():
    rng=np.random.default_rng(41)
    for n in range(2,9):
        for _ in range(12):
            adj=[set() for _ in range(n)]
            for i in range(1,n):adj[i].add(i-1);adj[i-1].add(i)
            for i in range(n):
                for j in range(i):
                    if rng.random()<.3:adj[i].add(j);adj[j].add(i)
            f=rng.integers(-3,4,size=n).astype(float);cap=np.full((n,n),-np.inf)
            for i in range(n):
                cap[i,i]=f[i]
                for j in adj[i]:cap[i,j]=min(f[i],f[j])
            for k in range(n):
                for i in range(n):
                    for j in range(n):cap[i,j]=max(cap[i,j],min(cap[i,k],cap[k,j]))
            exact=max(min(f[i],f[j])-cap[i,j] for i in range(n) for j in range(i))
            assert m.maximum_barrier(f,adj)['maximum']==exact

def test_decomposition_shared_DNA_and_arbitrary_vectors():
    rng=np.random.default_rng(81)
    for _ in range(30):
        a,b,c,d=rng.normal(size=(4,8));terms,err=m.quad_terms(a,b,c,d)
        assert math.isclose(sum(terms[:4]),terms[4],abs_tol=1e-12)
    zero=np.zeros(8);dna=np.arange(8.)
    t,_=m.quad_terms(zero,dna,zero,dna)
    assert t[0]==0 and t[1]==t[4] and t[2]==t[3]==0

@pytest.mark.parametrize('fault',['receipt','prediction','counts','protocol','source'])
def test_validation_guard_before_D(tmp_path,monkeypatch,fault):
    a,_,_=fixture(tmp_path);m.train(a);a.training=a.output;a.expected_training_receipt_sha256=m.digest(a.training/'receipt.json');a.output=tmp_path/'validation'
    if fault=='receipt':a.expected_training_receipt_sha256='0'*64
    elif fault in ['prediction','counts']:
        f=a.training/('predictions.json.gz' if fault=='prediction' else 'counts.json.gz');f.write_bytes(f.read_bytes()+b'X')
    elif fault=='protocol':a.expected_protocol_sha256='0'*64
    else:a.workbook.write_bytes(a.workbook.read_bytes()+b'X')
    def forbidden(*a,**kw):raise AssertionError('D reader called before guard')
    monkeypatch.setattr(m,'read_selected',forbidden)
    with pytest.raises(ValueError):m.validate(a)
    assert not a.output.exists()

def test_successful_pipeline_preserves_all_scenarios(tmp_path):
    a,_,_=validated(tmp_path);v=m.read_json(a.output/'validation.json')
    assert set(v['all_summaries'])=={'0.25','0.5','1.0'}
    for c,parts in v['all_summaries'].items():
        for part,s in parts.items():
            assert s['class_count']==1 and len(s['decomposition_aggregate_terms'])==45
            assert len(s['scenario_residual_pooled_cosines'])==45
            assert len(s['leave_one_stratum_out'])==1
            assert s['leave_one_stratum_out'][0]['remaining_cosine'] is None
    assert v['exposure_decoded_columns']['D']==17 and set(v['exposure_decoded_columns'])=={'A','D','G'}

def test_portable_CLI_without_git_checkout(tmp_path):
    a,_,_=fixture(tmp_path);script=tmp_path/'portable.py';shutil.copyfile(HERE/'gill_multiclass.py',script)
    base=[sys.executable,str(script)];common=['--protocol',str(a.protocol),'--inventory',str(a.inventory),'--workbook',str(a.workbook),'--expected-protocol-sha256',a.expected_protocol_sha256]
    subprocess.run(base+['train']+common+['--output',str(a.output)],check=True,capture_output=True)
    receipt=m.digest(a.output/'receipt.json')
    subprocess.run(base+['validate']+common+['--output',str(tmp_path/'validation'),'--training',str(a.output),'--expected-training-receipt-sha256',receipt],check=True,capture_output=True)
    assert m.read_json(tmp_path/'validation/receipt.json')['D_opened_after_training_gate']
    failed=subprocess.run([sys.executable,'-O',str(script),'--help'],capture_output=True)
    assert failed.returncode!=0

def test_raw_integer_rank_preserved_above_float_resolution():
    a=[2**60,2**60+1,2**60+2];assert list(m.ranks(a))==[0.,1.,2.]
    assert m.spearman(a,list(reversed(a)))==-1

@pytest.mark.parametrize('fault',['missing','duplicate','formula'])
def test_selected_membership_and_formula_failures(tmp_path,fault):
    a,classes,keys=fixture(tmp_path)
    with ZipFile(a.workbook) as z:x=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    sd=x.find(m.TAG+'sheetData')
    if fault=='missing':sd.remove(list(sd)[1])
    elif fault=='duplicate':list(sd)[2].find(m.TAG+'c').find('.//'+m.TAG+'t').text=keys[0]
    else:
        c=next(c for c in list(sd)[1] if c.attrib['r'].startswith('B'));ET.SubElement(c,m.TAG+'f').text='1+2'
    with ZipFile(a.workbook,'w') as z:z.writestr('xl/worksheets/sheet1.xml',ET.tostring(x))
    with pytest.raises(ValueError):m.read_selected(a.workbook,keys,m.TRAIN)

def test_connected_nonclique_codon_factor_projection():
    keys=sorted(''.join(v) for v in product('AC',['TTA','TTG','CTT','CTC','CTA','CTG'],['AAA','AAC'],['T']))
    edges=sum(sum(a!=b for a,b in zip(x,y))==1 for i,x in enumerate(keys) for y in keys[:i])
    x={'observed_keys':keys,'expected_vertices':24,'complete':True,'expected_component_sizes':[24],'expected_edges':edges,'codon_main_effect_rank':8,'codon_residual_dimensions':16}
    d=m.design(x);f=np.array([i%7 for i in range(24)],dtype=float);p,q=m.project(f,d)
    cols=[[1.]*24]
    for a,b in m.SEGMENTS:
        dom=sorted({k[a:b] for k in keys});cols += [[float(k[a:b]==v) for k in keys] for v in dom[1:]]
    X=np.array(cols).T
    assert np.linalg.matrix_rank(X)==8
    assert np.allclose(p,X@np.linalg.lstsq(X,f,rcond=None)[0],atol=1e-12)
    assert math.isclose(m.dot(p,q),0,abs_tol=1e-11)

def test_additive_clique_barrier_precision_stress():
    g=item('T','x','validation','AA');d=m.design(g);rng=np.random.default_rng(943)
    for scale in [1e-9,1.0,1e9]:
        for offset in [0.0,1e12,-1e12]:
            for _ in range(20):
                f=np.full(8,offset)
                for groups in d['groups']:
                    coefficients=rng.normal(size=len(groups))*scale
                    for idx,c in zip(groups,coefficients):f[idx]+=c
                fitted,_=m.project(f,d);b=m.maximum_barrier(fitted,d['adjacency'])
                assert not b['numerically_resolved_positive']
                assert b['maximum']>=0 and b['precision_tolerance']>0

def test_tiny_positive_is_retained_and_flagged_separately():
    d=m.design(item('T','x','validation','AA'));f=np.zeros(8);delta=16*np.finfo(float).eps
    f[0]=f[7]=delta;b=m.maximum_barrier(f,d['adjacency'])
    assert b['maximum']==delta and b['unresolved_tiny_positive']
    assert not b['numerically_resolved_positive']
    f[0]=f[7]=1.0;b=m.maximum_barrier(f,d['adjacency'])
    assert b['maximum']==1 and b['numerically_resolved_positive']
    assert not b['unresolved_tiny_positive']

def test_barrier_summary_precision_counts_and_shape_label(tmp_path):
    a,_,_=validated(tmp_path);result=m.read_json(a.output/'validation.json')
    assert 'shape' in result['centered_loss_interpretation']
    assert 'not out-of-sample level' in result['centered_loss_interpretation']
    assert 'not measurement uncertainty' in result['barrier_precision_rule']
    for parts in result['all_summaries'].values():
        for summary in parts.values():
            for architecture in summary['barriers_by_architecture'].values():
                for name in ['scenario_barriers','scenario_codon_barriers']:
                    for r in architecture[name]:
                        assert r['raw_positive_classes']==r['numerically_resolved_positive_classes']+r['unresolved_tiny_positive_classes']
