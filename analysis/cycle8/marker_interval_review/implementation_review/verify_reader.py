"""Synthetic-only selective-reader and orchestration challenge."""
if not __debug__:
    raise RuntimeError('Optimized execution is forbidden.')

from collections import Counter
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import copy
import hashlib
import itertools
import json
import os
import subprocess
import sys
import types
import zipfile
from xml.sax.saxutils import escape

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PILOT = HERE.parents[1] / 'cprime_pilot'
NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
RUNNER_PIN = 'bfba04ed71fe695978327680e54fca45b197da42809db0632fe935a4e95d3cd7'
CORE_PIN = '9d0bd0f3363e158e4df1cb028c7922b0bf94284bc6169138f894b072871c4cf9'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name, path, pin):
    assert sha(path) == pin
    result = types.ModuleType(name)
    result.__file__ = str(path)
    sys.modules[name] = result
    exec(compile(path.read_bytes(), str(path), 'exec'), result.__dict__)
    return result


def text(value):
    if isinstance(value,F): return str(Decimal(value.numerator)/Decimal(value.denominator))
    return str(value)


def worksheet(path, eligible, overrides=None, omit=None, extra=False, shared_probe=False):
    overrides = overrides or {}
    omit = set(omit or ())
    headers = {'A':'combined_var_id','B':'var_id','E':'exon_number','F':'syn_var','G':'cell_line','H':'RPM_DNA','I':'RPM_RNA','J':'cprime_score'}
    def cell(col,n,value,kind='inline'):
        if value is None:return ''
        if kind == 'formula':return f'<c r="{col}{n}"><f>PROTECTED_FORMULA()</f><v>POISON</v></c>'
        if kind == 'numeric':return f'<c r="{col}{n}"><v>{escape(text(value))}</v></c>'
        if kind == 'shared':return f'<c r="{col}{n}" t="s"><v>{value}</v></c>'
        return f'<c r="{col}{n}" t="inlineStr"><is><t>{escape(text(value))}</t></is></c>'
    rows = ['<row r="1">'+''.join(cell(k,1,v) for k,v in headers.items())+cell('K',1,'ignored','formula')+'</row>']
    raw=[]
    for index,(state,context,repeat) in enumerate(itertools.product(sorted(eligible['states']),eligible['cell_contexts'],range(3)),2):
        info=eligible['states'][state]
        n=sorted(eligible['states']).index(state)
        q=F(n+1,1)+F(repeat,10) if context==eligible['cell_contexts'][0] else F(500-n,1)+F(repeat,10)
        dna=(F(1),F(2),F(4))[repeat]
        values={'A':state,'B':info['var_id'],'E':info['exon_number'],'F':info['syn_var'],'G':context,'H':dna,'I':dna*q,'J':q}
        values.update(overrides.get(index,{}))
        if index not in omit:
            content=''
            for col,value in values.items():
                if shared_probe and col=='A' and state==sorted(eligible['states'])[0]:content+=cell(col,index,1,'shared')
                else:content+=cell(col,index,value,'numeric' if col in {'E','H','I','J'} else 'inline')
            rows.append(f'<row r="{index}">'+content+cell('K',index,'ignored','formula')+'</row>')
            raw.append((index,values))
    # An unselected body row contains poisoned outcome formulas.
    rows.append('<row r="9999">'+cell('A',9999,'NOT_SELECTED')+cell('H',9999,'ignored','formula')+cell('J',9999,'ignored','formula')+'</row>')
    if extra:
        index,values=raw[0]
        rows.append('<row r="9998">'+''.join(cell(col,9998,v,'numeric' if col in {'E','H','I','J'} else 'inline') for col,v in values.items())+'</row>')
    with zipfile.ZipFile(path,'w') as z:
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{NS}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Supplementary_Table_3_Multi_exo" r:id="rId1"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{NS}"><sheetData>'+''.join(rows)+'</sheetData></worksheet>')
        if shared_probe:
            z.writestr('xl/sharedStrings.xml',f'<sst xmlns="{NS}"><si><t>PROTECTED_UNSELECTED</t></si><si><t>{escape(sorted(eligible["states"])[0])}</t></si></sst>')
    return raw


def direct_bounds(pairs):
    positive={x for p in pairs for x in p if x is not None and x>0}
    positive |= {min(positive)/2} if positive else {F(1)}
    choices=[[x for x in positive if x>=l and (u is None or x<=u)] for l,u in pairs]
    vals=[max(v)/min(v) for v in itertools.product(*choices)]
    return min(vals), None if len(pairs)>1 and any(l==0 or u is None for l,u in pairs) else max(vals)


def cycle_by_transitive_closure(fields):
    vertices=sorted(next(iter(fields.values())))
    n=len(vertices)
    reach=[[any(d[vertices[i]][1] is not None and d[vertices[i]][1]<d[vertices[j]][0] for d in fields.values()) if i!=j else False for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):reach[i][j] |= reach[i][k] and reach[k][j]
    return any(reach[i][i] for i in range(n))


def main():
    core=module('interval_core',PILOT/'interval_core.py',CORE_PIN)
    run=module('fixture_runner',PILOT/'run_pilot.py',RUNNER_PIN)
    eligibility=json.loads((PILOT/'eligibility.json').read_text())
    assert len(eligibility['states'])==66 and len(eligibility['groups'])==30
    assert Counter(map(len,eligibility['groups'].values()))=={2:24,3:6}
    assert eligibility['exons']==[17,25,27,39,42]
    mini={'states':{'a':{'var_id':'edit','syn_var':'ma','exon_number':1},'b':{'var_id':'edit','syn_var':'mb','exon_number':1}},'groups':{'edit':['a','b']},'cell_contexts':['A','B'],'exons':[1]}
    tested=[]
    guards=[]
    with TemporaryDirectory(dir=HERE,prefix='invented_') as tmp:
        tmp=Path(tmp)
        f=tmp/'fixture.xlsx'
        worksheet(f,eligibility)
        rows=run.read_rows(f,eligibility['states'],eligibility['cell_contexts'])
        result=run.analyze(rows,eligibility)
        assert len(rows)==396 and len(result['state_cells'])==132 and len(result['sensitivities'])==4
        assert result['status']=='COMPLETE_DESCRIPTIVE_AUDIT'
        observed={(o['state'],o['cell']):(o['lower'],o['upper']) for o in result['state_cells']}
        for sensitivity in result['sensitivities']:
            gamma=sensitivity['gamma']
            assert len(sensitivity['group_cell_results'])==60 and len(sensitivity['exon_shared_order_results'])==5
            for record in sensitivity['group_cell_results']:
                pairs=[observed[v,record['cell']] for v in record['markers']]
                pairs=[(l/gamma,None if u is None else u*gamma) for l,u in pairs]
                assert (record['minimum_fold'],record['maximum_fold'])==direct_bounds(pairs)
            for record in sensitivity['exon_shared_order_results']:
                ids=[v for v in eligibility['states'] if int(eligibility['states'][v]['exon_number'])==record['exon']]
                fields={c:{v:(observed[v,c][0]/gamma,observed[v,c][1]*gamma) for v in ids} for c in eligibility['cell_contexts']}
                assert (record['status']=='infeasible')==cycle_by_transitive_closure(fields)
        assert run.analyze(list(reversed(rows)),eligibility)['sensitivities']==result['sensitivities']
        tested.append('full66state_132cell_396row_fixture_all4widenings_all60bounds_all5exons_and_row_permutation')
        # Metadata validity must be checked before trying poisoned outcome conversion.
        worksheet(f,mini,{2:{'G':'FOREIGN_CONTEXT','I':'POISON'}})
        try:run.read_rows(f,mini['states'],mini['cell_contexts'])
        except ValueError as e:assert 'before outcome decoding' in str(e)
        else:raise AssertionError('Foreign context accepted')
        tested.append('identity_before_outcome_conversion')
        # Deliberately trap semantic expansion of an unused shared-string entry.
        worksheet(f,mini,shared_probe=True)
        original=run.ET.fromstring
        class SI:
            def __init__(self,element,index):self.element,self.index=element,index
            def itertext(self):
                if self.index==0:raise AssertionError('Unselected shared string materialized')
                return self.element.itertext()
        class SST:
            def __init__(self,root):self.root=root
            def findall(self,*args):return [SI(e,i) for i,e in enumerate(self.root.findall(*args))]
        def controlled(data,*args,**kwargs):
            root=original(data,*args,**kwargs)
            return SST(root) if root.tag=='{'+NS+'}sst' else root
        run.ET.fromstring=controlled
        try:assert len(run.read_rows(f,mini['states'],mini['cell_contexts']))==12
        finally:run.ET.fromstring=original
        tested.append('unselected_shared_string_header_body_formula_and_row_never_decoded')
        for value in (F(0),F(-1)):
            worksheet(f,mini,{2:{'H':F(1),'I':F(1),'J':value}})
            rs=run.read_rows(f,mini['states'],mini['cell_contexts'])
            assert rs[0]['ratio']==1 and rs[0]['score_difference']==value-1 and rs[0]['source_score_consistent'] is False
            assert run.analyze(rs,mini)['status']=='SOURCE_REPRODUCTION_FAILED'
        tested.append('zero_and_negative_finite_score_gross_disagreement_fail_source_reproduction')
        for q,score,okay in [(F(1),F('1.00000001'),True),(F(1),F('1.0000000100000001'),False),(F(100),F('100.000001'),True),(F(100),F('100.00000100000001'),False)]:
            worksheet(f,mini,{2:{'H':F(1),'I':q,'J':score}})
            rs=run.read_rows(f,mini['states'],mini['cell_contexts'])
            assert rs[0]['source_score_consistent'] is okay
        tested.append('file_consistency_exact_tolerance_boundary_no_scientific_epsilon')
        for changes in ({'I':None},{'I':F(0)},{'I':'Inf'},{'J':'NaN'}):
            worksheet(f,mini,{2:changes})
            rs=run.read_rows(f,mini['states'],mini['cell_contexts'])
            rr=run.analyze(rs,mini)
            assert len(rr['state_cells'])==4 and len(rr['source_rows'])==12
            bad=next(o for o in rr['state_cells'] if o['state']=='a' and o['cell']=='A')
            assert not bad['complete'] and (bad['lower'],bad['upper'])==(F(0),None)
            assert rr['sensitivities'][0]['group_cell_results'][0]['maximum_fold'] is None
        tested.append('blank_zero_nonfinite_measurement_unknown_state_retained')
        for kwargs in ({'omit':{2}},{'extra':True}):
            worksheet(f,mini,**kwargs)
            try:run.analyze(run.read_rows(f,mini['states'],mini['cell_contexts']),mini)
            except ValueError:pass
            else:raise AssertionError('Metadata row count violation accepted')
        tested.append('missing_or_extra_metadata_row_integrity_failure')
        worksheet(f,mini,{2:{'H':F(1),'I':F(4),'J':F(4)},3:{'H':F(2),'I':F(2),'J':F(1)},4:{'H':F(4),'I':F(1),'J':F(1,4)}})
        rs=run.read_rows(f,mini['states'],mini['cell_contexts'])
        assert [r['ratio'] for r in rs[:3]]==[F(4),F(1),F(1,4)]
        tested.append('within_row_RNA_DNA_pairing_preserved')
        # Isolated CLI layout uses invented values and no source workbook.
        layout=tmp/'layout';pilot=layout/'cprime_pilot';gate=layout/'cprime_gate'
        pilot.mkdir(parents=True);gate.mkdir()
        for name in ['run_pilot.py','interval_core.py','protocol.md','prepare_metadata.py','test_synthetic.py','synthetic_receipt.json']:
            (pilot/name).write_bytes((PILOT/name).read_bytes())
        worksheet(gate/'supp_Supplementary Data 8.xlsx',eligibility)
        for name in ['selected_assay_metadata.json','search_sequence_metadata.json','measurement_contract.json']:
            (gate/name).write_text('{"synthetic_only":true}\n')
        fake=copy.deepcopy(eligibility)
        fake['input_sha256']={name:sha(gate/name) for name in eligibility['input_sha256']}
        (pilot/'eligibility.json').write_text(json.dumps(fake))
        pins={name:sha(pilot/name) for name in ['run_pilot.py','interval_core.py','protocol.md','prepare_metadata.py','eligibility.json','test_synthetic.py','synthetic_receipt.json']}
        (pilot/'freeze.json').write_text(json.dumps({'files':pins}))
        auth=layout/'auth.json'
        auth.write_text(json.dumps({'permit_outcomes':True,'freeze_sha256':sha(pilot/'freeze.json')}))
        env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'}
        def cli(output,flags=(),auth_path=auth):
            return subprocess.run([sys.executable,*flags,str(pilot/'run_pilot.py'),'--authorization',str(auth_path),'--out',str(output)],capture_output=True,text=True,env=env)
        good=layout/'success'
        completed=cli(good)
        assert completed.returncode==0,completed.stderr
        assert json.loads((good/'receipt.json').read_text())['source_rows']==396
        guards.append({'case':'authorized_synthetic_only_CLI','passed':True})
        for flags in [('-O',),('-OO',)]:
            target=layout/('optimized'+''.join(flags))
            r=cli(target,flags)
            assert r.returncode!=0 and 'Optimized execution is not permitted' in r.stderr and not target.exists()
            guards.append({'case':''.join(flags),'rejected_before_output':True})
        r=cli(good);assert r.returncode!=0 and 'FileExistsError' in r.stderr
        guards.append({'case':'existing_output','rejected':True})
        denied=layout/'denied.json';denied.write_text('{"permit_outcomes":false}')
        target=layout/'denied_output';r=cli(target,auth_path=denied)
        assert r.returncode!=0 and 'Missing matching outcome authorization' in r.stderr and not target.exists()
        guards.append({'case':'authorization_denied','rejected_before_output':True})
        altered=layout/'wrong_freeze.json';altered.write_text('{"permit_outcomes":true,"freeze_sha256":"wrong"}')
        target=layout/'wrong_freeze_output';r=cli(target,auth_path=altered)
        assert r.returncode!=0 and 'Missing matching outcome authorization' in r.stderr and not target.exists()
        guards.append({'case':'wrong_freeze_hash','rejected_before_output':True})
        core_path=pilot/'interval_core.py';original_core=core_path.read_bytes();canary=layout/'core_imported'
        core_path.write_bytes(('from pathlib import Path\nPath('+repr(str(canary))+').write_text("executed")\n').encode()+original_core)
        target=layout/'altered_core_output';r=cli(target)
        assert r.returncode!=0 and 'Frozen input changed: interval_core.py' in r.stderr and not canary.exists() and not target.exists()
        core_path.write_bytes(original_core)
        guards.append({'case':'modified_scientific_module','rejected_before_import_canary_and_output':True})
        workbook=gate/'supp_Supplementary Data 8.xlsx';original_workbook=workbook.read_bytes();workbook.write_bytes(original_workbook+b'toy-change')
        target=layout/'altered_source_output';r=cli(target)
        assert r.returncode!=0 and 'Source input changed' in r.stderr and not target.exists()
        guards.append({'case':'changed_synthetic_source','rejected_before_output':True})
    assert sha(PILOT/'run_pilot.py')==RUNNER_PIN and sha(PILOT/'interval_core.py')==CORE_PIN
    report={'status':'pass','source_runner_sha256':RUNNER_PIN,'core_sha256':CORE_PIN,'synthetic_checks':tested,'guard_checks':guards,
            'actual_source_workbook_opened':False,'production_analysis_run':False,'only_cli_run':'Isolated copied layout with entirely invented 396-row workbook',
            'shared_source_files_modified':False}
    (HERE/'reader_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
