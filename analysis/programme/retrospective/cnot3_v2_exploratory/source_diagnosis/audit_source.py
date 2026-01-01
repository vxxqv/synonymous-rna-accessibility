"""Independent bounded source diagnosis. No alternate cohort, grouping or models."""
import sys
if not __debug__: raise RuntimeError('Optimization is refused before input access')
sys.dont_write_bytecode=True
from collections import Counter
import argparse
import csv
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
import hashlib
import json
import math
from pathlib import Path
import posixpath
import re
import xml.etree.ElementTree as ET
import xml.sax
from xml.sax.handler import ContentHandler,feature_external_ges,feature_external_pes
import zipfile

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
WORKBOOK=PARENT/'data_eligibility/Carneiro_gene_level_TableS2.xlsx'
EXECUTION=PARENT/'execution/reproduction_v1'
PINNED_WORKBOOK='a338e7211fcb0ad6bed0b1bfdb57433a5a43266f3a8d93e473b1bf50c583d712'
PINNED_EXPOSURE='b8bae6cc9e71f7740f54a304b6182c108fa2e69d6fc577660e2522ec86490159'
PINNED_COHORT='cfd4cededb5be6719795fa88262793ca580eb267649127cdf9cf300f75d9df1e'
FIELDS={'A':'GENE_ID','B':'HL_WT','C':'HL_CNOT3KO','D':'HL_FC_LOG2','I':'LOG2_FC_IP','J':'SE_IP','K':'STAT_IP'}
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))

class Cells(ContentHandler):
    def __init__(self,rows):
        super().__init__(); self.allowed_rows=rows; self.cells={}; self.selected=None; self.active=None; self.parts=[]; self.tags=[]
    def startElement(self,name,attrs):
        self.tags.append(name)
        if name=='c':
            coordinate=attrs['r']; match=re.fullmatch(r'([A-Z]+)(\d+)',coordinate)
            if not match: raise ValueError('Malformed cell coordinate')
            col,row=match[1],int(match[2])
            self.selected=(row,col) if row in self.allowed_rows and col in FIELDS else None
            if self.selected is not None:
                if self.selected in self.cells: raise ValueError('Duplicate selected cell')
                self.cells[self.selected]={'type':attrs.get('t','n'),'style':int(attrs.get('s','0')),'value_raw':None,'coordinate':coordinate}
        if self.selected is not None and name=='f': raise ValueError('Selected formula refused')
        if self.selected is not None and name in ('v','t'):
            self.active=name; self.parts=[]
    def characters(self,content):
        # Never buffer or semantically decode unselected cells.
        if self.selected is not None and self.active is not None: self.parts.append(content)
    def endElement(self,name):
        if self.selected is not None and self.active==name:
            value=''.join(self.parts)
            if name=='t' and self.cells[self.selected]['type']=='inlineStr':
                self.cells[self.selected]['value_raw']=(self.cells[self.selected]['value_raw'] or '')+value
            elif name=='v': self.cells[self.selected]['value_raw']=value
            self.active=None; self.parts=[]
        if name=='c': self.selected=None; self.active=None
        self.tags.pop()

class Strings(ContentHandler):
    def __init__(self,needed): super().__init__(); self.needed=needed; self.index=-1; self.active=False; self.in_text=False; self.parts=[]; self.values={}
    def startElement(self,name,attrs):
        if name=='si': self.index+=1; self.active=self.index in self.needed; self.parts=[]
        if name=='t': self.in_text=True
    def characters(self,content):
        if self.active and self.in_text: self.parts.append(content)
    def endElement(self,name):
        if name=='t': self.in_text=False
        if name=='si':
            if self.active: self.values[self.index]=''.join(self.parts)
            self.active=False

def sax_parse(stream,handler):
    parser=xml.sax.make_parser()
    for feature in (feature_external_ges,feature_external_pes):
        try: parser.setFeature(feature,False)
        except xml.sax.SAXNotSupportedException: pass
    parser.setContentHandler(handler); parser.parse(stream)

def independent_cells(rows):
    with zipfile.ZipFile(WORKBOOK) as archive:
        workbook=ET.fromstring(archive.read('xl/workbook.xml'))
        relationships=ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        selected=[node for node in workbook.iter('{'+N+'}sheet') if node.attrib['name']=='TableS2']
        if len(selected)!=1: raise ValueError('Expected one data worksheet')
        rid=selected[0].attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        mapping={node.attrib['Id']:node.attrib['Target'] for node in relationships}
        target=mapping[rid]
        member=target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join('xl',target))
        if not member.startswith('xl/') or '..' in member.split('/'): raise ValueError('Unexpected worksheet target')
        collector=Cells(rows)
        with archive.open(member) as stream: sax_parse(stream,collector)
        needed={int(cell['value_raw']) for cell in collector.cells.values() if cell['type']=='s' and cell['value_raw'] is not None}
        strings=Strings(needed)
        with archive.open('xl/sharedStrings.xml') as stream: sax_parse(stream,strings)
        if needed!=set(strings.values): raise ValueError('Selected shared string missing')
        for cell in collector.cells.values():
            raw=cell['value_raw']
            cell['decoded']=strings.values[int(raw)] if cell['type']=='s' and raw is not None else raw
        styles=ET.fromstring(archive.read('xl/styles.xml'))
        formats={0:'General',1:'0',2:'0.00',9:'0%',10:'0.00%',11:'0.00E+00'}
        formats.update({int(x.attrib['numFmtId']):x.attrib['formatCode'] for x in styles.iter('{'+N+'}numFmt')})
        xfs=styles.find('{'+N+'}cellXfs')
        for cell in collector.cells.values():
            fid=int(xfs[cell['style']].attrib.get('numFmtId','0'))
            cell['number_format_id']=fid; cell['number_format']=formats.get(fid,'uninterpreted_builtin_'+str(fid))
        return collector.cells

def number(cell):
    if cell is None or cell['type'] in ('e','b'): return None
    try: value=Decimal(cell['decoded'])
    except (TypeError,ValueError,ArithmeticError): return None
    return value if value.is_finite() else None

def significand10(value):
    return value.quantize(Decimal(1).scaleb(value.adjusted()-9),rounding=ROUND_HALF_EVEN)

def logratio(ko,wt): return (ko/wt).ln()/Decimal(2).ln()

def lexical_quantum(cell):
    value=number(cell)
    return Decimal(1).scaleb(value.as_tuple().exponent) if value is not None else None

def display_quantum(cell):
    match=re.fullmatch(r'0(?:\.(0+))?',cell['number_format'])
    return Decimal(1).scaleb(-len(match[1] or '')) if match else None

def interval_compatible(wt,ko,stored,qwt,qko,qd):
    if qwt is None or qko is None or qd is None: return None
    if wt-qwt/2<=0 or ko-qko/2<=0: return None
    lower=logratio(ko-qko/2,wt+qwt/2); upper=logratio(ko+qko/2,wt-qwt/2)
    return {'compatible':lower<=stored+qd/2 and upper>=stored-qd/2,'log2_ratio_lower':str(lower),'log2_ratio_upper':str(upper),
            'stored_log2_lower':str(stored-qd/2),'stored_log2_upper':str(stored+qd/2)}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out',required=True)
    args=parser.parse_args(); out=Path(args.out).resolve()
    if not out.is_relative_to(HERE): raise ValueError('This bounded audit output must stay in its owned directory')
    if out.exists(): raise FileExistsError(out)
    inputs={'workbook':sha(WORKBOOK),'exposure':sha(EXECUTION/'exposure.json'),'cohort':sha(EXECUTION/'cohort_rows.json'),'result':sha(EXECUTION/'result.json')}
    if inputs['workbook']!=PINNED_WORKBOOK or inputs['exposure']!=PINNED_EXPOSURE or inputs['cohort']!=PINNED_COHORT: raise ValueError('Authorized inputs changed')
    exposure=load(EXECUTION/'exposure.json')
    second=[stage for stage in exposure['stages'] if stage['columns']==['B','C','D']]
    if len(second)!=1 or not isinstance(second[0]['row_selector'],list): raise ValueError('No unique frozen candidate-row selector')
    candidate_rows=set(second[0]['row_selector'])-{1}
    cells=independent_cells(candidate_rows|{1})
    for col,name in FIELDS.items():
        if cells[(1,col)]['decoded']!=name: raise ValueError('Source header changed')
    paired=[]; excluded=[]; mismatches=[]; errors=[]; round10_errors=[]; valid_ids=[]; stat_mismatches=[]; raw_compatible=0; display_compatible=0; number_formats=Counter(); raw_decimal_places=Counter(); independent_computed={}
    with localcontext() as context:
        context.prec=80
        for excel_row in sorted(candidate_rows):
            record={col:cells.get((excel_row,col)) for col in FIELDS}
            gene=record['A']['decoded'] if record['A'] else None
            if not isinstance(gene,str) or not re.fullmatch(r'ENSG\d{11}',gene) or gene in valid_ids: raise ValueError('Invalid/duplicate candidate ID')
            values={col:number(record[col]) for col in FIELDS if col!='A'}
            if any(values[x] is None for x in ('I','J','K')): raise ValueError('Frozen candidate is not finite IP')
            valid_ids.append(gene)
            se=values['J']; stat_expected=values['I']/se if se!=0 else None
            stat_tol=Decimal('1e-10')+Decimal('1e-8')*max(abs(values['K']),abs(stat_expected)) if stat_expected is not None else None
            if stat_expected is None or se<0 or abs(values['K']-stat_expected)>stat_tol: stat_mismatches.append(gene)
            reasons=[]
            for col in ('B','C','D'):
                if values[col] is None: reasons.append('missing_or_nonfinite_'+FIELDS[col])
            for col in ('B','C'):
                if values[col] is not None and values[col]<=0: reasons.append('nonpositive_'+FIELDS[col])
            if reasons:
                excluded.append({'Excel_row':excel_row,'GeneID':gene,'reasons':reasons,'HL_WT_cell_type':record['B']['type'] if record['B'] else 'absent','HL_KO_cell_type':record['C']['type'] if record['C'] else 'absent','stored_ratio_cell_type':record['D']['type'] if record['D'] else 'absent'})
                continue
            wt,ko,stored=(values[x] for x in ('B','C','D')); computed=logratio(ko,wt); residual=stored-computed
            independent_computed[gene]=computed
            tolerance=Decimal('1e-10')+Decimal('1e-8')*max(abs(stored),abs(computed))
            rounded_wt,rounded_ko=significand10(wt),significand10(ko)
            rounded_computed=logratio(rounded_ko,rounded_wt); rounded_error=stored-rounded_computed
            errors.append(abs(residual)); round10_errors.append(abs(rounded_error))
            paired.append({'Excel_row':excel_row,'GeneID':gene,'HL_WT':float(wt),'HL_CNOT3KO':float(ko),'LOG2_FC_IP':float(values['I']),'SE_IP':float(se)})
            for col in ('B','C','D'):
                number_formats[(col,record[col]['number_format'])]+=1
                raw_decimal_places[(col,-values[col].as_tuple().exponent)]+=1
            if abs(residual)>tolerance:
                lexical=interval_compatible(wt,ko,stored,lexical_quantum(record['B']),lexical_quantum(record['C']),lexical_quantum(record['D']))
                display=interval_compatible(wt,ko,stored,display_quantum(record['B']),display_quantum(record['C']),display_quantum(record['D']))
                raw_compatible+=bool(lexical and lexical['compatible']); display_compatible+=bool(display and display['compatible'])
                entry={'Excel_row':excel_row,'GeneID':gene,'HL_WT_XML':record['B']['decoded'],'HL_KO_XML':record['C']['decoded'],'HL_FC_LOG2_XML':record['D']['decoded'],
                       'computed_log2_KO_over_WT_80digit':str(computed),'stored_minus_computed':str(residual),'absolute_error':str(abs(residual)),
                       'frozen_tolerance':str(tolerance),'absolute_error_over_tolerance':str(abs(residual)/tolerance),
                       'HL_WT_numFmt':record['B']['number_format'],'HL_KO_numFmt':record['C']['number_format'],'ratio_numFmt':record['D']['number_format'],
                       'raw_lexical_rounding_compatible':None if lexical is None else lexical['compatible'],
                       'raw_lexical_interval':json.dumps(lexical),'display_rounding_compatible':None if display is None else display['compatible'],
                       'display_rounding_interval':json.dumps(display),'HL_WT_round10':str(rounded_wt),'HL_KO_round10':str(rounded_ko),
                       'computed_after10sig':str(rounded_computed),'stored_minus_after10sig':str(rounded_error)}
                mismatches.append(entry)
    # Only now compare independently reconstructed quantities to the frozen author outputs.
    author=load(EXECUTION/'result.json'); exported=load(EXECUTION/'cohort_rows.json')
    independent_ids={r['GeneID'] for r in paired}; author_ids=set(author['cohort_ids'])
    expected_failures={r['GENE_ID'] for r in author['half_life_ratio_diagnostics'] if not r['consistent']}
    if independent_ids!=author_ids or {r['GeneID'] for r in mismatches}!=expected_failures: raise AssertionError('Substantive cohort/mismatch disagreement')
    export_map={r['GeneID']:r for r in exported}
    for r in paired:
        if any(r[key]!=export_map[r['GeneID']][key] for key in ('HL_WT','HL_CNOT3KO','LOG2_FC_IP','SE_IP')): raise AssertionError('Exported selected values differ')
    decimal_vs_author=max(abs(Decimal(str(q['computed_log2_ratio']))-independent_computed[q['GENE_ID']]) for q in author['half_life_ratio_diagnostics'])
    out.mkdir()
    with (out/'ratio_mismatches.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(mismatches[0]) if mismatches else ['GeneID']); writer.writeheader(); writer.writerows(mismatches)
    with (out/'candidate_exclusions.json').open('x',encoding='utf-8') as stream: json.dump(excluded,stream,indent=2); stream.write('\n')
    summary={'status':'independent_source_diagnosis_complete','frozen_commit':'cc9cb18954026b9fe6199afe9820920823aaa767',
             'independent_parser':'SAX selected-cell collector; Decimal80digit ratio arithmetic; no author functions imported',
             'scope':'Frozen7359finite-IP candidate rows only, columnsA-D/I-K; no F/G or other endpoints decoded',
             'candidate_universe_provenance':'Frozen exposure stage2 row selector; no independent candidate search outside already exposed rows',
             'finite_IP_candidates_verified':len(candidate_rows),'paired_cohort_independent':len(paired),'published_cohort_n':5481,'cohort_ID_set_matches_frozen':True,
             'candidate_exclusion_count':len(excluded),'exclusion_reason_counts':dict(Counter(reason for r in excluded for reason in r['reasons'])),
             'cohort_ratio_failures':len(mismatches),'failed_ID_set_matches_frozen':True,'max_absolute_raw_ratio_error':str(max(errors)),
             'minimum_failed_absolute_error':str(min(Decimal(r['absolute_error']) for r in mismatches)),
             'maximum_failure_over_frozen_tolerance':str(max(Decimal(r['absolute_error_over_tolerance']) for r in mismatches)),
             'raw_last_decimal_rounding_compatible_failures':raw_compatible,'display_rounding_compatible_failures':display_compatible,
             'display_rounding_bounds_available_failures':sum(r['display_rounding_compatible'] is not None for r in mismatches),
             'maximum_native_author_vs_Decimal_computed_ratio_difference':str(decimal_vs_author),
             'observed_number_formats':[{'column':col,'format':fmt,'count':count} for (col,fmt),count in sorted(number_formats.items())],
             '10sig_diagnostic':'Explicit post-exposure source hypothesis supplied by root; not a new gate or replacement endpoint',
             '10sig_ROUND_HALF_EVEN_agreement_tolerance':'2e-13','10sig_all_cohort_rows_within_tolerance':sum(x<=Decimal('2e-13') for x in round10_errors),
             '10sig_max_absolute_error_80digit':str(max(round10_errors)),
             '10sig_generation_procedure_proven':False,'IP_statistic_failures_in_authorized_candidates':stat_mismatches,
             'published_missing_gene_identified':False,'exported_raw_numeric_values_match_frozen':True,
             'input_sha256_before':inputs,'input_sha256_after':{'workbook':sha(WORKBOOK),'exposure':sha(EXECUTION/'exposure.json'),'cohort':sha(EXECUTION/'cohort_rows.json'),'result':sha(EXECUTION/'result.json')},
             'decoded_cells':len(cells),'decoded_columns':list(FIELDS),'decoded_data_rows':len(candidate_rows),'protected_or_new_endpoint_columns_decoded':[],
             'code_sha256':sha(Path(__file__)),'scientific_protocol_or_failed_gate_modified':False,'new_models_or_grouping':False}
    if summary['input_sha256_before']!=summary['input_sha256_after']: raise AssertionError('Source inputs changed')
    with (out/'summary.json').open('x',encoding='utf-8') as stream: json.dump(summary,stream,indent=2); stream.write('\n')
    print(json.dumps({k:summary[k] for k in ('paired_cohort_independent','candidate_exclusion_count','exclusion_reason_counts','cohort_ratio_failures','max_absolute_raw_ratio_error','10sig_all_cohort_rows_within_tolerance','10sig_max_absolute_error_80digit','observed_number_formats')}))

if __name__=='__main__': main()
