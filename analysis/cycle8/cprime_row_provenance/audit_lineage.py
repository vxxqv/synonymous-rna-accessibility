"""Source lineage and identity-only row-to-candidate-library audit.
No phenotype, RPM, count, notebook-output or missingness value is accessed.
"""
from pathlib import Path
import ast,csv,hashlib,json,re,zipfile,xml.etree.ElementTree as ET
from collections import Counter
from trace_notebook_source import code_sources

P=Path(__file__).resolve().parent;G=P.parent/'cprime_gate';E=P.parent/'cprime_pilot/eligibility.json'
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')

def permitted_identity_cells(path):
    """Decode only A/B/E/F/G cell payloads; all others remain raw bytes."""
    allowed={'A':'combined_var_id','B':'var_id','E':'exon_number','F':'syn_var','G':'cell_line'}
    with zipfile.ZipFile(path) as z:
        book=ET.fromstring(z.read('xl/workbook.xml'))
        sheet=next(s for s in book.find('s:sheets',NS) if s.attrib['name']=='Supplementary_Table_3_Multi_exo')
        rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        target=next(r.attrib['Target'] for r in rels if r.attrib['Id']==rid)
        member=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        strings=z.read('xl/sharedStrings.xml') if 'xl/sharedStrings.xml' in z.namelist() else b''
        string_spans=list(re.finditer(rb'<si(?:\s[^>]*)?>.*?</si>',strings,re.S))
        def decode(cell):
            el=ET.fromstring(cell);kind=el.attrib.get('t')
            value=next((a.text for a in el if a.tag.rsplit('}',1)[-1]=='v'),None)
            if kind=='s':
                item=ET.fromstring(string_spans[int(value)].group())
                return ''.join(e.text or '' for e in item.iter() if e.tag.rsplit('}',1)[-1]=='t')
            if kind=='inlineStr':return ''.join(e.text or '' for e in el.iter() if e.tag.rsplit('}',1)[-1]=='t')
            return value
        rows={}
        # Inspect cell reference first, then decode only the five allowed columns.
        for match in re.finditer(rb'<c\s[^>]*?(?:/>|>.*?</c>)',z.read(member),re.S):
            ref=re.search(rb'\br="([A-Z]+)(\d+)"',match.group())
            if not ref:raise ValueError('Cell reference missing')
            col=ref.group(1).decode('ascii');row=int(ref.group(2))
            if col not in allowed or row==1:continue
            rows.setdefault(row,{})[allowed[col]]=decode(match.group())
        return rows,{'sheet':sheet.attrib['name'],'xml_member':member,'allowed_columns':allowed,
                     'external_link_members':[n for n in z.namelist() if n.startswith('xl/externalLinks/')],
                     'defined_names':[n.attrib for n in book.findall('s:definedNames/s:definedName',NS)]}

def main():
    if (P/'lineage_receipt.json').exists():raise FileExistsError('Fresh receipt required')
    e=json.loads(E.read_text(encoding='utf8'));states=e['states'];tree=json.loads((G/'github_tree.json').read_text(encoding='utf8'))
    files={r['path']:r for r in tree['tree'] if r['type']=='blob'}
    checks={}
    paths={'multi_exon_screen/analysis/analyze.ipynb':'author_analyze.ipynb',
           'multi_exon_screen/count_edits_in_fastq/count_ss.py':'author_count_ss.py',
           'multi_exon_screen/count_edits_in_fastq/count_edits_fast.sh':'author_count_edits_fast.sh',
           'multi_exon_screen/count_edits_in_fastq/merge_output.py':'author_merge_output.py'}
    for author,local in paths.items():
        b=(G/local).read_bytes();blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        if blob!=files[author]['sha']:raise ValueError('Pinned cached source mismatch')
        checks[author]={'git_blob_sha1':blob,'sha256':sha(G/local)}
    cells=dict(code_sources((G/'author_analyze.ipynb').read_bytes()))
    export_cell=next((i,s) for i,s in cells.items() if '/rpm_new.tsv' in s)
    assignments=[a for a in ast.parse(export_cell[1]).body if isinstance(a,ast.Assign)]
    columns=next(ast.literal_eval(a.value.slice) for a in assignments if isinstance(a.value,ast.Subscript))
    schema=json.loads((G/'supplement_schema_headers_only.json').read_text(encoding='utf8'))
    headers=next(s['row1_text_only'] for s in schema if s['file']=='supp_Supplementary Data 8.xlsx' and s['sheet']=='Supplementary_Table_3_Multi_exo')
    if columns!=[x['text'] for x in headers]:raise ValueError('Source export/header mismatch')
    identities,workbook=permitted_identity_cells(G/'supp_Supplementary Data 8.xlsx')
    selected={r:x for r,x in identities.items() if x.get('combined_var_id') in states and x.get('cell_line') in e['cell_contexts']}
    if len(selected)!=396:raise ValueError('Wrong selected-row identity count')
    groupcounts=Counter((x['combined_var_id'],x['cell_line']) for x in selected.values())
    if len(groupcounts)!=132 or set(groupcounts.values())!={3}:raise ValueError('Wrong metadata triplets')
    with (G/'ena_samples.tsv').open(encoding='utf8',newline='') as f:ena=list(csv.DictReader(f,delimiter='\t'))
    libraries={}
    for r in ena:
        m=re.fullmatch(r'(DNA|RNA)_E(17|25|27|34|39|42)_(K562|KBM7)_R([123])',r['library_name'])
        if m:
            key=(m[1],int(m[2]),m[3],int(m[4]))
            if key in libraries:raise ValueError('Ambiguous archived library label')
            libraries[key]=r
    if len(libraries)!=72:raise ValueError('Wrong full multi-exon metadata inventory')
    pairs=[];pair_lookup={}
    for exon in sorted(e['exons']):
        for cell in e['cell_contexts']:
            for replicate in (1,2,3):
                key=f'E{exon}/{cell}/R{replicate}'
                dna=libraries[('DNA',exon,cell,replicate)];rna=libraries[('RNA',exon,cell,replicate)]
                record={'source_label_pair':key,'exon':exon,'cell_line':cell,'source_replicate_label':replicate,
                        **{f'{mod}_{field}':r[field] for mod,r in [('DNA',dna),('RNA',rna)] for field in ('library_name','run_accession','sample_accession')},
                        'pairing_basis':'same filename replicate/exon/cell key as author pivot; not proof of independent culture identity',
                        'total_reads_DNA':None,'total_reads_RNA':None,'RPM_total_status':'processed FLASH-count denominator not deposited in inspected metadata'}
                pairs.append(record);pair_lookup[key]=record
    row_candidates=[]
    for row,x in sorted(selected.items()):
        exon=int(float(x['exon_number']))
        if exon!=int(states[x['combined_var_id']]['exon_number']):raise ValueError('Eligibility exon mismatch')
        row_candidates.append({'worksheet_row':row,**x,'authoritative_replicate_assignment':None,
          'candidate_source_label_pairs':[f'E{exon}/{x["cell_line"]}/R{i}' for i in (1,2,3)],
          'status':'unresolved_three_way_row_assignment; no ordinal pairing assumed'})
    dump(P/'labelled_library_pairs.json',pairs);dump(P/'selected_row_candidate_map.json',row_candidates)
    receipt={'status':'PARTIAL_LINEAGE_RECOVERED_ROW_PAIRING_AND_RPM_TOTALS_UNRESOLVED',
      'author_commit':'30125415304664c33370cd7801f820a345834e40','source_blob_checks':checks,
      'selected_rows':len(selected),'selected_state_cell_triplets':len(groupcounts),'selected_states':len(states),'selected_intended_edits':e['groups_count'],
      'selected_exons':e['exons'],'source_labelled_RNA_DNA_pairs':len(pairs),'source_labelled_libraries':len(pairs)*2,
      'authoritative_row_replicate_assignments':0,'authenticated_RPM_total_values':0,
      'export_cell_zero_based':export_cell[0],'export_schema_columns':len(columns),'export_schema_exactly_matches_supplement':True,
      'omitted_fields':['replicate','file','counts_complete','total_reads'],
      'missing_in_pinned_repository':[n for n in ('wide_annotated_new.tsv','rpm_new.tsv','Supplementary_Table_3_Multi_exo','Supplementary Data 8.xlsx') if not any(n in s for s in files)],
      'workbook_metadata':workbook,'new_phenotype_cells_accessed':0,'notebook_outputs_decoded':False,'source_code_executed':False,
      'input_sha256':{str(x.relative_to(P.parent)):sha(x) for x in [E,G/'github_tree.json',G/'ena_samples.tsv',G/'supp_Supplementary Data 8.xlsx',G/'supplement_schema_headers_only.json']},
      'output_sha256':{n:sha(P/n) for n in ('selected_row_candidate_map.json','labelled_library_pairs.json')},'audit_code_sha256':sha(__file__)}
    dump(P/'lineage_receipt.json',receipt);print(json.dumps({k:v for k,v in receipt.items() if k not in ('source_blob_checks','input_sha256','output_sha256','workbook_metadata')},indent=2))

if __name__=='__main__':main()
