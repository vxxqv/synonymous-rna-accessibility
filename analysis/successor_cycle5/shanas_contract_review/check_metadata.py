"""Independent check of already selected seven-field metadata only."""
if not __debug__:raise RuntimeError('Optimization refused')
from pathlib import Path
from collections import Counter
import hashlib,json,ast
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'shanas_contract_pilot'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,msg):
 if not x:raise RuntimeError(msg)
allowed={'library index','subset','gene name','exonlength','exonstart in sequence','exonend in sequence','sequence'}
datasets=[]
for name in ['selected_library_metadata.json','selected_workbook_metadata.json']:
 rows=json.loads((SOURCE/name).read_text(encoding='utf-8'))
 require(len(rows)==9608,'Wrong metadata row count')
 for row in rows:
  require(set(row)<=allowed and allowed-{'gene name'}<=set(row),'Extra/protected or missing metadata field')
  require(set(row['sequence'].upper())<=set('ACGT'),'Non-DNA design sequence')
  require(len(row['sequence'])==210,'Unexpected construct fragment length')
 by_id={r['library index']:r for r in rows}
 require(len(by_id)==9608,'Duplicate metadata IDs')
 datasets.append(by_id)
require(set(datasets[0])==set(datasets[1]),'CSV/XLSX ID mismatch')
for id in datasets[0]:
 for key in allowed:
  require((datasets[0][id].get(key) or '')==(datasets[1][id].get(key) or ''),'Metadata mismatch '+id+' '+key)
rows=list(datasets[0].values())
native=[r for r in rows if r['subset']=='cassette_filtered']
counts={'rows':len(rows),'unique_library_IDs':len(datasets[0]),'unique_full210_sequences':len({r['sequence'].upper() for r in rows}),'unique_barcodes_18_30':len({r['sequence'][18:30].upper() for r in rows}),'unique_variable147_sequences_45_192':len({r['sequence'][45:192].upper() for r in rows}),'cassette_filtered_rows':len(native),'cassette_filtered_unique_genes':len({r['gene name'] for r in native}),'cassette_filtered_unique_variable_sequences':len({r['sequence'][45:192].upper() for r in native}),'cassette_filtered_gene_multiplicities':dict(Counter(Counter(r['gene name'] for r in native).values())),'unlabeled_gene_rows':sum(not r.get('gene name') for r in rows),'exon_length_mod3_counts':dict(Counter(int(r['exonlength'])%3 for r in rows))}
evidence=json.loads((SOURCE/'evidence.json').read_text(encoding='utf-8'))
for key,value in counts.items():
 require(json.loads(json.dumps(value))==evidence['sequence_metadata_counts'][key],'Count mismatch '+key)
require(evidence['synonymous_edge_count'] is None and evidence['synonymous_square_count'] is None,'Unknown graphs changed to number')
tree=ast.parse((SOURCE/'author_code/ce.pipeline.py').read_text(encoding='utf-8'))
constants={}
for node in ast.walk(tree):
 if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ['const0','const1','const2']:
  require(isinstance(node.value,ast.Constant),'Nonliteral reference component')
  constants[node.targets[0].id]=node.value.value.upper()
require([len(constants[k]) for k in ['const0','const1','const2']]==[717,288,833],'Unexpected source constants')
require(30+288+(210-45)+833==1316,'Amplicon arithmetic error')
require(constants['const1'][58:60]=='GT' and constants['const2'][729:731]=='AG','Canonical splice-site coordinate mismatch')
require(45+273==318 and 162+273==435 and 30+288+(210-45)+731==1214,'Coordinate arithmetic mismatch')
result={'verdict':'METADATA_IDENTITIES_COUNTS_AND_REFERENCE_ARITHMETIC_AGREE','counts':counts,'numeric_phenotype_fields_read':False,'assay_coverage_missingness_read':False,'author_code_executed':False,'synonymous_graph_counts':'not computed; not zero','reference':'1316nt mapping amplicon, not authenticated full coding construct','source_hashes':{x:sha(SOURCE/x) for x in ['selected_library_metadata.json','selected_workbook_metadata.json','evidence.json','selective_readers.py','author_code/ce.pipeline.py','author_code/notebook_code_only.py']},'audit_source_sha256':sha(Path(__file__))}
out=HERE/'metadata_checks.json'
if out.exists():raise RuntimeError('No overwrite')
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'verdict':result['verdict'],'receipt_sha256':sha(out)},indent=2))
