"""Narrow source-reference and selected promoter-metadata check, no outcomes."""
if not __debug__:raise RuntimeError('Optimization refused')
from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
WORK=HERE.parents[2]
PARENT=HERE.parent/'shanas_contract_pilot'
SOURCE=PARENT/'s4_followup'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,msg):
 if not x:raise RuntimeError(msg)
pins={'addendum.md':'3909a7a0e0fc0eb63199af11a006ff0ef226ecb5165f9cf6418ae19b8755c14d','closure_receipt.json':'43f7a0cf6bfc39b629cb90bc6424e32d90c92a18d0415c8b2dcc224ea8232496','delivery_manifest.json':'41a52d48501114beb4bf92f122626b90a44efb41eefb5702c56356d99789c6cd','s4_sequence_schema_summary.json':'266fab1f0cd73c0752c14fb251a2f21870f1157ca8697b4dbf621b376266d7ca'}
for k,v in pins.items():require(sha(SOURCE/k)==v,'Source pin changed '+k)
receipt=json.loads((SOURCE/'closure_receipt.json').read_bytes())
for k,v in receipt['frozen_parent_artifacts_unchanged'].items():require(sha(PARENT/k)==v,'Parent changed '+k)
require(sha(HERE/'review.md')=='fdd164d933c691dc542bf1510525428daad1d2b2ad08ad2c74ed967c6aa61340','Prior review changed')
require(sha(HERE/'receipt.json')=='a65187daabefd01eec12eae4ca90779e0f05b7629d611edaae438207c18d1964','Prior receipt changed')
primary=WORK/'work/joint_function_cycle5/successor_methods_review/shanas2026.xml'
require(sha(primary)=='9d05575f8002c5478ded9cba8677500d0e85476a2dc3d73435be6f9716444d20','Primary source changed')
root=ET.parse(primary).getroot()
media=[x for x in root.iter('media') if x.get('id')=='MOESM6']
require(len(media)==1,'Ambiguous supplementary reference')
require(''.join(media[0].itertext()).strip()=='Supplementary Data 4','Supplement label mismatch')
filename=media[0].get('{http://www.w3.org/1999/xlink}href')
require(filename=='41467_2026_76370_MOESM6_ESM.xlsx','Supplement filename mismatch')
inventory=json.loads((SOURCE/'official_inventory_receipt.json').read_bytes())
entry=[x for x in inventory['files'] if x['key']=='PMC13534650.1/'+filename]
require(len(entry)==1,'Official inventory mismatch')
workbook=SOURCE/'official_SupplementaryData4.xlsx'
require(workbook.stat().st_size==int(entry[0]['size'])==11971,'File size mismatch')
require(sha(workbook)==receipt['source_workbook_sha256']=='9acf21d4b939451248dbebb911fc3769228a079ea64da6404db49779eeac316e','Workbook byte mismatch')
require(hashlib.md5(workbook.read_bytes()).hexdigest()==entry[0]['etag'].strip('"')==receipt['source_workbook_md5'],'Single-object ETag mismatch')
selected=json.loads((SOURCE/'selected_s4_metadata.json').read_text())
summary=json.loads((SOURCE/'s4_sequence_schema_summary.json').read_text())
require(len(selected)==len(summary['records'])==11,'Row count mismatch')
require(len({x['promoter construct'] for x in selected})==11,'Duplicate construct labels')
index={x['construct']:x for x in summary['records']}
ambiguous={}
for row in selected:
 require(set(row)=={'promoter construct','sequence'},'Unexpected/protected metadata field')
 seq=row['sequence'].upper();label=row['promoter construct'];s=index[label]
 require(set(seq)<=set('ACGTN'),'Unexpected sequence character')
 require(len(seq)==s['sequence_length'],'Length disagreement')
 require(seq.count('N')==s['ambiguous_N_count'],'Ambiguity disagreement')
 if seq.count('N'):ambiguous[label]=seq.count('N')
require(ambiguous=={'pACTB':1},'Ambiguous-base summary mismatch')
require(min(len(x['sequence']) for x in selected)==46 and max(len(x['sequence']) for x in selected)==1703,'Length range mismatch')
require('\u2014' not in (HERE/'s4_addendum.md').read_text(encoding='utf-8'),'Punctuation constraint')
out={'verdict':'S4_ACCESS_RESOLVED_SAME_PROTEIN_GATE_UNCHANGED','blocking_findings':[],'primary_cross_reference_verified':True,'official_inventory_file_bytes_verified':True,'selected_metadata_rows':11,'selected_columns':['promoter construct','sequence'],'ambiguous_N':ambiguous,'annotated_full_CDS_supplied':False,'synonymous_graph_counts':'not authenticated or computed; not zero','source_pins':pins,'source_workbook_sha256':sha(workbook),'selected_metadata_sha256':sha(SOURCE/'selected_s4_metadata.json'),'parent_report_evidence_selector_unchanged':receipt['frozen_parent_artifacts_unchanged'],'prior_review_and_receipt_unchanged':True,'addendum_sha256':sha(HERE/'s4_addendum.md'),'audit_code_sha256':sha(Path(__file__)),'phenotype_or_missingness_access':False,'new_retrieval_or_shared_edit':False,'limits':'Metadata/reference audit only. No assertion that absent exact GFP anchor exhaustively rules out coding; no measured overlap or validation result computed.'}
target=HERE/'s4_addendum_receipt.json'
if target.exists():raise RuntimeError('No overwrite')
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'verdict':out['verdict'],'addendum_sha256':out['addendum_sha256'],'receipt_sha256':sha(target)},indent=2))
