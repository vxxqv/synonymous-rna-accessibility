"""Authorized first-record source-format audit, never outcome/coverage analysis.

Read one physical first line per immutable count file. Numeric values are not
reported. Full-file SHA256/Git blob verification is byte provenance only.
"""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CONTRACT = ROOT/'fas_pilot'/'source_contract.json'
INPUTS = ROOT/'fas_inputs'
WT = 'GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
ANCESTOR = 'GATCCAGATCTAACTTGCTGTGGTTGTGTCTCCTGCTTCTCCCGATTCTAGTAATTGTTTGGG'
PATTERN = re.compile('GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG')
cfg = json.loads(CONTRACT.read_text(encoding='utf8'))
results = []
for item in cfg['files']:
    path = INPUTS/item['local_name']
    b=path.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    if blob != item['git_blob']:
        raise RuntimeError('Pinned source blob changed')
    with path.open('rb') as stream:
        first=stream.readline()
    fields=first.decode('utf8-sig').rstrip('\r\n').split('\t')
    is_header=fields==['Genotype','Occurrences','MutationCount']
    data=(len(fields)==3 and PATTERN.fullmatch(fields[0]) is not None
          and fields[1].isascii() and fields[1].isdigit()
          and fields[2].isascii() and fields[2].isdigit())
    results.append({'label':item['label'],'local_name':item['local_name'],
        'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'git_blob':blob,
        'expected_git_blob_matches':True,'first_line_sha256':hashlib.sha256(first).hexdigest(),
        'first_line_field_count':len(fields),'first_line_exact_header':is_header,
        'first_line_valid_declared_data_record':bool(data),
        'first_record_is_native_human_WT':bool(data and fields[0]==WT),
        'first_record_is_ancestor':bool(data and fields[0]==ANCESTOR),
        'first_record_mutation_annotation_matches_ancestor':bool(data and int(fields[2])==sum(a!=b for a,b in zip(fields[0],ANCESTOR))),
        'count_value':'REDACTED_NOT_DISPLAYED',
        'later_records_inspected':False})
receipt={'scope':'First physical line format/identity only, no phenotype analysis or coverage inspection',
         'contract_sha256':hashlib.sha256(CONTRACT.read_bytes()).hexdigest(),
         'files':results,'all_twelve_first_lines_are_data':len(results)==12 and all(r['first_line_valid_declared_data_record'] for r in results),
         'all_twelve_first_records_are_native_WT':len(results)==12 and all(r['first_record_is_native_human_WT'] for r in results),
         'numeric_counts_disclosed_by_this_script':False,
         'earlier_exposure':'Root reported BR1 and TR1 WT first-row counts; not newly used for selection here.',
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'first_record_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps({'all_twelve_first_lines_are_data':receipt['all_twelve_first_lines_are_data'],
                  'all_twelve_first_records_are_native_WT':receipt['all_twelve_first_records_are_native_WT'],
                  'labels':[r['label'] for r in results],
                  'counts_redacted':True,'later_records_inspected':False},indent=2))
