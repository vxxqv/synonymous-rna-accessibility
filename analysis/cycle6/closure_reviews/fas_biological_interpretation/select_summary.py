"""Inspect authorized top-level summary/certificate only, skip outcome arrays."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'fas_results'/'full_result_v1.json'
WHITELIST = {'status','scope','vertices','edges','finite_envelopes','unknown_envelopes',
             'codon_design_columns','nucleotide_design_columns','primary_codon_null',
             'conservative_barrier_pairs','freeze_sha256','limitations'}
raw=SOURCE.read_bytes();text=raw.decode('utf8');decoder=json.JSONDecoder()
def whitespace(i):
    while i<len(text) and text[i].isspace():i+=1
    return i
def skip(i):
    start=i;i=whitespace(i)
    if text[i] in '[{':
        depth=0;quote=False;escape=False
        for j in range(i,len(text)):
            c=text[j]
            if quote:
                if escape:escape=False
                elif c=='\\':escape=True
                elif c=='"':quote=False
            elif c=='"':quote=True
            elif c in '[{':depth+=1
            elif c in ']}':
                depth-=1
                if depth==0:return j+1
        raise ValueError('Unclosed JSON structure')
    if text[i]=='"':
        escape=False
        for j in range(i+1,len(text)):
            c=text[j]
            if escape:escape=False
            elif c=='\\':escape=True
            elif c=='"':return j+1
        raise ValueError('Unclosed JSON string')
    j=i
    while j<len(text) and text[j] not in ',}':j+=1
    return j
i=whitespace(0)
if text[i]!='{':raise ValueError('Top-level object required')
i+=1;selected={};skipped=[]
while True:
    i=whitespace(i)
    if text[i]=='}':break
    key,end=decoder.raw_decode(text,i);i=whitespace(end)
    if text[i]!=':':raise ValueError('Missing colon')
    i=whitespace(i+1);end=skip(i)
    if key in WHITELIST:selected[key]=json.loads(text[i:end])
    else:skipped.append(key)
    i=whitespace(end)
    if text[i]==',':i+=1
    elif text[i]!='}':raise ValueError('Unexpected separator')
(HERE/'authorized_summary.json').write_text(json.dumps(selected,indent=2)+'\n',encoding='utf8')
receipt={'source':str(SOURCE),'source_sha256':hashlib.sha256(raw).hexdigest(),
         'selected_keys':sorted(selected),'skipped_top_level_keys':skipped,
         'read_boundary':'Unselected values lexically skipped without JSON decoding or display. No count/ratio/bottleneck rows inspected.',
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'summary_exposure_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(selected,indent=2))
