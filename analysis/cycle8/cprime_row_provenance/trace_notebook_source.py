"""Pinned public notebook retrieval; only code-source JSON values are decoded.

Notebook outputs, execution counts and metadata values are lexically skipped.
No author code is executed. No measurement workbook cells are accessed.
"""
from pathlib import Path
import hashlib,json,re,urllib.request

P=Path(__file__).resolve().parent
G=P.parent/'cprime_gate'
COMMIT='30125415304664c33370cd7801f820a345834e40'

def ws(b,i):
    while i<len(b) and b[i] in b' \t\r\n':i+=1
    return i

def end(b,i):
    i=ws(b,i)
    if b[i]==34:
        j=i+1
        while j<len(b):
            if b[j]==92:j+=2
            elif b[j]==34:return j+1
            else:j+=1
        raise ValueError('Unterminated JSON string')
    if b[i] in b'[{':
        close=93 if b[i]==91 else 125;j=i+1
        while True:
            j=ws(b,j)
            if b[j]==close:return j+1
            if b[j] in b',:':j+=1
            else:j=end(b,j)
    j=i
    while j<len(b) and b[j] not in b',]} \r\n\t':j+=1
    return j

def members(b,i):
    i=ws(b,i)
    if b[i]!=123:raise ValueError('Object required')
    i+=1
    while True:
        i=ws(b,i)
        if b[i]==125:return
        k=end(b,i);key=json.loads(b[i:k]);i=ws(b,k)
        if b[i]!=58:raise ValueError('Object colon')
        start=ws(b,i+1);stop=end(b,start)
        yield key,start,stop
        i=ws(b,stop)
        if b[i]==44:i+=1
        elif b[i]!=125:raise ValueError('Object separator')

def code_sources(b):
    cells=next((a,z) for k,a,z in members(b,0) if k=='cells')
    i=ws(b,cells[0]);assert b[i]==91;i+=1;index=0
    while True:
        i=ws(b,i)
        if b[i]==93:return
        fields={k:(a,z) for k,a,z in members(b,i)}
        a,z=fields['cell_type'];kind=json.loads(b[a:z])
        if kind=='code':
            a,z=fields['source'];source=json.loads(b[a:z])
            yield index,''.join(source) if isinstance(source,list) else source
        i=ws(b,end(b,i));index+=1
        if b[i]==44:i+=1
        elif b[i]!=93:raise ValueError('Cell separator')

def main():
    output=P/'gmm_source_only.py';receipt=P/'gmm_retrieval.json'
    if output.exists() or receipt.exists():raise FileExistsError('Fresh source trace required')
    tree=json.loads((G/'github_tree.json').read_text(encoding='utf8'))['tree']
    path='multi_exon_screen/analysis/gmm.ipynb'
    expected=next(x for x in tree if x['path']==path)
    url=f'https://raw.githubusercontent.com/Boettcher-Lab/cPRIME/{COMMIT}/{path}'
    request=urllib.request.Request(url,headers={'User-Agent':'Source-provenance-audit/1.0'})
    with urllib.request.urlopen(request,timeout=35) as response:
        raw=response.read(1_000_001);status=response.status
    if len(raw)>1_000_000:raise ValueError('Bounded source size exceeded')
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if blob!=expected['sha'] or len(raw)!=expected['size']:raise ValueError('Pinned source blob mismatch')
    text='\n\n'.join(f'# AUTHOR NOTEBOOK CELL {i}\n{s}' for i,s in code_sources(raw))
    (P/'gmm_pinned_unread_outputs.ipynb').write_bytes(raw)
    output.write_text(text+'\n',encoding='utf8')
    record={'url':url,'status':status,'author_commit':COMMIT,'git_blob_sha1':blob,
            'raw_sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),
            'source_only_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
            'outputs_and_metadata_JSON_values_decoded':False,'author_code_executed':False,
            'selective_parser_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    receipt.write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps(record,indent=2))
    for i,line in enumerate(text.splitlines(),1):
        if re.search(r'Supplementary|to_excel|to_csv|read_csv|read_excel|replicate|sort_|drop\(|concat',line,re.I):
            print(i,line)

if __name__=='__main__':main()
