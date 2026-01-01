"""Bounded source schema audit. No phenotype rows or notebook outputs inspected."""
from pathlib import Path
import urllib.request,hashlib,gzip,json
HERE=Path(__file__).resolve().parent
tree=json.loads((HERE/'source_tree.json').read_text())
base='https://raw.githubusercontent.com/MustoeLab/Publications/'+tree['sha']+'/'
paths=['Dao_MPRAsplicing/splice_prediction/ptreseq_mpra_seq.fa.gz',
       'Dao_MPRAsplicing/data/barcode_indexed.txt',
       'Dao_MPRAsplicing/process_raw_count_to_splicing_efficiency.ipynb']
receipt=[]
for path in paths:
    rec=next(x for x in tree['tree'] if x['path']==path)
    req=urllib.request.Request(base+path,headers={'User-Agent':'Research-source-audit'})
    with urllib.request.urlopen(req,timeout=45) as f: raw=f.read()
    sha=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if sha!=rec['sha']:raise ValueError('Blob mismatch')
    dest=HERE/path.split('/')[-1];dest.write_bytes(raw)
    receipt.append({'path':path,'source_commit':tree['sha'],'source_blob':sha,'sha256':hashlib.sha256(raw).hexdigest()})
    if path.endswith('.ipynb'):
        doc=json.loads(raw)
        code='\n\n'.join(''.join(c['source']) for c in doc['cells'] if c['cell_type']=='code')
        (HERE/'processing_code_only.txt').write_text(code,encoding='utf-8')
    elif path.endswith('.fa.gz'):
        records={}
        for line in gzip.decompress(raw).decode().splitlines():
            if line.startswith('>'):key=line[1:];records[key]=''
            else:records[key]+=line.strip()
        (HERE/'fasta_metadata.json').write_text(json.dumps({'records':len(records),'lengths':sorted(set(map(len,records.values()))),'id_examples':list(records)[:6]},indent=2)+'\n')
(HERE/'acquisition_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'files':len(receipt),'fasta':json.loads((HERE/'fasta_metadata.json').read_text()),'phenotype_values_read':False}))
