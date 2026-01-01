"""Audit intended reference objects only; never decode any measured outcome."""
from pathlib import Path
from collections import Counter
import gzip, hashlib, json, re
from Bio.Seq import Seq

HERE=Path(__file__).resolve().parent
path=HERE/'ptreseq_mpra_seq.fa.gz'
receipt=json.loads((HERE/'acquisition_receipt.json').read_text())
expected=next(r for r in receipt if r['path'].endswith('/'+path.name))
assert hashlib.sha256(path.read_bytes()).hexdigest()==expected['sha256']
records={}; key=None
for line in gzip.decompress(path.read_bytes()).decode('ascii').splitlines():
    if line.startswith('>'):
        key=line[1:];assert key not in records;records[key]=''
    else:
        assert key is not None;records[key]+=line.strip()
assert len(records)==6500 and set(map(len,records.values()))=={1045}
reference=next(iter(records.values()))
assert reference.startswith('ATG')
stop=next(i+3 for i in range(0,len(reference)-2,3) if reference[i:i+3] in {'TAA','TAG','TGA'})
protein=str(Seq(reference[:stop]).translate())
assert protein.endswith('*') and '*' not in protein[:-1]
same=all(seq[:stop]==reference[:stop] for seq in records.values())
first_vary=next((i for i in range(len(reference)) if len({s[i] for s in records.values()})>1),None)
groups=Counter(re.sub(r'_[0-9]+$','',name) for name in records)
result={'source':'author full intended transcript FASTA; not observed isoforms',
 'records':len(records),'reference_orf_amino_acids':len(protein)-1,'reference_orf_end_exclusive':stop,
 'all_records_same_reference_ORF':same,'reference_variation_first_index':first_vary,
 'reporter_suffix_groups':len(groups),'barcode_group_multiplicity_counts':dict(Counter(groups.values())),
 'non10_groups':{k:v for k,v in groups.items() if v!=10},'phenotypes_inspected':False,'spliced_products_not_inferred':True}
prior=json.loads((HERE/'reference_object_audit.json').read_text())
assert json.loads(json.dumps(result))==prior, 'Independent reference audit differs'
(HERE/'reference_object_replay.json').write_text(json.dumps({'status':'PASS','source_sha256':expected['sha256'],'result':result},indent=2)+'\n')
print(json.dumps({'status':'PASS','records':len(records),'groups':len(groups),'reference_ORF_aa':len(protein)-1}))
