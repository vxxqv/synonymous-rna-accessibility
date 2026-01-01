"""Recover pinned public inputs without running external analysis code."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);a=ap.parse_args()
    root=a.work.resolve();root.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[2]
    manifests=[json.loads((repo/'data/programme'/name).read_text()) for name in ['measured_source_manifest.json','measured_endpoint_source_manifest.json']]
    for item in [x for manifest in manifests for x in manifest['files']]:
        dest=(root/item['work_relative_path']).resolve()
        if not dest.is_relative_to(root):raise ValueError('Input path escapes source root')
        if dest.exists():
            if hashlib.sha256(dest.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('Existing source differs: '+str(dest))
            print('Verified existing',item['work_relative_path']);continue
        with urlopen(Request(item['url'],headers={'User-Agent':'measured-function-reproducibility/1.0'}),timeout=60) as r:
            data=r.read(item['bytes']+1_000_000)
        variants=[data]
        if dest.suffix.lower() in ['.fa','.fasta']:
            normal=data.replace(b'\r\n',b'\n').rstrip(b'\n')
            variants.extend(normal+b'\n'*count for count in range(1,5))
        accepted=next((v for v in variants if hashlib.sha256(v).hexdigest()==item['sha256']),None)
        if accepted is None:raise ValueError('Pinned source hash mismatch; do not substitute a newer input: '+item['work_relative_path'])
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(accepted)
        print('Recovered and verified',item['work_relative_path'])
if __name__=='__main__':main()
