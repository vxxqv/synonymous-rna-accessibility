"""Bounded public metadata retrieval. Never retrieves sequencing reads or outcomes."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request
import urllib.error

HERE = Path(__file__).resolve().parent
REQUESTS = [
 ('mordstein_cloud_inventory.xml','https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix=PMC7181179'),
 ('mittal_cloud_inventory.xml','https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix=PMC6112741'),
 ('ena_study.xml','https://www.ebi.ac.uk/ena/browser/api/xml/PRJNA596086'),
 ('SAMN13612106.xml','https://www.ebi.ac.uk/ena/browser/api/xml/SAMN13612106'),
 ('SAMN13612105.xml','https://www.ebi.ac.uk/ena/browser/api/xml/SAMN13612105'),
 ('ena_runs.json','https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA596086&result=read_run&fields=run_accession,experiment_accession,sample_accession,study_accession,library_name,experiment_title&format=json'),
 ('ena_analyses.json','https://www.ebi.ac.uk/ena/portal/api/search?result=analysis&query='+urllib.parse.quote('study_accession="PRJNA596086"')+'&fields=analysis_accession,study_accession&format=json'),
 ('mittal_supplementary_files.response','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6112741/supplementaryFiles'),
]

def fetch(item):
    name,url=item
    record={'url':url,'file':name,'requested_at_utc':datetime.now(timezone.utc).isoformat(),'limit_bytes':5_000_000}
    path=HERE/name
    if path.exists():raise FileExistsError(path)
    try:
        try:
            response=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'PublicResearchMetadataAudit/1.0'}),timeout=30)
        except urllib.error.HTTPError as error:
            response=error
        with response:
            body=response.read(5_000_001)
            record.update(status=response.status,final_url=response.url,content_type=response.headers.get('Content-Type'))
        if len(body)>5_000_000:
            record['stopped']='size_limit';return record
        path.write_bytes(body)
        record.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
    except Exception as error:
        record['error']=str(error)
    return record

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(fetch,REQUESTS))
    with (HERE/'retrieval_receipts.json').open('x',encoding='utf-8') as f:
        json.dump(results,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(results,indent=2))
