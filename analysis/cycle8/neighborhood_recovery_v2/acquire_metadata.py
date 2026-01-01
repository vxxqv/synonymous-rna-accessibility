"""Bounded public repository searches. Metadata only, no data downloads."""
from pathlib import Path
import hashlib,json,requests
P=Path(__file__).resolve().parent
cases=[('tp53_figshare_doi',{'resource_doi':'10.1158/1541-7786.MCR-17-0245','page_size':100}),('tp53_figshare_title',{'search_for':'A Massively Parallel Fluorescence Assay to Characterize the Effects of Synonymous Mutations on TP53 Expression','page_size':100}),('schmitz_public_datasets',{'search_for':'Schmitz','item_type':3,'page_size':100})]
ledger=[]
for name,payload in cases:
    out=P/(name+'.json')
    if out.exists():raise FileExistsError('No duplicate request '+name)
    url='https://api.figshare.com/v2/articles/search'
    try:
        r=requests.post(url,json=payload,timeout=25)
        out.write_bytes(r.content)
        item={'id':name,'method':'POST public read-only search','url':url,'payload':payload,'status':r.status_code,'final_url':r.url,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()}
        if r.status_code==200:
            data=r.json();item['records']=[{k:x.get(k) for k in ('id','title','url_public_html','url','resource_doi','resource_title','defined_type_name')} for x in data]
        else:item['error_body_file']=out.name
    except requests.RequestException as e:item={'id':name,'url':url,'payload':payload,'error':str(e),'no_retry':True}
    ledger.append(item)
(P/'repository_search_receipts.json').write_text(json.dumps(ledger,indent=2)+'\n')
for item in ledger:print(json.dumps(item,ensure_ascii=False))
