"""Creates and runs invented preprocessed RF fixtures, no source archive access."""
from pathlib import Path
import csv,json,hashlib,random,subprocess
import selective_reader as sr
from test_implementation import synthetic_row

P=Path(__file__).resolve().parent
OUT=P/'synthetic_RF_v1'


def main():
    if OUT.exists():raise FileExistsError('Fresh synthetic RF directory required')
    prepared=OUT/'prepared';prepared.mkdir(parents=True)
    rng=random.Random(1208)
    fields=list(sr.SELECTED)+['Sample','ROI','source_line','record_id']
    for roi in (1,2):
        rows=[];truth=[]
        for i in range(240):
            pos=i+1 if roi==1 else i+301
            r=synthetic_row(sample=f'COMT_{107+roi}_train',pos=pos,dp=1000,cao=20 if i%2 else 1,line=i+2)
            r['ROI']='ROI_'+str(roi)
            r['UP_REF_NT']='ACGT'[i%4];r['DOWN_REF_NT']='ACGT'[(i//4)%4]
            for n in sr.NUMERIC[1:]:r[n]=rng.uniform(5,40)
            rows.append(r)
            truth.append({'#CHROM':'toy','POS':pos,'REF':'A','ALT':'C','INFO':f'CAO={20 if i%2 else 0};CAF={.02 if i%2 else 0}'})
        sr.write_tsv(prepared/f'simulation_{roi}.tsv',rows,fields)
        sr.write_tsv(prepared/f'truth_{roi}.tsv',truth,['#CHROM','POS','REF','ALT','INFO'])
        gdna=[37,38,39,112,113,114,115] if roi==1 else [116,117,118,119]
        observed=[]
        for sample in gdna:
            for i,r in enumerate(rows[:40]):
                new=dict(r);new['Sample']=f'COMT_{sample}_R';new['record_id']=new['Sample']+':'+str(i+2);observed.append(new)
        sr.write_tsv(prepared/f'rows_{roi}.tsv',observed,fields)
    out=OUT/'RF'
    proc=subprocess.run(['C:/Program Files/R/R-4.6.1/bin/Rscript.exe','--vanilla',str(P/'apply_source_rf.R'),str(P.parent),str(prepared),str(out)],capture_output=True,text=True,timeout=300)
    (OUT/'run_log.txt').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf8')
    if proc.returncode:raise RuntimeError('Synthetic RF adapter failed; see preserved log')
    checks=0;expected={}
    for roi in (1,2):
        with (out/f'predictions_{roi}.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        grouped={}
        for r in rows:
            prob=float(r['probability_TRUE'])
            if not 0<=prob<=1:raise AssertionError('probability domain')
            if r['Prediction']!=str(prob>=.5).upper():raise AssertionError('class cutoff')
            key=(r['Sample'],r['VAR_ID'],r['TYPE'])
            if r['Prediction']=='FALSE':grouped[key]=grouped.get(key,0)+1
            checks+=2
        counts={}
        for (sample,vartype,typ),n in grouped.items():
            if typ=='SNP' and n==1:counts[vartype]=counts.get(vartype,0)+1
        expected[roi]={v:n for v,n in counts.items() if n>=(4 if roi==1 else 2)}
    with (out/'rejections.tsv').open() as f:actual=list(csv.DictReader(f,delimiter='\t'))
    byroi={roi:{r['VAR_ID']:int(r['false_gdna_samples']) for r in actual if int(r['ROI'])==roi} for roi in (1,2)}
    if byroi!=expected:raise AssertionError('Independent source rejection aggregation mismatch')
    checks+=1
    report={'status':'PASS_SYNTHETIC_ONLY','checks':checks,'invented_source_rows':440,'invented_simulation_rows':480,
            'real_source_rows_read':0,'source_algorithm_rejection_logic_independently_recomputed':True,
            'RF_adapter_sha256':sr.sha(P/'apply_source_rf.R'),'fixture_source_sha256':sr.sha(Path(__file__)),
            'outputs':{str(p.relative_to(OUT)).replace('\\','/'):sr.sha(p) for p in OUT.rglob('*') if p.is_file()}}
    (OUT/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='outputs'},indent=2))


if __name__=='__main__':main()
