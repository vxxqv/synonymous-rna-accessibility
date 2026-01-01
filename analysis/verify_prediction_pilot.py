"""Independent numerical replay of every frozen synthetic prediction and outcome."""
from pathlib import Path
import json,csv,hashlib,math
import RNA
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/expansion/barrier_pilot'
def require(x,message):
    if not x:raise ValueError(message)
def main():
    predictions=json.loads((OUT/'predictions.json').read_text());validation=json.loads((OUT/'validation.json').read_text());table=list(csv.DictReader((OUT/'prediction_table.csv').open()))
    require(hashlib.sha256((OUT/'predictions.json').read_bytes()).hexdigest()==validation['summary']['prediction_sha256'],'frozen prediction digest')
    for key,file in [('protocol_sha256','analysis/expansion/barrier_pilot_freeze.json'),('implementation_sha256','analysis/expansion/barrier_prediction_pilot.py')]:require(predictions[key]==hashlib.sha256((ROOT/file).read_bytes()).hexdigest(),key)
    rows={r['scaffold']:r for r in validation['rows']};checks=0
    for r in predictions['rows']:
        require(r['maximum_prior90nt_identity']<.8,'sequence-related gate')
        for variant in ('targeted','control'):
            require(str(Seq(r['baseline']).translate())==str(Seq(r[variant]).translate()),'independent translation')
            dif=[i for i,(a,b) in enumerate(zip(r['baseline'],r[variant])) if a!=b];require(dif==[33 if variant=='targeted' else 14],'mutation definition')
        require(sum(b in 'GC' for b in r['targeted'])==sum(b in 'GC' for b in r['control']),'GC matched control')
        for variant in ('baseline','targeted','control'):
            p=r['predictions'][variant];v=rows[r['scaffold']]['validation'][variant];energies=[];scores=[]
            for i,s in enumerate(p['sequences']):
                require(str(Seq(s).translate())==r['peptide'],'genotype peptide')
                md=RNA.md();md.temperature=37.;md.dangles=2;md.noLP=0;fc=RNA.fold_compound(s,md);energy=fc.eval_structure('.'*32+'(.....)'+'.'*51);energies.append(1/(1+math.exp(energy/(.00198720425864083*310.15))));fc.pf();scores.append(fc.bpp()[33][39]);checks+=1
            for x,y in zip(energies,p['scores']):require(abs(x-y)<1e-12,'isolated prediction replay')
            for x,y in zip(scores,v['scores']):require(abs(x-y)<1e-12,'ensemble replay')
            for field,expected in ((energies,p),(scores,v)):
                lifetime=max(0,min(field[0],field[2])-field[1]);require(abs(lifetime-expected['filtration']['longest_lifetime'])<1e-12,'closed-form lifetime');require(abs(.5*lifetime-expected['filtration']['integrated_excess'])<1e-12,'closed-form integral')
            t=next(x for x in table if x['scaffold']==r['scaffold'] and x['variant']==variant);require(abs(float(t['observed_L'])-v['filtration']['longest_lifetime'])<1e-12,'prediction-table field')
    require(len(table)==72 and checks==216,'complete pilot coverage')
    print(json.dumps({'status':'passed','independently_replayed_genotypes':checks,'complete_prediction_table_rows':len(table),'interpretation':'numerical replay, not experimental validation'}))
if __name__=='__main__':main()
