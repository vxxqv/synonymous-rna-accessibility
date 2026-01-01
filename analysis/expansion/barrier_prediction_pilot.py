"""Two-stage, frozen synthetic barrier pilot; no measured biology is implied."""
from pathlib import Path
import sys,json,hashlib,argparse,itertools,csv,math
import numpy as np
import RNA
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).parent))
from accessibility import CODONS,translate,adjacency,compatible,components
from retention_filtration import filtration
from smoothness_controls import roughness
OUT=ROOT/'results/expansion/barrier_pilot'
GRAPH=[[1],[0,2],[1]];PAIR=(32,38);RT=.00198720425864083*310.15
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,m):
    if not x:raise ValueError(m)
def model():
    md=RNA.md();md.temperature=37.;md.dangles=2;md.noLP=0
    return md
def scores(seq,ensemble=False):
    values=[];energies=[];sequences=[]
    structure='.'*PAIR[0]+'('+'.'*(PAIR[1]-PAIR[0]-1)+')'+'.'*(len(seq)-PAIR[1]-1)
    for a,b in [('A','U'),('G','U'),('G','C')]:
        s=list(seq);s[PAIR[0]]=a;s[PAIR[1]]=b;s=''.join(s);sequences.append(s)
        fc=RNA.fold_compound(s,model())
        if ensemble:
            fc.pf();values.append(float(fc.bpp()[PAIR[0]+1][PAIR[1]+1]))
        else:
            e=float(fc.eval_structure(structure));energies.append(e);values.append(1/(1+math.exp(e/RT)))
    require(len(set(translate(s) for s in sequences))==1,'peptide mismatch')
    require(all(compatible(s,[PAIR]) for s in sequences),'pair mismatch')
    require(adjacency(sequences)==GRAPH,'complete three-state path')
    f=filtration(GRAPH,np.array(values),full=True)
    closed=max(0.,min(values[0],values[2])-values[1])
    require(abs(f['longest_lifetime']-closed)<1e-12,'independent closed-form lifetime check')
    require(abs(f['integrated_excess']-.5*closed)<1e-12,'independent integral check')
    return {'sequences':sequences,'scores':values,'isolated_hairpin_energies':energies,'filtration':f}
def contexts():
    rng=np.random.default_rng(2689);sense=sorted(c for aa in CODONS for c in CODONS[aa]);chosen=[];audit=[]
    prior=json.loads((ROOT/'data/native_windows.json').read_text());windows=set()
    for b in prior:
        windows.update(b['context'][j:j+90] for j in range(len(b['context'])-89))
    windows=sorted(windows)
    for index in range(24):
        for attempt in range(100):
            codons=list(rng.choice(sense,30));forward=index%2==0
            codons[4]='GGC' if forward else 'GGU';codons[10:13]=['AAG','CUA' if forward else 'UUA','UUU'];base=''.join(codons)
            previous=[x['baseline'] for x in chosen]
            maximum=max((sum(a==b for a,b in zip(base,w))/90 for w in windows+previous),default=0)
            audit.append({'scaffold':index+1,'attempt':attempt,'maximum_prior90nt_identity':maximum,'accepted':maximum<.8})
            if maximum<.8:break
        else:raise ValueError('No sequence passing prespecified identity gate')
        target=list(base);target[33]='U' if forward else 'C';target=''.join(target)
        control=list(base);control[14]='U' if forward else 'C';control=''.join(control)
        for variant in (target,control):
            require(translate(base)==translate(variant),'intervention not synonymous')
            require(sum(a!=b for a,b in zip(base,variant))==1,'mutation number')
        require(sum(c in 'GC' for c in target)==sum(c in 'GC' for c in control),'GC match')
        chosen.append({'scaffold':f'S{index+1:02d}','baseline':base,'targeted':target,'control':control,'peptide':translate(base),'direction':'C_to_U' if forward else 'U_to_C','targeted_position0':33,'control_position0':14,'pair0':list(PAIR),'maximum_prior90nt_identity':maximum})
    return chosen,audit
def prepare():
    require(not (OUT/'predictions.json').exists(),'Prediction file exists; preserve the original freeze')
    OUT.mkdir(parents=True,exist_ok=True);chosen,audit=contexts();rows=[]
    for c in chosen:
        c['predictions']={v:scores(c[v]) for v in ('baseline','targeted','control')};rows.append(c)
    out={'scope':'new synthetic computational predictions before ensemble evaluation','protocol_sha256':sha(Path(__file__).with_name('barrier_pilot_freeze.json')),'implementation_sha256':sha(Path(__file__)),'ViennaRNA':RNA.__version__,'rows':rows,'selection_audit':audit}
    (OUT/'predictions.json').write_text(json.dumps(out,indent=2,allow_nan=False))
    print(json.dumps({'scaffolds':len(rows),'frozen_predictions_sha256':sha(OUT/'predictions.json'),'ensemble_partition_functions_run':0},indent=2))
def ci(values):
    values=np.array(values);rng=np.random.default_rng(2689);draws=values[rng.integers(0,len(values),(10000,len(values)))].mean(axis=1)
    return {'mean':float(values.mean()),'lower':float(np.quantile(draws,.025)),'upper':float(np.quantile(draws,.975)),'scaffold_clusters':len(values)}
def validate():
    p=OUT/'predictions.json';pred=json.loads(p.read_text());require(pred['protocol_sha256']==sha(Path(__file__).with_name('barrier_pilot_freeze.json')),'protocol changed');require(pred['implementation_sha256']==sha(Path(__file__)),'implementation changed')
    import subprocess
    frozen=subprocess.check_output(['git','show','HEAD:results/expansion/barrier_pilot/predictions.json'],cwd=ROOT)
    require(frozen==p.read_bytes(),'Predictions must be committed before ensemble evaluation')
    rows=[];table=[];covered=[];contrasts=[]
    for r in pred['rows']:
        q={k:v for k,v in r.items() if k!='predictions'};q['predictions']=r['predictions'];q['validation']={};bpred=r['predictions']['baseline']['filtration']['longest_lifetime']
        for variant in ('baseline','targeted','control'):
            v=scores(r[variant],True);q['validation'][variant]=v;obs=v['filtration']['longest_lifetime'];estimate=r['predictions'][variant]['filtration']['longest_lifetime']
            perms=[];matched=[];values=np.array(v['scores']);target=roughness([(0,1),(1,2)],values)
            for perm in itertools.permutations(values):
                f=filtration(GRAPH,perm);perms.append(f['integrated_excess'])
                if abs(roughness([(0,1),(1,2)],perm)-target)<=.05*max(target,1e-14)+1e-14:matched.append(f['integrated_excess'])
            # Three noncollinear states admit an exact two-coordinate additive fit.
            x=np.array([[1,1,0],[1,0,0],[1,0,1]],float);coef=np.linalg.solve(x,values)
            require(np.max(np.abs(x@coef-values))<1e-12,'saturated additive benchmark')
            q['validation'][variant]['references']={'permutation_mean_I_E':float(np.mean(perms)),'roughness5percent_mean_I_E':float(np.mean(matched)),'accepted_label_permutations':len(matched),'additive_residual_df':0,'additive_coefficients':coef.tolist()}
            table.append({'scaffold':r['scaffold'],'variant':variant,'sequence':r[variant],'peptide':r['peptide'],'direction':r['direction'],'pair_i0':PAIR[0],'pair_j0':PAIR[1],'predicted_L':estimate,'observed_L':obs,'observed_I_E':v['filtration']['integrated_excess'],'predicted_barrier_ge001':estimate>=.01,'observed_barrier_ge001':obs>=.01,'AU_probability':v['scores'][0],'GU_probability':v['scores'][1],'GC_probability':v['scores'][2],'exact_roughness_reference_I_E':float(np.mean(matched))})
        bobs=q['validation']['baseline']['filtration']['longest_lifetime']
        for variant in ('targeted','control'):
            dp=r['predictions'][variant]['filtration']['longest_lifetime']-bpred;do=q['validation'][variant]['filtration']['longest_lifetime']-bobs
            result={'scaffold':r['scaffold'],'variant':variant,'predicted_change':dp,'observed_change':do,'covered':abs(dp)>=1e-4,'success':abs(do)>=1e-4 and dp*do>0}
            if result['covered']:covered.append(result)
            q.setdefault('intervention_outcomes',[]).append(result)
        contrasts.append(q['validation']['targeted']['filtration']['longest_lifetime']-q['validation']['control']['filtration']['longest_lifetime']);rows.append(q)
    error=[abs(r['predicted_L']-r['observed_L']) for r in table];zero=[r['observed_L'] for r in table]
    confusion={k:0 for k in ('TP','FP','TN','FN')}
    for r in table:confusion[('T' if r['predicted_barrier_ge001']==r['observed_barrier_ge001'] else 'F')+('P' if r['predicted_barrier_ge001'] else 'N')]+=1
    accuracy=sum(r['success'] for r in covered)/len(covered) if covered else None
    summary={'scope':'prospective synthetic computational feasibility pilot; no biological validation','prediction_sha256':sha(p),'scaffolds':len(rows),'landscapes':len(table),'partition_functions':len(table)*3,'barriers_ge001':sum(r['observed_barrier_ge001'] for r in table),'predicted_barriers_ge001':sum(r['predicted_barrier_ge001'] for r in table),'predicted_lifetime_MAE':float(np.mean(error)),'zero_lifetime_MAE':float(np.mean(zero)),'confusion':confusion,'covered_interventions':len(covered),'all_interventions':48,'direction_successes':sum(r['success'] for r in covered),'direction_accuracy':accuracy,'targeted_minus_control_lifetime':ci(contrasts),'observed_mean_I_E':float(np.mean([r['observed_I_E'] for r in table])),'mean_exact_roughness_reference_I_E':float(np.mean([r['exact_roughness_reference_I_E'] for r in table])),'barrier_generic_WC_over_GU_accuracy':sum(r['observed_barrier_ge001'] for r in table)/len(table),'additive_fits_saturated':len(table),'expansion_gates':{'coverage_at_least4':len(covered)>=4,'MAE_better_than_zero':float(np.mean(error))<float(np.mean(zero)),'direction_accuracy_at_least075':accuracy is not None and accuracy>=.75,'distinguishes_additive_explanation':False},'interpretation':'The three-state design cannot distinguish an additive field from any folding field. Do not expand this motif into a folding-specific engineering claim.'}
    (OUT/'validation.json').write_text(json.dumps({'summary':summary,'rows':rows},indent=2,allow_nan=False))
    with (OUT/'prediction_table.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=('prepare','validate'));args=a.parse_args();prepare() if args.stage=='prepare' else validate()
