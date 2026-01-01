"""Frozen Q3 pilot: prepare predictors first; evaluate development before holdout.

MFE structures are declared predictor features. Full partition functions are the
computational outcomes, and are never used by prepare(). No measured biology.
"""
from pathlib import Path
from itertools import product, permutations, combinations
import argparse
import csv
import hashlib
import json
import math
import sys
import time

import numpy as np
import RNA
from Bio.Seq import Seq
from Bio.Data import CodonTable

from design_controls import cube, adjacency, endpoint_lifetime, integral, mobius

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUT=ROOT/'results/programme/cube_pilot_v1'
PROTOCOL=HERE/'cube_pilot_protocol_v1.json'
POINTS=cube(3)
GRAPH=adjacency(POINTS)
ENDPOINTS=(0,6)  # 000,110, frozen before predictor and outcome calculations
RT=.00198720425864083*310.15
PAIR={'AU','UA','GC','CG','GU','UG'}
INPUTS=[PROTOCOL,Path(__file__),HERE/'design_controls.py',ROOT/'data/native_windows.json',
        ROOT/'results/expansion/barrier_pilot/predictions.json',
        ROOT/'results/programme/design_gate_v1/results.json']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition,message):
    if not condition:
        raise ValueError(message)


def engine():
    require(RNA.__version__=='2.7.2','Frozen engine version required')
    md=RNA.md();md.temperature=37.;md.dangles=2;md.noLP=0
    return md


def peptide(seq):
    return str(Seq(seq.replace('U','T')).translate(table=1))


def genotype_sequences(scaffold,variant):
    rows=[]
    for x in POINTS:
        seq=list(scaffold[variant])
        for p,bit in zip(scaffold['variable_positions0'],x):
            seq[p]='G' if bit else 'A'
        rows.append(''.join(seq))
    require(len(set(map(peptide,rows)))==1,'Peptide must be fixed')
    require(adjacency([tuple(s) for s in rows])==GRAPH,'Complete one-nucleotide cube required')
    require(all(s[i]+s[j] in PAIR for s in rows for i,j in scaffold['pairs0']),'Compatibility changed')
    return rows


def windows(strings,length):
    strings=[s[j:j+length] for s in strings for j in range(len(s)-length+1)]
    if not strings:
        return np.empty((0,length),dtype=np.uint8)
    return np.frombuffer(''.join(strings).encode('ascii'),dtype=np.uint8).reshape(-1,length)


def maximum_identity(seq,references,length):
    if not len(references):
        return 0.
    return max(float((references==np.frombuffer(seq[j:j+length].encode('ascii'),dtype=np.uint8)).mean(axis=1).max())
               for j in range(len(seq)-length+1))


def contexts():
    rng=np.random.default_rng(8191)
    sense=sorted(c.replace('T','U') for c in CodonTable.unambiguous_dna_by_id[1].forward_table)
    native=json.loads((ROOT/'data/native_windows.json').read_text())
    native_strings=[x['context'] for x in native]
    old=json.loads((ROOT/'results/expansion/barrier_pilot/predictions.json').read_text())
    old_strings=[r[v] for r in old['rows'] for v in ('baseline','targeted','control')]
    native100=windows(native_strings,100)
    prior90=windows(old_strings,90)
    chosen=[];attempts=[]
    for k in range(24):
        split='development' if k<12 else 'heldout'
        starts=[7,17,27] if k<12 else [8,20,32]
        forward=k%2==0
        for attempt in range(100):
            codons=list(map(str,rng.choice(sense,40)))
            codons[3]='AUU' if forward else 'AUC';codons[4]='UCU'
            for s in starts:
                codons[s:s+4]=['AAA','CUU' if forward else 'CUC','UUU','UUU']
            base=''.join(codons)
            target_site=3*starts[k%3]+5
            control_site=11
            target=list(base);target[target_site]='C' if forward else 'U'
            control=list(base);control[control_site]='C' if forward else 'U'
            variants=[base,''.join(target),''.join(control)]
            previous=windows([r[v] for r in chosen for v in ('baseline','targeted','control')],100)
            refs=np.vstack([native100,previous])
            max100=max(maximum_identity(s,refs,100) for s in variants)
            max90=max(maximum_identity(s,prior90,90) for s in variants)
            accepted=max100<.8 and max90<.8
            attempts.append({'scaffold':k+1,'attempt':attempt+1,'maximum100nt_identity':max100,
                             'maximum_old90nt_identity':max90,'accepted':accepted})
            if accepted:
                break
        else:
            raise ValueError('Identity gate failed after 100 attempts')
        for s in variants:
            require(peptide(s)==peptide(base),'Intervention is not synonymous')
        for s,p in zip(variants[1:],[target_site,control_site]):
            require(sum(a!=b for a,b in zip(s,base))==1,'One change per intervention')
            require(base[p-1]==s[p-1]=='U' and base[p+1]==s[p+1]=='U','Matched U/U flanks')
        require(sum(s in 'GC' for s in variants[1])==sum(s in 'GC' for s in variants[2]),'GC matching')
        starts_nt=[3*s for s in starts]
        chosen.append({'id':f'Q{k+1:02d}','split':split,'baseline':variants[0],
                       'targeted':variants[1],'control':variants[2],'peptide':peptide(base),
                       'variable_positions0':[s+2 for s in starts_nt],
                       'pairs0':[(s+i,s+10-i) for s in starts_nt for i in range(3)],
                       'unit_pairs0':[[(s+i,s+10-i) for i in range(3)] for s in starts_nt],
                       'target_position0':target_site,'control_position0':control_site,
                       'direction':'U_to_C' if forward else 'C_to_U',
                       'maximum100nt_identity':max100,'maximum_old90nt_identity':max90})
    return chosen,attempts


def structure(n,pairs):
    chars=['.']*n
    for i,j in pairs:
        chars[i]='(';chars[j]=')'
    return ''.join(chars)


def pairs_of(dbn):
    stack=[];pairs=[]
    for j,c in enumerate(dbn):
        if c=='(':
            stack.append(j)
        elif c==')':
            require(bool(stack),'Unbalanced structure')
            pairs.append((stack.pop(),j))
    require(not stack,'Unbalanced structure')
    return pairs


def summarize_field(scores):
    return {'scores':list(map(float,scores)),
            'B_primary':float(endpoint_lifetime(GRAPH,scores,*ENDPOINTS)),
            'I_E':float(integral(GRAPH,scores)),
            'mobius':{k:float(v) for k,v in mobius(POINTS,scores).items()},
            'all_endpoint_lifetimes':[{'a':i,'b':j,'B':float(endpoint_lifetime(GRAPH,scores,i,j))}
                                      for i,j in combinations(range(8),2)]}


def predict(scaffold,variant):
    sequences=genotype_sequences(scaffold,variant)
    subset_structures=[]
    for included in POINTS:
        pairs=[p for include,unit in zip(included,scaffold['unit_pairs0']) if include for p in unit]
        subset_structures.append(structure(120,pairs))
    mfes=[];independent=[];isolated=[]
    for seq in sequences:
        fc=RNA.fold_compound(seq,engine())
        dbn,energy=fc.mfe();mfes.append(dbn)
        energies=[float(fc.eval_structure(structure(120,unit))) for unit in scaffold['unit_pairs0']]
        isolated.append(energies)
        independent.append(sum(1/(1+math.exp(e/RT)) for e in energies)/3)
    library=sorted(set(subset_structures+mfes))
    pair_lists=[pairs_of(s) for s in library]
    retained=[len(set(p).intersection(map(tuple,scaffold['pairs0'])))/9 for p in pair_lists]
    scores=[];energy_table=[]
    for seq in sequences:
        fc=RNA.fold_compound(seq,engine())
        energies=[float(fc.eval_structure(dbn)) if all(seq[i]+seq[j] in PAIR for i,j in pp) else None
                  for dbn,pp in zip(library,pair_lists)]
        minimum=min(e for e in energies if e is not None)
        weights=[math.exp(-(e-minimum)/RT) if e is not None else 0. for e in energies]
        scores.append(sum(w*r for w,r in zip(weights,retained))/sum(weights))
        energy_table.append(energies)
    additive=summarize_field(independent)
    require(additive['B_primary']<1e-12,'Independent-hairpin field must have no barrier')
    require(all(abs(v)<1e-12 for k,v in additive['mobius'].items() if k.count('1')>=2),
            'Independent-hairpin field must be additive')
    return {'finite_competitor':summarize_field(scores),'independent_hairpins':additive,
            'MFE_structures':mfes,'structure_library':library,'library_energies_kcal':energy_table,
            'isolated_hairpin_energies_kcal':isolated,'sequences':sequences}


def prepare():
    require(not (OUT/'predictions.json').exists(),'Frozen predictions already exist')
    gate=json.loads((ROOT/'results/programme/design_gate_v1/results.json').read_text())
    require(gate['ranks']['Q3_pairwise']==7 and gate['ranks']['Q3_additive']==4,'Design gate')
    start=time.perf_counter();scaffolds,attempts=contexts()
    for row in scaffolds:
        row['predictions']={v:predict(row,v) for v in ('baseline','targeted','control')}
        for v in ('targeted','control'):
            a=row['predictions']['baseline']['isolated_hairpin_energies_kcal']
            b=row['predictions'][v]['isolated_hairpin_energies_kcal']
            require(np.max(np.abs(np.array(a)-np.array(b)))<1e-8,'Loop/control independent energies must match')
    result={'scope':'Predictor features only; no full partition functions',
            'input_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in INPUTS},
            'software':{'Python':sys.version.split()[0],'ViennaRNA':RNA.__version__,'NumPy':np.__version__},
            'partition_functions_run':0,'MFE_predictor_calculations':24*3*8,
            'seconds':time.perf_counter()-start,'scaffolds':scaffolds,'selection_attempts':attempts}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'predictions.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'scaffolds':24,'partition_functions_run':0,'seconds':result['seconds'],
                      'prediction_sha256':sha(OUT/'predictions.json')}))


def correlation_reference(scores):
    interior=[i for i in range(8) if i not in ENDPOINTS]
    edge_pairs=[(i,j) for i,neighbours in enumerate(GRAPH) for j in neighbours if i<j]
    def rough(q):
        return sum((q[i]-q[j])**2 for i,j in edge_pairs)/len(edge_pairs)
    observed=rough(scores);values=[];eligible=[];distinct=set()
    for permutation in permutations(interior):
        q=list(scores)
        for target,source in zip(interior,permutation):
            q[target]=scores[source]
        B=float(endpoint_lifetime(GRAPH,q,*ENDPOINTS));values.append(B)
        if abs(rough(q)-observed)<=.05*observed+1e-15:
            eligible.append(B);distinct.add(tuple(q))
    return {'permutations':len(values),'eligible':len(eligible),'distinct_eligible_fields':len(distinct),
            'unconditional_mean_B':float(np.mean(values)),
            'roughness_conditioned_mean_B':float(np.mean(eligible)) if eligible else None,
            'informative_coverage':len(eligible)>=20,'observed_edge_roughness':observed}


def evaluate(scaffold,variant):
    predicted=scaffold['predictions'][variant]
    scores=[]
    for seq in predicted['sequences']:
        fc=RNA.fold_compound(seq,engine());fc.pf();prob=fc.bpp()
        scores.append(sum(float(prob[i+1][j+1]) for i,j in scaffold['pairs0'])/9)
    observed=summarize_field(scores)
    anchor={i:scores[i] for i in (0,1,2,4)}
    additive=[anchor[0]+sum(x[k]*(anchor[1<<(2-k)]-anchor[0]) for k in range(3)) for x in POINTS]
    masked=[3,5,6,7]
    observed.update({'correlation_reference':correlation_reference(scores),
                     'four_anchor_predictions':additive,
                     'four_anchor_out_of_range':sum(not 0<=q<=1 for q in additive),
                     'four_anchor_masked_score_MAE':float(np.mean([abs(additive[i]-scores[i]) for i in masked])),
                     'competitor_masked_score_MAE':float(np.mean([abs(predicted['finite_competitor']['scores'][i]-scores[i]) for i in masked]))})
    return observed


def bootstrap(values):
    values=np.asarray(values,dtype=float)
    rng=np.random.default_rng(8192)
    means=values[rng.integers(0,len(values),(10000,len(values)))].mean(axis=1)
    return {'mean':float(values.mean()),'lower95':float(np.quantile(means,.025)),
            'upper95':float(np.quantile(means,.975)),'scaffold_clusters':len(values)}


def validate(split):
    destination=OUT/(split+'.json')
    require(not destination.exists(),'Preserve existing outcomes')
    pred=json.loads((OUT/'predictions.json').read_text())
    for relative,digest in pred['input_sha256'].items():
        require(sha(ROOT/relative)==digest,'Frozen input changed: '+relative)
    if split=='heldout':
        dev=json.loads((OUT/'development.json').read_text())
        require(all(dev['summary']['gates'].values()),'Development gate failed; heldout remains unopened')
    start=time.perf_counter();rows=[];table=[];errors=[];zero=[];deltas=[];score_differences=[]
    covered=[];positive=[]
    for scaffold in pred['scaffolds']:
        if scaffold['split']!=split:
            continue
        obs={v:evaluate(scaffold,v) for v in ('baseline','targeted','control')}
        row={'id':scaffold['id'],'outcomes':obs}
        predicted_B={v:scaffold['predictions'][v]['finite_competitor']['B_primary'] for v in obs}
        actual_B={v:obs[v]['B_primary'] for v in obs}
        pD=predicted_B['targeted']-predicted_B['control'];oD=actual_B['targeted']-actual_B['control']
        is_covered=abs(pD)>=.001
        success=is_covered and abs(oD)>=.0001 and pD*oD>0
        row['intervention']={'predicted_D':pD,'observed_D':oD,'covered':is_covered,'direction_success':success,
                             'target_change_from_baseline':actual_B['targeted']-actual_B['baseline'],
                             'control_change_from_baseline':actual_B['control']-actual_B['baseline']}
        rows.append(row);deltas.append(oD)
        if is_covered:covered.append(success)
        positive.append(any(q>=.005 for q in actual_B.values()))
        errors.append(float(np.mean([abs(predicted_B[v]-actual_B[v]) for v in obs])))
        zero.append(float(np.mean(list(actual_B.values()))))
        score_differences.append(float(np.mean([obs[v]['four_anchor_masked_score_MAE']-obs[v]['competitor_masked_score_MAE'] for v in obs])))
        for v in obs:
            table.append({'scaffold':scaffold['id'],'split':split,'landscape':v,
                          'predicted_B':predicted_B[v],'observed_B':actual_B[v],
                          'absolute_B_error':abs(predicted_B[v]-actual_B[v]),
                          'predicted_three_way':scaffold['predictions'][v]['finite_competitor']['mobius']['111'],
                          'observed_three_way':obs[v]['mobius']['111'],
                          'four_anchor_masked_score_MAE':obs[v]['four_anchor_masked_score_MAE'],
                          'competitor_masked_score_MAE':obs[v]['competitor_masked_score_MAE'],
                          'reference_eligible':obs[v]['correlation_reference']['eligible'],
                          'predicted_D':pD,'observed_D':oD,'intervention_covered':is_covered,'direction_success':success})
    accuracy=sum(covered)/len(covered) if covered else None
    gates={'positive_scaffolds_at_least4':sum(positive)>=4,'covered_contrasts_at_least4':len(covered)>=4,
           'B_MAE_below_zero_baseline':np.mean(errors)<np.mean(zero),
           'direction_accuracy_at_least075':accuracy is not None and accuracy>=.75}
    summary={'split':split,'scaffolds':len(rows),'landscapes':len(table),'partition_functions':len(table)*8,
             'seconds':time.perf_counter()-start,'positive_scaffolds_B_ge0005':sum(positive),
             'positive_landscapes_B_ge0005':sum(t['observed_B']>=.005 for t in table),
             'predicted_positive_landscapes_B_ge0005':sum(t['predicted_B']>=.005 for t in table),
             'B_MAE':float(np.mean(errors)),'zero_baseline_B_MAE':float(np.mean(zero)),
             'paired_zero_minus_predictor_B_MAE':bootstrap(np.array(zero)-np.array(errors)),
             'paired_anchor_minus_competitor_masked_score_MAE':bootstrap(score_differences),
             'covered_intervention_contrasts':len(covered),'direction_successes':sum(covered),
             'direction_accuracy':accuracy,'intervention_D':bootstrap(deltas),
             'reference_informative_landscapes':sum(r['outcomes'][v]['correlation_reference']['informative_coverage'] for r in rows for v in r['outcomes']),
             'gates':{k:bool(v) for k,v in gates.items()},
             'heldout_opened':split=='heldout'}
    if split=='heldout':
        median=float(np.median([r['outcomes'][v]['B_primary'] for r in dev['rows'] for v in r['outcomes']]))
        constant=[float(np.mean([abs(median-r['outcomes'][v]['B_primary']) for v in r['outcomes']])) for r in rows]
        improvement=bootstrap(np.array(constant)-np.array(errors))
        summary.update({'discovery_median_B':median,'constant_baseline_B_MAE':float(np.mean(constant)),
                        'paired_constant_minus_predictor_MAE':improvement})
        summary['gates'].update({k:bool(v) for k,v in {'MAE_ten_percent_better_than_zero':np.mean(errors)<=.9*np.mean(zero),
                                'MAE_ten_percent_better_than_constant':np.mean(errors)<=.9*np.mean(constant),
                                'positive_lower95_zero_improvement':summary['paired_zero_minus_predictor_B_MAE']['lower95']>0,
                                'positive_lower95_constant_improvement':improvement['lower95']>0}.items()})
    result={'prediction_sha256':sha(OUT/'predictions.json'),'protocol_sha256':sha(PROTOCOL),
            'summary':summary,'rows':rows,
            'interpretation':'Deterministic folding-model pilot; no measured biology. All scaffolds and unsuccessful interventions are retained. Bootstrap intervals describe generated-context variability.'}
    destination.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (OUT/(split+'_table.csv')).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','development','heldout'])
    stage=parser.parse_args().stage
    if stage=='prepare':prepare()
    else:validate(stage)
