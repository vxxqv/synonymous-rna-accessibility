"""Run a frozen, complete finite-scenario quantifier follow-up."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,gzip,hashlib,importlib.util,json

def main():
    if not __debug__:raise RuntimeError('Optimized Python prohibited')
    p=argparse.ArgumentParser();p.add_argument('--freeze',required=True,type=Path)
    p.add_argument('--result',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    if args.out.exists():raise FileExistsError('Refusing existing result')
    freeze=json.loads(args.freeze.read_text())
    for name,sha in freeze['files'].items():
        if hashlib.sha256((args.freeze.parent/name).read_bytes()).hexdigest()!=sha:raise ValueError('Changed frozen source: '+name)
    raw=args.result.read_bytes()
    if args.result.suffix=='.gz':raw=gzip.decompress(raw)
    if hashlib.sha256(raw).hexdigest()!=freeze['primary_result_sha256']:raise ValueError('Changed primary result')
    d=json.loads(raw);v=d['vertices_data']
    if len(v)!=72 or [z['vertex'] for z in v]!=list(range(72)):raise ValueError('Unexpected frozen vertices')
    if any(len(z['scenario_ratios'])!=27 for z in v):raise ValueError('Incomplete scenario array')
    sequences=[z['sequence'] for z in v]
    edges=[(i,j) for i,j in combinations(range(72),2) if sum(a!=b for a,b in zip(sequences[i],sequences[j]))==1]
    if len(edges)!=204:raise ValueError('Graph mismatch')
    fields=[[None if z['scenario_ratios'][r] is None else F(z['scenario_ratios'][r]) for z in v] for r in range(27)]
    spec=importlib.util.spec_from_file_location('scenario_followup',args.freeze.parent/'scenario_geometry/scenario_geometry.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    result=core.analyze(fields,edges)
    if result['pair_count']!=2556:raise ArithmeticError('Lost endpoint pairs')
    if not result['unknown_cells']:
        if len(d['bounded_pairs'])!=2556:raise ValueError('Incomplete primary pair inventory')
        for pair,old in zip(result['pairs'],d['bounded_pairs']):
            if (pair['s'],pair['t'])!=(old['a'],old['b']):raise ArithmeticError('Pair order mismatch')
            for key,oldkey in [('E','endpoint_lower'),('C_fix','bottleneck_lower'),('C_upper','bottleneck_upper')]:
                if pair[key]!=F(old[oldkey]):raise ArithmeticError('Primary bound changed')
        counts={label:sum(z['nonempty'] for p in result['pairs'] for z in p['intervals'] if z['label']==label) for label in core.LABELS}
        summary={'category_nonempty_pair_counts':counts,
                 'all_scenario_obstruction_pairs':sum(p['all_scenario_obstruction']['nonempty'] for p in result['pairs']),
                 'upper_box_obstruction_pairs':sum(p['upper_box_obstruction']['nonempty'] for p in result['pairs']),
                 'interpretation':'Counts may overlap because different threshold ranges have different classifications. Not independent confirmations or probabilities.'}
    else:summary={'status':'bounded_unknown','exact_classification_not_emitted':True}
    result['summary']=summary;result['primary_result_sha256']=hashlib.sha256(raw).hexdigest()
    result['scope']='Exploratory exposed-source sensitivity with all27 dependent joint scenarios and common absolute thresholds; original primary unchanged.'
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(core.json_ready(result),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':main()
