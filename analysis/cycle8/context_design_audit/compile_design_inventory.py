"""Compile explicit unknown-coordinate and contrast inventories from metadata."""
import argparse, collections, hashlib, json
from pathlib import Path
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    m=json.loads((a.inputs/'mapping.json').read_text());c=json.loads((a.inputs/'certificates.json').read_text())
    r=json.loads((a.inputs/'result.json').read_text());byex=collections.defaultdict(set)
    for row in m:byex[row['exon']].update(x['pos'] for x in row['changes'])
    inventory=[]
    for row in m:
        coords=sorted(byex[row['exon']]);observed=row['observed_bases']
        inventory.append({'id':row['id'],'exon':row['exon'],'selected_exon_coordinates':coords,
          'bases_at_selected_coordinates':[observed.get(str(p)) for p in coords],
          'unknown_selected_coordinates':[p for p in coords if str(p) not in observed],
          'all_other_outside_word_coordinates':'unknown',
          'descriptor_present_alleles':row['alleles']})
    groups=collections.defaultdict(list)
    for i,row in enumerate(m):groups[row['intended_allele']].append(i)
    contrast_support={i for v in c['left_nullspace'] for i,x in enumerate(v) if x!='0'}
    contrasts=[]
    for group,inds in sorted(groups.items()):
        anchor=inds[0]
        for i in inds[1:]:
            contrasts.append({'group':group,'minus_state':m[anchor]['id'],'plus_state':m[i]['id'],
              'interpretation':'Difference between two combined intended-edit/marker local states, not an isolated intended-edit effect.',
              'both_in_sole_additive_constraint':anchor in contrast_support and i in contrast_support})
    loo=[{'id':row['id'],'within_context_leave_one_state_out_prediction_estimable_under_descriptor_additivity':i in contrast_support,
          'reason':'Nonzero coefficient in exact left-null relation' if i in contrast_support else 'No exact additive relation includes this row'} for i,row in enumerate(m)]
    mats={'strict_shared':c['matrix']+c['matrix']}
    X=c['matrix'];n=len(X[0]);pairlabels=c['interaction_pair_labels']
    Z=[x+[int(set(p).issubset(row['alleles'])) for p in pairlabels] for x,row in zip(X,m)]
    mats['context_specific_additive']=[x+[0]*n for x in X]+[[0]*n+x for x in X]
    z=len(Z[0]);mats['context_specific_all_observed_pairs']=[x+[0]*z for x in Z]+[[0]*z+x for x in Z]
    model_rows=[{'cell':cell,'id':row['id'],'exon':row['exon']} for cell in ['K562','KBM7'] for row in m]
    out={'status':'metadata_only','source_sha256':sha(Path(__file__)),
      'input_sha256':{p.name:sha(p) for p in [a.inputs/'mapping.json',a.inputs/'certificates.json',a.inputs/'result.json']},
      'states':len(m),'states_with_unknown_selected_coordinates':sum(bool(x['unknown_selected_coordinates']) for x in inventory),
      'within_group_contrast_basis_count':len(contrasts),'leave_one_state_out_estimable':len(contrast_support),
      'leave_one_state_out_nonestimable':len(m)-len(contrast_support),'outcomes_read':False}
    a.out.mkdir(parents=True)
    for name,obj in [('coordinate_inventory.json',inventory),('marker_contrast_inventory.json',contrasts),('leave_one_state_out_estimability.json',loo),('nested_matrices.json',{'rows':model_rows,'matrices':mats,'ranks':r['context_ranks']})]:
        (a.out/name).write_text(json.dumps(obj,indent=2)+'\n')
    out['outputs_sha256']={p.name:sha(p) for p in a.out.iterdir()}
    (a.out/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['input_sha256','outputs_sha256']}))
if __name__=='__main__':main()
