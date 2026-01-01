"""Interpretation-only extraction from authorized existing results and sequences."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import hashlib
import json

if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
root = Path.cwd()
out = root / 'work/research_cycle7/biological_interpretation'
dest = out / 'evidence.json'
if dest.exists():
    raise FileExistsError(dest)
paths = [
    'work/research_cycle7/results/robustness.json',
    'work/research_cycle7/results/scenario_geometry.json',
    'outputs/synonymous-rna-accessibility/analysis/cycle7/robustness_protocol.md',
    'outputs/synonymous-rna-accessibility/analysis/cycle7/scenario_geometry/protocol.md',
    'outputs/synonymous-rna-accessibility/analysis/cycle7/freeze.json',
    'outputs/synonymous-rna-accessibility/analysis/cycle6/fas_pilot/source_contract.json',
    'outputs/Research_Cycle_6/FAS_Results_and_Decision.md',
    'work/research_cycle7/recovery_identifiability/report.md',
    'work/research_cycle7/recovery_independent_review/review.md',
    'work/research_cycle7/robustness_result_review/review.md',
    'work/research_cycle7/scenario_result_review/observed_run/verification.json',
    'work/research_cycle6/fas_successor_metadata_gate/membership_receipt.json',
]
hashes = {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
a = json.loads((root/paths[0]).read_text())
d = json.loads((root/paths[1]).read_text())
seq = json.loads((root/paths[5]).read_text())['sequences']
epsilon = F(a['full']['additive']['parameter_exact'])
t = F(a['full']['multiplicative_squared']['parameter_exact'])
assert epsilon == F(519286149121,31050521652448)
assert t == F(9128860950384,8609574801263)
feasible = [x['omitted_vertex'] for x in a['all_single_vertex_deletions'] if x['certificate']['status']=='feasible']
assert feasible == [13,42,46]
prep = []
for x in a['all_single_preparation_deletions']:
    e = F(x['results']['additive']['parameter_exact'])
    m = F(x['results']['multiplicative_squared']['parameter_exact'])
    assert e >= epsilon and m >= t
    prep.append({'kind':x['kind'],'zero_based_index':x['omitted_zero_based'],
                 'epsilon':str(e),'gamma_squared':str(m),
                 'both_unchanged':e==epsilon and m==t})
sets = {label:{(p['s'],p['t']) for p in d['pairs'] if any(z['label']==label and z['nonempty'] for z in p['intervals'])}
        for label in ['fixed_route','adaptive_routes_only','mixed_scenario_connectivity','joint_obstruction_box_connected','upper_box_obstruction']}
assert {k:len(v) for k,v in sets.items()} == d['summary']['category_nonempty_pair_counts']
membership = []
for v in [13,16,42,46,47]:
    changes=[{'position1':i+1,'WT':x,'variant':y} for i,(x,y) in enumerate(zip(seq[46],seq[v])) if x!=y]
    membership.append({'vertex_zero_based':v,'distance':len(changes),'substitutions':changes,
                       'theoretical_radius2':len(changes)<=2,'observed_2016_coverage':None})
with localcontext() as c:
    c.prec=35
    gamma=(Decimal(t.numerator)/Decimal(t.denominator)).sqrt()
    approximate={'epsilon':str(Decimal(epsilon.numerator)/Decimal(epsilon.denominator)),
                 'gamma':str(gamma),'upward_percent':str(100*(gamma-1)),
                 'downward_percent':str(100*(1-1/gamma))}
evidence={'scope':'Interpretation and summary arithmetic, not independent full numeric verification',
          'input_hashes':hashes,'epsilon_exact':str(epsilon),'gamma_squared_exact':str(t),
          'gamma_enclosure':a['full']['multiplicative_squared']['gamma_enclosure'],
          'approximate':approximate,'feasible_row_omissions':feasible,
          'WT_identity_vertex':46,'preparation_omissions':prep,
          'unchanged_omission_radii_count':sum(x['both_unchanged'] for x in prep),
          'scenario_summary':d['summary'],'adaptive_mixed_overlap':len(sets['adaptive_routes_only']&sets['mixed_scenario_connectivity']),
          'adaptive_mixed_union':len(sets['adaptive_routes_only']|sets['mixed_scenario_connectivity']),
          'new_post_outcome_target_membership':membership,
          'new_target_sequence_only':True,'no_new_source_outcomes_read':True,
          'no_source_acquisition':True,'historical_git_reauthenticated':False,
          'report_sha256':hashlib.sha256((out/'report.md').read_bytes()).hexdigest(),
          'extractor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert '\u2014' not in (out/'report.md').read_text()
dest.write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'report_sha256':evidence['report_sha256'],'evidence_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
                  'new_target_radius2_members':sum(x['theoretical_radius2'] for x in membership),
                  'status':'COMPLETE_INTERPRETATION_AUDIT'}))
