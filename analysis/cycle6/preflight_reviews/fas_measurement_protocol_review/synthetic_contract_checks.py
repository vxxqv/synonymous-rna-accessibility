"""Small exact checks of proposed measurement algebra, no real data."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
output = {'WT': 40, 'v': 20, 'other': 140}
inputs = [{'WT': 20, 'v': 10, 'other': 70},
          {'WT': 30, 'v': 20, 'other': 150},
          {'WT': 15, 'v': 25, 'other': 60}]
ratios = []
for inp in inputs:
    esv = F(output['v'],sum(output.values())) / F(inp['v'],sum(inp.values()))
    esw = F(output['WT'],sum(output.values())) / F(inp['WT'],sum(inp.values()))
    direct = F(output['v'],output['WT']) * F(inp['WT'],inp['v'])
    if esv/esw != direct:
        raise RuntimeError('Separate-input total cancellation failed')
    ratios.append(direct)
med = lambda xs: sorted(xs)[1]
median_source_ratio = F(output['v'],output['WT']) * med([F(i['WT'],sum(i.values())) for i in inputs]) / med([F(i['v'],sum(i.values())) for i in inputs])
if not min(ratios) <= median_source_ratio <= max(ratios):
    raise RuntimeError('Ratio-of-medians enclosure failed')

# Nondegenerate Q2: constant inclusion, genotype-dependent included-RNA recovery.
# The same q bias can persist across every preparation/normalization scenario.
inclusion = dict.fromkeys(['00','01','10','11'],F(1,2))
recovery = {'00':F(1), '01':F(1,5), '10':F(1,5), '11':F(1)}
observed = {v:inclusion[v]*recovery[v] for v in inclusion}
# A common nondecreasing response of a+b*x+c*y would require c<0
# from y00>y01, and c>0 from y11>y10, an impossibility.
contradictory_orders = observed['00'] > observed['01'] and observed['11'] > observed['10']
if not contradictory_orders or len(set(inclusion.values())) != 1:
    raise RuntimeError('Known-answer assay-confounding example malformed')

receipt = {
 'scope':'Synthetic exact algebra only; no observed source values or calibration literals',
 'separate_TR_total_cancellation':True,
 'synthetic_relative_scores':[str(x) for x in ratios],
 'median_normalization_relative_score_enclosed':True,
 'median_relative_score':str(median_source_ratio),
 'global_library_total_needed_for_separate_TR_relative_WT':False,
 'q2_constant_true_inclusion_recovery_bias_counterexample':{
   'inclusion':{k:str(v) for k,v in inclusion.items()},
   'recovery':{k:str(v) for k,v in recovery.items()},
   'observed':{k:str(v) for k,v in observed.items()},
   'codon_additive_common_nondecreasing_model_rejected_in_observed_score':True,
   'true_inclusion_interaction':False,
   'argument':'y00>y01 forces c<0; y11>y10 forces c>0 in z=a+b*x+c*y.'},
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'synthetic_checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(receipt,indent=2))
