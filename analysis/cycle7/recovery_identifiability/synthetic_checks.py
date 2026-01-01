"""Synthetic exact arithmetic only, no source outcomes or model fits."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
y=[F(1),F(3,2),F(1,5),F(7)]
p=F(1,3);c=1/max(y);recoveries=[c*v for v in y]
intensities=[p*q for q in recoveries]
reconstructed=[z/intensities[0] for z in intensities]
if reconstructed!=y or not all(0<q<=1 for q in recoveries):
    raise RuntimeError('Constant-inclusion arbitrary-order construction failed')

# Scaling convention: all-coordinate lambda=2 can explain point(1,4)
# under the constant-field class, but anchored kappa=2 cannot.
lam=F(2);original=[F(1),F(4)];fprime=[F(2),F(2)]
d=[a/b for a,b in zip(original,fprime)]
f=[v/fprime[0] for v in fprime];b=[fprime[0]*v for v in d]
if b[0]!=1 or [u*v for u,v in zip(b,f)]!=original or max(b)/min(b)!=lam**2:
    raise RuntimeError('All-coordinate/pairwise-spread construction failed')
anchored_two_cannot_fit_constant=(original[1]/lam>1)
if not anchored_two_cannot_fit_constant:
    raise RuntimeError('WT-bound versus spread counterexample malformed')

# Per-preparation recovery is not the same as one persistent correction.
scenario_factors=[F(2),F(1,2)]
persistent_constant_inclusion_possible=len(set(scenario_factors))==1
if persistent_constant_inclusion_possible:
    raise RuntimeError('Persistent/preparation-specific counterexample malformed')

# An absolute epsilon does not imply a fixed multiplicative recovery bound.
eps=F(1,100);smallf=F(1,10**6);large_bias=(smallf+eps)/smallf
if abs((smallf+eps)-smallf)!=eps or large_bias<=1000:
    raise RuntimeError('Additive versus multiplicative example malformed')

receipt={'scope':'Synthetic exact arithmetic corroboration, not proof or observed-data analysis',
 'constant_inclusion_arbitrary_positive_order':{'target':list(map(str,y)),'reconstructed':list(map(str,reconstructed)),'primitive_recovery_probabilities_valid':True},
 'all_coordinate_expansion_vs_persistent_pairwise_spread':{'lambda':str(lam),'normalized_persistent_factors':list(map(str,b)),'spread':str(max(b)/min(b)),'anchored_factor_two_not_enough_for_this_toy_constant_model':True},
 'preparation_specific_adjustment_not_one_persistent_factor':True,
 'additive_example':{'absolute_error':str(eps),'true_score':str(smallf),'multiplicative_factor':str(large_bias)},
 'all_checks_pass':True,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'synthetic_checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(receipt,indent=2))
