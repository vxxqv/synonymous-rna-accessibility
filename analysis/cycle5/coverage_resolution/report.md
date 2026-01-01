# Frozen hypothetical coverage stress test

The selected design can lose coverage rapidly when missing fields are scattered across genotypes. Under the stipulated5% probability of independently missing each required field, only3.09% of development draws and3.81% of reserved draws meet the descriptive minimum. Prediction evaluation meets its minimum in2.19% and3.81%, respectively. These are results of a hypothetical missingness model, not estimates of missingness in the source experiment or evidence that the actual pilot will fail.

The original47 squares,94 diagonal pairs,181 shared genotype IDs,32 related seed groups and46 globally masked IDs were fixed. Prediction retains the original22 development and24 reserved structural cases, with21 and24 distinct target IDs. No case, corner, split, threshold or parameter was changed after simulation. Only the pinned split and mask metadata were read. No phenotype archive or outcomes were opened.

## Design and complete scenario results

`design_v1.json` was frozen before synthetic draws at SHA256 `90a3e1cdeb425cbaf482df4de0ecac3c701510b0fa95e6836918428752ba5363`. Each scenario uses10,000 fixed-seed draws, for130,000 draws total. Shared genotype IDs receive one shared missingness realization across all their squares. Coverage requires at least10 squares and8 related seed groups, separately in each arm.

- **A, field missingness:** each of seven required fields is independently missing for each genotype with probability p. Complete cases need all28 square fields. An individual square therefore survives with probability(1-p)^28, although squares sharing genotypes have dependent survival.
- **B, group absence:** each related seed group is unavailable with probability p. All its genotypes disappear together.
- **C, masked-target absence:** each globally masked genotype is unavailable with probability p. All unmasked genotypes remain complete. The structural multi-mask exclusion stays excluded from prediction.

The table reports the percentage of10,000 draws meeting coverage. “Known” means enough predictions can be made from the three known corners. “Evaluated” additionally requires the withheld target. “Descriptive” includes every original square with four complete corners.

| Model | p | Dev descriptive | Dev known | Dev evaluated | Reserved descriptive | Reserved known | Reserved evaluated |
|---|---:|---:|---:|---:|---:|---:|---:|
|A|0|100.00|100.00|100.00|100.00|100.00|100.00|
|A|0.01|99.90|99.99|99.90|99.99|100.00|99.99|
|A|0.05|3.09|17.96|2.19|3.81|27.18|3.81|
|A|0.10|0.00|0.01|0.00|0.00|0.01|0.00|
|A|0.25|0.00|0.00|0.00|0.00|0.00|0.00|
|B|0|100.00|100.00|100.00|100.00|100.00|100.00|
|B|0.10|100.00|100.00|100.00|100.00|100.00|100.00|
|B|0.25|99.18|99.15|99.15|99.25|99.25|99.25|
|B|0.50|58.06|57.92|57.92|59.77|59.77|59.77|
|C|0|100.00|100.00|100.00|100.00|100.00|100.00|
|C|0.10|100.00|100.00|100.00|100.00|100.00|100.00|
|C|0.25|99.91|100.00|99.90|99.98|100.00|99.98|
|C|0.50|73.15|100.00|71.65|83.66|100.00|83.66|

The percentages are Monte Carlo estimates, except for deterministic endpoint properties. Zero observed successes at a nonendpoint rate do not prove zero true probability. Likewise100% observed success at a nonendpoint rate does not prove certainty. Every scenario, including failures, is retained.

## Target resolution and conditional coverage

At A,p=0.05, the mean evaluated cases are5.2349 in development and5.7139 in reserved validation. The corresponding distinct targets are5.1144 and5.7139. The2.5th,50th and97.5th percentiles of evaluated case count are2,5,9 and2,6,10; the development distinct-target percentiles are2,5,9. These are simulated count-distribution summaries, not uncertainty intervals for a biological effect.

Prediction generation alone is not enough. Conditional on known-corner coverage passing at A,p=0.05, evaluation coverage passes in219/1796 draws, or12.1938%, in development and381/2718, or14.0177%, in the reserved arm. Conditional on descriptive coverage passing, the corresponding results are219/309, or70.8738%, and381/381, or100%. Reserved descriptive and evaluated cases are identical by design because all24 reserved squares are structurally prediction-eligible. Development has one prespecified prediction exclusion.

At C,p=0.50, all22 and24 known-corner predictions remain available in every draw, but an average10.9920 and12.0238 predictions are unscorable after target unsealing. Mean distinct evaluated targets are10.5019 and11.9762. Target finiteness must therefore remain an evaluation-stage result. It cannot be used to choose a more favorable set of predictions beforehand. Under whole-group absence B there are no predicted-but-unscorable cases because a group's known corners and target are available or absent together.

`results.json` includes every arm/scenario's square, group and distinct-target count distributions; means and sample standard deviations; quantiles at0,0.025,0.25,0.5,0.75,0.975,1 using the prespecified inverted-CDF convention; unconditional pass rates; both-arm joint pass rates; and both conditional pass rates with their exact conditioning denominators. A conditional denominator of zero gives null, not a zero rate. The compressed count files retain every draw's count and pass indicators without storing any biological values.

## Monte Carlo error and interpretation

The reported rate MCSE is sqrt(p_hat*(1-p_hat)/n), using10,000 for unconditional rates or the observed conditioning denominator for conditional rates. A count-mean MCSE is the sample standard deviation divided by sqrt(10,000). At A,p=0.05, the evaluation-rate MCSEs are0.1464 and0.1914 percentage points. The conditional evaluation-given-known MCSEs are0.7721 and0.6659 percentage points. A plug-in zero MCSE for zero observed events is not a claim that an event is impossible. These quantities describe finite simulation integration error under stipulated probabilities, not biological variation, measurement uncertainty or uncertainty in the unknown real missingness mechanism.

The scenario contrasts cannot be interpreted as an advantage of losing whole groups: they impose different missingness structures and rates on different units. No model was estimated or selected from phenotype data. In particular, these independent random-missingness models cannot correct selective assay failure, low-expression censoring, outcome-dependent missingness or shared assay covariance. Meeting coverage does not guarantee both classification labels, useful margins, reliable errors or successful prediction. This is not a phenotype power calculation.

The design consequence is to retain strict coverage reporting and the staged prediction seal, then let actual authorized coverage determine whether the frozen analysis can answer its question. The simulations provide no justification to relax the minimum, impute missing targets, change the split or replace a failed pilot with a favorable scenario.

## Verification and reproducibility

All six deterministic p=0/p=1 controls pass. At p=0, development has23 descriptive squares and22 prediction cases/21 distinct targets; reserved has24 in all three counts. At p=1, descriptive and evaluation counts are zero for every family. In C,p=1, every structural prediction is still made and every target is unscorable. A separately written Python set-based reference agrees on832 generated draws, making23,296 scalar count/coverage comparisons. Stage containment and prediction accounting are also checked for every draw: evaluation is a subset of both prediction and descriptive eligibility, and predicted count equals evaluated plus unscorable count.

The frozen implementation uses CPython3.12.14, NumPy2.5.4 and PCG64, with one documented seed per scenario and batches of1,000. `simulation_v1/execution_freeze.json` pins the source hash, exact input hashes and iteration orders before random draws. `simulation_v1/receipt.json` pins all generated artifacts. The result JSON SHA256 is `63acd3854727125af664322f7de0a4ebd0db9014b07e7cc6850ed91030c40a6e`.

```text
python simulate_coverage.py --design design_v1.json --output NEW_RESULT_DIRECTORY
```

The script refuses changed inputs, a changed design hash, a mismatched runtime, disabled assertions or an existing output directory. It has no phenotype ZIP reader. All existing pilot code and protocols remain unchanged.
