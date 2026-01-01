# Original-CT PAIR-MaP validation

This branch uses the original processed data of Mustoe et al. (2019), *PNAS*, DOI [10.1073/pnas.1905491116](https://doi.org/10.1073/pnas.1905491116). The versioned data record [10.6084/m9.figshare.9978596.v1](https://doi.org/10.6084/m9.figshare.9978596.v1) is explicitly CC0. The source archive SHA-256 is `8320f005b563333690b27078a6e744879b713f602f1371a61bcf369c945ca4e9`.

Only the rpsM/rpsB source subsets are included. Per-file checksums are in `data/pairmap/source_manifest.json`. Original CT partner identities are retained, including crossings in the manually revised reference structures. No pseudoknot energies are computed.

## Coordinates and evidence

All source-table coordinates are 1-based inclusive. Computational windows and the compatibility engine use zero-based positions; `variable_start0` identifies that convention explicitly. Codon phase in the position map is 1, 2, or 3.

Both full source sequences map uniquely and exactly to the study-matched GenBank reference U00096.2, annotation dated 26 February 2013:

| Locus | Source length | Genomic interval | Strand | Coding positions in source |
|---|---:|---|---|---|
| rpsM / b3298 | 241 nt | 3440392–3440632 | minus | 140–241 |
| rpsB / b0169 | 204 nt | 189713–189916 | plus | 162–204 |

The rpsM fragment covers 34 complete coding codons. The rpsB fragment covers 14 complete coding codons and one incomplete terminal codon. Initiating codons remain fixed; incomplete codons are never varied. The eligible variable windows do not overlap any other annotated CDS on either strand.

PAIR-MaP records identify interactions between three-nucleotide windows. They are duplex-sensitive evidence, not direct identification of every individual base pair. Principal coding interactions are rpsM windows (182,209), (184,207), (188,202), and (189,201) in both conditions, and rpsB (170,179) in-cell. Every listed interaction corresponds to a complete three-base-pair antiparallel segment in the respective deposited CT. The source statistics and the model correspondence are retained separately.

## Frozen computation and outcome

Selection considered every complete five-codon block after the fixed initiator, retained at least 80% finite HQ-profile coverage in both conditions, capped the complete Cartesian sequence space at 1,024 variants, and required at least three target pairs touching the variable block in each condition. Full 241/204-nt contexts remain fixed outside the block. This yields 20 windows: 17 rpsM and 3 rpsB, with 6,800 Cartesian variants per condition before compatibility filtering.

All 40 condition/window graphs were enumerated completely. ViennaRNA 2.7.2 used 37 degrees Celsius, dangles 2, and noLP 0, without transferring wild-type experimental restraints to mutants. Filtering used mean target-pair probability thresholds 0.1, 0.25, 0.5, and 0.75 on exactly matched compatibility endpoints. This is a model of pair retention, not measured mutant fitness or folding kinetics.

Of 160 filtered graphs, 150 are nonempty and 10 empty. No nonempty graph has additional fragmentation beyond compatibility. All 17 rpsM compatibility spaces are connected; all three rpsB spaces have same-component probability 5/9. The condition-specific targets give the same compatibility probabilities in these selected windows. These are two loci with overlapping windows; the 20 windows are not independent biological replicates and no population significance claim is made.

A defined scan found no antiparallel run of at least three base pairs with a constant partner-register shift of one, two, or three nucleotides between the deposited condition CTs. Consequently these loci do not validate the restricted coding-register theorem. The revised rpsM reference contains a supported kissing-loop interaction omitted by the automated condition CTs; its crossing pair table is preserved as a separate structural hypothesis.

## Reproduction

From the repository root, with the repository environment active:

```text
python analysis/pairmap/fetch_sources.py
python analysis/pairmap/pairmap_audit.py
python analysis/pairmap/pairmap_register.py
python analysis/pairmap/select_pairmap_windows.py
python analysis/pairmap/run_pairmap_windows.py
python analysis/pairmap/summarize.py
```

The runner uses two workers and resumes complete existing result files. `results/pairmap/` retains the complete genotypes, graph edges, component sizes, folds, energies, and model pair probabilities. The flat result source table and summary are in `data/pairmap/`. Input and engine details are in `run_manifest.json`.
