# Cycle 6: completed frozen FAS pilot

The primary result is exactly certified ordinal infeasibility, with 72 finite envelopes, zero unknown vertices and zero conservative barrier pairs among 2,556 endpoint pairs. The corrected point has 73 barrier pairs. Independent implementations reconstruct all measurements and all 74,124 primary and secondary pair records. See `../../results/cycle6/decision.md` for the complete scope and retained limitations.

The full 72-state native included-protein class is fixed before deposited count access. The study measures relative included-RNA enrichment. It does not establish protein output, a splicing-specific mechanism, population confidence intervals or independent validation.

The primary test asks whether any categorical additive latent score and any common nondecreasing response can fit the empirical preparation envelopes. Feasible and infeasible results require exact rational certificates. Accessibility bounds retain unknown states. The two additive encodings are equivalent on this domain and are checked as such.

Run `python fas_pilot/verify_synthetic.py` and `python fas_pilot/verify_count_reader.py` from this directory. Both use synthetic data only and reject optimized Python. Dependencies are NumPy, SciPy and SymPy; the native-reference eligibility replay additionally uses Biopython. Tested versions are Python 3.12.14, SciPy 1.18.1 and SymPy 1.14.0.

The twelve source files, their immutable Git blob hashes and the fixed genotype list are in `fas_pilot/source_contract.json`. They originate at the pinned public `lehner-lab/Scaling_Law` commit. The analysis accepts an input directory containing those exact file basenames, validates every blob, and refuses to overwrite its output:

`python fas_pilot/run_headerless.py --directory ../../data/cycle6/fas_counts --output NEW_RESULT.json`

The script checks both freezes before opening counts. The original core failed at its header guard before row parsing. The separately reviewed format-only adapter preserves all twelve original headerless blobs and the first reference record. The median-input point is a corrected headerless source-formula reconstruction, not authenticated historical reproduction. No scientific rule changed. All failures and reviews remain in `format_review/`.

No deposited outcomes were examined when the original protocol snapshot was committed. Published literature and five illustrative source-README counts were already exposed. The format diagnosis additionally exposed the first native-reference rows. This is a prospective analytical freeze within a retrospective, adaptively selected source study.

For independent offline replay from any directory, run `python PATH_TO_REPO/analysis/cycle6/result_reviews/replay.py --repo PATH_TO_REPO --out NEW_DIRECTORY_OUTSIDE_REPO`. The standard-library-only replay verifies 32 pinned inputs, all exact numerical comparisons, 30 certificates and every pair bound without Git metadata or network access. Its README specifies the complete archive inventory and guards. The full frozen result is stored losslessly as `results/cycle6/full_result_v1.json.gz`.

Figure source is `render_figure.py`; pass the compressed result with `--result` and a separate export directory with `--out`. Manuscript and figure exports are maintained outside this computational repository. The figure reports all 72 variants and all 27 scenarios. Source studies retain their own terms and are identified by DOI and immutable Git blobs in the contracts and audits.
