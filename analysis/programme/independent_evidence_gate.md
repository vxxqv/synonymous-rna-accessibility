# Independent measured evidence: bounded eligibility decision

No experimental outcome analysis was conducted in this research cycle. Article-level findings encountered during source discovery are recorded as known information. Numerical mutant fitness cells were not extracted or examined in the new GFP source audit.

## Yeast GFP/URA3 source

Primary study: Zhang et al., *Mistranslation suppresses mistranscription in eukaryotes*, Nature Communications 17, 3181 (2026), https://doi.org/10.1038/s41467-026-69969-x . The article identifies NCBI BioProject PRJNA1347074 and a public analysis repository. These measurements concern function/fitness, not direct RNA base-pair occupancy.

The author repository exists; its default branch is `master`. A first query of the nonexistent `main` branch returned 404. This was resolved by reading repository metadata, not treated as dataset unavailability. Pinned source commit: `2c3f85db7703d1f717accc351901930eb243c818`. Its `data.zip` is 516,072 bytes, SHA-256 `a06e3ba6ecceecdf67abcd1b82c00643178b9cf2c09bdd78967167a44dfe75f3`. It contains per-figure RData summaries. The repository does not declare a code license in the inspected metadata; none of its source code or data files is redistributed here.

The primary article links both the source workbook (51 sheets; 3,771,260 bytes; SHA-256 `ea6e6843227339de32c3baa37a20739a4cc39aea60745027e0de36c32183c51c`) and a distinct Supplementary Data 1–7 workbook. The latter is 7,660,811 bytes, SHA-256 `7b612a56351fd4815397df67325fef5e94e7adf3910303a03114dc9aacb4a1b6`. Both were inspected for dimensions, text headers and genotype labels only. Source workbooks remain outside the computational repository.

Exact supplementary URL:
https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69969-x/MediaObjects/41467_2026_69969_MOESM3_ESM.xlsx

Supplement S7 has 50,690 distinct genotype labels: 4,068 single-substitution and 46,622 double-substitution labels. All labels parsed; none contains three substitutions. The columns distinguish GFP/URA3, two promoters and several media. These are conditions, not independent genes or studies.

**Decisive Q3 gap.** A complete three-site binary cube contains a genotype bearing a nonreference allele at each of its three varied positions. It therefore needs at least one genotype with at least three substitutions relative to the reference, even if some alleles at a site are both nonreference. A table restricted to single and double substitutions cannot contain a complete three-site cube. This excludes S7 as complete Q3 validation material without reading its numerical fitness values. It does not prove that all raw sequencing material lacks higher mutants.

**What remains feasible.** This source may support a separately frozen two-site functional study. Before doing so, recover the exact full construct sequences and gene-specific position maps, establish wild-type normalization and replicate/error semantics, check all four genotypes and nonmissing measurements in each prespecified condition, and compare appropriate global-epistasis and expression baselines. Shared single mutants and one reporter gene create dependence. No structural-mechanism claim follows without an appropriate structural assay. The numerical outcomes remain unexamined by this project and can remain held out while that eligibility work is completed.

## Other sources and limits

Lind and Andersson (2013), https://doi.org/10.1371/journal.pone.0063373 , reports 18 selected synonymous single mutants and nine doubles in rpsT, including failed compensatory predictions. Its published outcome summary was already read. A reanalysis would be retrospective, selected, largely one-gene evidence; it cannot supply unseen three-way validation or measured structural occupancy.

The earlier FoldAtlas recovery audit remains valid: the historical structure/probing database export was not recovered, although code and historical annotation sources were found. The earlier Cambray source audit establishes accessible reporter phenotypes and some observed synonymous edges, but not experimental target-pair occupancies or complete construct closure. These are specific gaps, not assertions that all public biological evidence is unavailable.

The independent biological validation branch is therefore incomplete for the present Q3 mechanism. The next evidence route should explicitly choose between a two-site measured functional question and a new complete-cube structural/functional experiment. It must not present a folding-engine comparison or synthetic-context holdout as biological validation.
