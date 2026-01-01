# General-graph exact verification

From the repository root, with NumPy installed:

```text
python analysis/expansion/verification/run_verification.py --outdir results/expansion/new_verification_run
```

Use a fresh output directory. The driver records source hashes before/after,
software versions, commands, exit statuses, runtimes, code revision and result
hashes. It does not rerun the native folding study. A changed input invalidates
the run. The authoritative release receipt is
`results/expansion/verification_release/manifest.json`. The older
`results/expansion/verification_final/manifest.json` is retained as a historical
pilot receipt, before the explicit optimized-Python guard was installed.
Every verifier entry point rejects `python -O` before writing results.
The separate feasibility release receipt is
`results/expansion/feasibility_verification_release_v2/manifest.json`.

`arbitrary_pair_graph_oracle.py` covers every nonempty instance among
13,107,360 possible one/two-codon peptide and simple-pair-graph combinations.
Every compatible genotype x is visited, and every subgraph of its allowed-pair
edges is visited; this proves coverage of all nonempty instances. Empty cases
are counted by subtraction. These are abstract compatibility graphs: the
universe includes crossing edges, pairs within codons and arbitrary degree.
It imposes no hairpin minimum or physical multi-partner bonding model.

`leucine_hardcore_volume.py` checks projected domains and Leu-only activity-two
weights, plus an independent rational-polynomial tree-poset counter.
`review_score_controls.py` checks the conditional vertex null, score-filtration
ties/barcodes/area against an independent BFS implementation, and the additive
superlevel-connectivity theorem on all four-vertex graph/sign cases.

The public compiler is `../graph_compiler.py`; it emits witness/predecessor
certificates and separately reports resource-limited path/volume counts.
`../validate_graph_compiler.py` compares it with complete genotype graphs.

## Quantitative figure-source fields

`arbitrary_pair_graph_results.json` provides:

- `distance_counts`: connected/separated unordered endpoint-pair counts by
  nucleotide Hamming distance;
- `mutation_order_fraction_counts`: exact valid-order fraction among k!
  direct mutation orders, accumulated over connected endpoint pairs;
- `distance_order_height_counts`: joint exact counts by distance, valid/all
  orders and dependency height;
- `peptide_component_summaries`:420 one/two-amino-acid rows. Filter
  `codons == 2` for the complete20×20 peptide heat map. Maximal component and
  effective-component counts are across every compatible simple pair graph;
  fractions are exact strings.

Endpoint-pair histograms weight every eligible instance and its unordered
endpoint pairs once. They are not uniform peptide counts or native frequencies.
Plot conditional distance strata rather than interpreting the abundant
distance-one cases as a biological effect. Whole-theorem proofs remain
separate from these finite assertions; successful tests do not establish
universal conclusions or novelty.
