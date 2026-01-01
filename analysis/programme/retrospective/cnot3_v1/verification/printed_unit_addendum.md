# Rendered figure-axis unit check

**Material finding: the reanalysis explicitly prints minutes.** I visually inspected the actual official Figure1C/D and Figure2C/D images, then individually inspected saved enlarged crops. All four lifetime axes print `(min)`. This closes the earlier gap concerning **what the published figures label**. It does not independently authenticate the original fitting/export/plot conversion that connects the stored TableS2 numbers to that label.

| Primary figure and panel | Exact visible axis label | Evidence retained |
| --- | --- | --- |
| Carneiro Figure1C | `mRNA half-life (min)` | `Carneiro_Figure1_C_axis_evidence.png` |
| Carneiro Figure1D | `ΔCNOT3 mRNA half-life (min)` | `Carneiro_Figure1_D_axis_evidence.png` |
| Carneiro Figure2C | `mRNA half-life (min)` | `Carneiro_Figure2_C_axis_evidence.png` |
| Carneiro Figure2D | `ΔCNOT3 mRNA half-life (min)` | `Carneiro_Figure2_D_axis_evidence.png` |

Here ΔCNOT3 denotes the perturbed condition. The displayed quantity in panelsD is its half-life, not a subtraction of two half-lives.

The official images were selected by their exact names from the previously authenticated primary PMC metadata, not by guessing media sequence numbers:

- [Reanalysis Figure1, g001](https://pmc-oa-opendata.s3.amazonaws.com/PMC13398649.1/FEB4-9999-0-g001.jpg?md5=95c0a9746cb44896f4351a342ede3d0e), SHA256`efeff4a0ef43e49822aa0a9adc9d5738243140a0c3b1b54b210d27db143e2520`, original709x858pixels.
- [Reanalysis Figure2, g004](https://pmc-oa-opendata.s3.amazonaws.com/PMC13398649.1/FEB4-9999-0-g004.jpg?md5=b4dd3fdea6199adfac66add2b981444f), SHA256`b315a89f9c580592ede3eeb0ea8dff46afe43bc71ee83abd00b1a52fde423cd9`, original709x899pixels.
- [Original Zhu Figure2, f0003](https://pmc-oa-opendata.s3.amazonaws.com/PMC11583848.1/nihms-2027590-f0003.jpg?md5=2c72f7e5001864cecc6da1d1a007ac8a), SHA256`5de865bfea0e27463ed04244f66adb1bb8764a953c9fbd7357aab92aed4d801a`, original1590x2040pixels.

Every retrieval returned200 and exactly matched its metadata MD5. Crop coordinates and hashes are in `crop_coordinates.json`; originals are unchanged.

## Original-study comparison and relation to TableS2

Original Zhu Figure2A labels `log2 FC_Half-life(sgCNOT3/sgNT)` and is a dimensionless half-life ratio. Figure2B labels its horizontal axis `Time after dox addition (hours)` and its vertical axis `Percent remaining mRNA`. That is a distinct reporter-assay panel. It does not calibrate the native SLAM-seq half-life columns. Original Figure2C/D are enrichment/rank panels, not absolute half-life axes.

The reanalysis's TableS2 documentation identifies `HL_WT` and `HL_CNOT3KO` as the original study's wild-type and CNOT3-knockout half-lives. Its Figure1C/D and Figure2C/D display those corresponding concepts and explicitly claim minutes. It is therefore justified to report **the reanalysis's published minute-labelled half-life scale**. It is not justified to say no published unit label exists.

The remaining narrower gap is numerical provenance: no fitting/export or plotting code was inspected here to establish whether TableS2 stores those units directly or whether a conversion is applied before plotting. The stored headers/description contain no unit, and the original methods give chase times in hours without documenting that conversion. Hours for sampling and minutes for a plotted fitted lifetime are not intrinsically contradictory; conversion can occur. This check neither establishes a label error nor silently applies a factor60 correction to the data.

For protocol purposes, retain a distinction between the verified printed unit claim and unauthenticated source-to-column calibration. A dimensionless KO/WT response or ranks of a reciprocal-HL contrast under one declared common positive scale remain usable for the limited descriptive question. A physically calibrated per-minute or per-hour degradation-rate magnitude requires the missing linkage or an explicitly declared reliance on the published label. A common scale conversion leaves ranks unchanged but changes rate magnitudes. Growth dilution, fit uncertainty and causal pathway identification remain unresolved by these images.

## Scope and exposure

Full published figure composites were displayed during label inspection, followed by the four relevant panel crops. Their aggregate boxes, distributions and annotations were visible; no values were digitized, estimated, fitted or imported into a model. The original composite also contains other published experimental panels, which were not used for this unit conclusion. No gene-level XLSX half-life, association, conversion-ratio or quality outcome was read. No spreadsheet, analysis, repository or manuscript file was changed. Only this unit-check directory was written. This bounded check is complete; no other figure or unit claim is certified.
