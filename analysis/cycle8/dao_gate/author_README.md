


# U-rich elements drive pervasive cryptic splicing in 3’ UTR massively parallel reporter assays 

Corresponding author: Anthony Mustoe - anthony.mustoe@bcm.edu

Last updated: 07/03/2025

### Analysis codes (organized by figure)
  - `process_raw_count_to_splicing_efficiency.ipynb`
      - `Input`: raw count by categories for each reporter (`ptreseq_raw_count/`). See **Method** for a description of each category.
      - `Output`: median observed splicing fraction (`ptreseq_splicing_quantification/`)
  - `F1_splicing_quantification.ipynb`       
  - `F1S_F2s_palindrome_and_splicing_reproducibility.ipynb`     
  - `F2_donor_acceptor_identification.ipynb` 
  - `F3_AU_rich_element_splicing.ipynb`      
  - `F4_modelling_splicing_impact.ipynb`     
  - `F4_ptreseq_splicing_reanalysis.ipynb`   
  - `F5_cryptic_splicing_mpra.ipynb`        
     
### `Data/`
  - tape station electropherogram for sequencing library
    - Full PTRE-seq library: `2023-08-04 - 12-17-54-D1000_Electropherogram.csv` 
    - Full PTRE-seq library (alternative library prep method): `2025-03-16 - DNA RNA gel_Electropherogram.csv`
    - "No spacer" library: `2020-11-16 - HELA no spacer_Electropherogram` & `2023-07-21 - DNA no spacer_Electropherogram`
  - `HELA-1_gap_sequences.txt.gz`: nucleotide sequences at boundaries of internal deletions in HeLa RNA sample
  - `single_reporter_validation_quantification.xlsx`: westernblot and RT-PCR quantification for single reporter validation
  - Maximum Entropy score of annotated human splice sites (Ensembl hg38 assembly)
      - `human_5UTR_5ss_MAXENT.txt`
      - `human_5UTR_3ss_MAXENT.txt`
      - `human_CDS_5ss_MAXENT.txt`
      - `human_CDS_3ss_MAXENT.txt`
      - `human_3UTR_5ss_MAXENT.txt`
      - `human_3UTR_3ss_MAXENT.txt`
      - `human_5ss_sequence.fa` and `human_3ss_sequence.fa` (sequence for all annotated human splice sites)
  - Measurements from original PTRE-seq publication
      - `sup2_readscount_plasmid.xlsx`: read count for input DNA plasmid library
      - `sup3_readscount_HeLa_total.xlsx`: read count for total RNA
      - `sup4_readscount_HeLa_polysome.xlsx`: read count for polysome-associated RNA
  - Nanopore sequencing of PTRE-seq DNA library
      - `plasmidsaurus_length_distribution.txt`
      - `long_read_seq_report.txt`

### `splice_prediction/`
 - Scripts for splicing prediction using spliceAI or Pangolin
      - Example usage:
   
   `python spliceai_script.py ptreseq_mpra_seq.fa ptreseq_full_transcript_prediction_best_sites.txt`
 - Sequences for 3'UTR MPRAs and other reporter genes
      - `ptreseq_mpra_seq.fa.gz`
      - `griesemer_mpra_seq.fa.gz`
      - `zhao_mpra_seq.fa.gz`
      - `siegel_mpra_seq.fa.gz`
      - `MapUTR_full_seq.fa.gz`
      - `reporter_genes_seq.fa`
 - SpliceAI prediction for 3'UTR MPRAs
      - `ptreseq_full_transcript_prediction_best_sites.txt`
      - `griesemer_spliceai_full_transcript_prediction_best_sites.txt`
      - `zhao_spliceai_full_transcript_prediction_best_sites.txt`
      - `siegel_spliceai_full_transcript_prediction_best_sites.txt`
      - `fu_spliceai_full_transcript_prediction_best_sites.txt`
  - XSTREME motif analysis for spliced reporters in 3'UTR MPRAs
      - `griesemer_motif_analysis.txt` and `griesemer_motif_analysis.html`
      - `siegel_motif_analysis.txt` and `siegel_motif_analysis.html`
      - `zhao_motif_analysis.txt` and `zhao_motif_analysis.html`
      - `fu_motif_analysis.txt` and `fu_motif_analysis.html`
  - Expression measurements from 3'UTR MPRA
      - `griesemer_expression_measurements.xlsx`
      - `siegel_expression_measurements_jurkat.csv`
      - `zhao_expression_measurements.xls`
      - `fu_expression_measurements_HEK293.txt` and `fu_expression_measurements_HELA.txt`
  - Pangolin prediction for 3'UTR MPRAs (not shown in manuscript)
      - `ptreseq_pangolin_prediction.txt`

### `PLOTS/` 
  - contain all figures saved as PDF


## Changelog
07/30/24
  - Initial commit: add analysis codes and data files to generate all figures for publication
    
02/18/25
  - Update codes following revision

05/20/25
  - Update codes following second revision

06/27/25
  - Update codes for final submission

07/03/25
  - Add scripts for performing splicing prediction

