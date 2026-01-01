BOUNDED SOURCE RECOVERY CLOSURE

FoldAtlas
---------
Both actual repository trees were inspected: mnori/foldatlas (original article link) and JIC-CSB/foldatlas (successor reached through the former data URL). This is stronger evidence than a homepage timeout.

Original README requires source files or a database copy. The repository does not embed a structure/probing database: static/downloads contains an ignore rule for foldatlas.sql.tar.gz; settings points CT models to external in_silico_structures/in_vivo_structures directories, reactivity counts to external files, and chromosome/GFF inputs to external source_data. Its importers parse these files; they do not regenerate the missing experimental data. The help page describes an SQL database containing sequences, CDS/UTR features, pair tables, BPPs, and replicate measurements.

The documented exact database URL http://www.foldatlas.com/static/downloads/foldatlas.sql.tar.gz now resolves to https://github.com/JIC-CSB/foldatlas and returns HTML, not gzip/SQL. The successor README still describes a live-dump import; its defaults point to that same URL, but the download/load tasks are commented out with a provisioning-problem note. No CT/FASTA/GFF/SQL data archive occurs in the inspected successor tree. Documented site_url http://foldatlas.jic.ac.uk fails DNS resolution in both HTTP/HTTPS in this session. The original HTTPS SQL URL times out. These are bounded endpoint observations, not a claim that no archival database copy exists anywhere.

One exact unconstrained example, AT5G65390.1/sid_207394, is embedded in original help source and recovered as foldatlas_embedded_example.dbn. It is not an experimentally guided alternative-conformer dataset and cannot substitute for the planned independent cohort.

Annotation is independently recoverable: the successor explicitly documents the UCL mirror. Exact URL:
http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/consolidated_annotation.Col_0.gff3.bz2
An8192-byte range returned206/application-x-bzip2 with BZh9 signature, total4887597bytes. Therefore do not describe historical CDS annotation as generally unavailable. The original Oxford host failed DNS; the documented UCL mirror works. Full GFF and genomic sequences have not been downloaded/mapped in this audit because the original structure/probing dataset remains the material gate.

Precise unresolved gate: obtain an authenticated exported SQL/CT/reactivity archive for the cited release, then verify transcript sequence identity and historical CDS coordinates. Code and one illustrative unconstrained model are available; an independent experimentally guided multi-transcript structure cohort was not recovered. No invented endpoint sweep, database installation, full reprocessing or new folding was attempted.

Cambray reporter source
-----------------------
Author OSF project a56vu is available. Its registered project license is BSD3-Clause, copyright GuillaumeCambray2018, obtained through node-linked license563c1cf88c5e4a3877f9e967. Preserve the license when redistributing covered material; the article's separate publisher terms are not replaced by the project license.

Supplementary_Data_15_Integrated_Phenotypic_Measurements.csv.zip is130122456bytes; https://osf.io/download/pm2qa/ supports byte ranges. A262144-byte compressed prefix was retrieved (206 Content-Range bytes0-262143/130122456) and its first member decompressed. This establishes actual accessible data, not only a dead inventory link. The first503complete records, all seed1, have96nt gs.sequence and103columns including seed/combi/rep; protein replicates/means; RNA decay/half-life; polysome; growth and sequence-design features.

Small nonrandom prefix pilot:503unique sequences,394distinct translated peptides,105raw Hamming1edges,29synonymous Hamming1edges with finite proteinmean at both endpoints,7two-step synonymous paths with Hamming2endpoints. Translation frame0 is supported by the author seed-generator's triplet code/stop rules and the author's RBS benchmark appending gs.sequence after an explicit upstream sequence endingATG. This pilot is not a global coverage estimate. Seed identity does not imply a fixed peptide. The observed small synonymous components have sizes1(452),2(16),3(5),4(1); these counts concern the restricted prefix and not the full library.

Author code15 supplies an explicit upstream context for its RBS benchmark, improving partial construct closure. It does not establish the complete measured downstream transcript or its isoforms/terminator. This is a precision gate, not grounds to call the reporter data unavailable.

What this dataset can answer: retrospective, benign reporter phenotype differences along observed synonymous edges; a later separately designed prediction/phenotypic-robustness question, with appropriate family-held-out validation and prior accessibility/initiation baselines. A whole-library graph audit could quantify how many intermediates are observed, because small measured edges already exist.

What it cannot directly answer: correctness of a compatibility theorem, actual retention of a specified base-pair target along a mutation route, or empirical thermodynamic-filter fragmentation. No experimentally resolved pair-target table was present in the inspected article-linked OSF inventory or downloaded author scripts. Predicted MFE/accessibility feature tables exist, so a computational structural target could be chosen later, but it would be a new model-based analysis and must not be called independent experimental structural validation. Protein output is not a pair-retention measurement. A full predictive expansion is therefore a different question from the present paper and should not be appended merely to increase ambition.

Evidence: repository_tree.txt, successor_repository_tree.txt, foldatlas_repository_evidence.json, successor_repository_evidence.json, recovery_manifest.json, recovery_final_manifest.json, successor_hosts_manifest.json, cambray_license.json, cambray_prefix_synonymous_pilot.json, cambray_design_manifest.json. No primary pipeline/manuscript files were changed.
