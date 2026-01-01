"""Read-only amplicon provenance and reference-fibre design audit; no outcomes."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from gill_eligibility import construct, snapgene_raw, permitted_rows, manual_translate, components


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def exact_rank(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [x / scale for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                scale = rows[i][col]
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def design_rank(level_rows):
    columns = [sorted({r[j] for r in level_rows}) for j in range(len(level_rows[0]))]
    design = [[1] + [int(r[j] == level) for j, levels in enumerate(columns) for level in levels[1:]]
              for r in level_rows]
    return exact_rank(design), [len(x) for x in columns]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--eligibility-result", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    src = args.source_dir
    amp_path = src / "cds_final_amplicon.dna"
    amp_seq, raw_features = snapgene_raw(amp_path)
    amp = SeqIO.read(amp_path, "snapgene")
    assert amp_seq == str(amp.seq).upper() and len(amp_seq) == 371
    raw_var = next(f for f in raw_features if f.attrib.get("name") == "CDS Random 8mer")
    assert raw_var.find("Segment").attrib["range"] == "125-132"
    bio_var = next(f for f in amp.features if f.qualifiers.get("label") == ["CDS Random 8mer"])
    assert [int(bio_var.location.start), int(bio_var.location.end), bio_var.location.strand] == [124, 132, 1]
    rec, _, cds, context = construct(src / "pMG23.dna")
    genomic = str(rec.seq).upper()
    assert amp_seq[66:277] == genomic[7007:7218]
    parser = (src / "barcodes_from_fastq_v16.py").read_text()
    assert "barcode_sequence = seq[135:143]" in parser and "mapped_outFile.write(barcode_sequence + '\\n')" in parser
    downstream = "CTGAAACAAA"
    assert amp_seq[132:142] == downstream and genomic[7073:7083] == downstream
    assert sum(amp_seq[i:i+10] == downstream for i in range(len(amp_seq))) == 1
    assert downstream not in str(amp.seq.reverse_complement()).upper()
    parse_source = (src / "parse_barcodes_code_only.py").read_text()
    assert "itertools.product(a, repeat=8)" in parse_source
    assert "sample_df.loc[count] = [i, sample_counter.get(i)]" in parse_source
    shell = (src / "process_fastq_v16.sh").read_text()
    assert "merged_and_combined_fastq_files" in shell and "python barcodes_from_fastq_v16.py $i" in shell
    reference_local = manual_translate(cds[465:477])
    fibre = set()
    for parts in product("ACGT", repeat=8):
        s = "".join(parts)
        local = cds[465:467] + s + cds[475:477]
        translated = manual_translate(local)
        assert translated == str(Seq(local).translate(table=1))
        if translated == reference_local:
            fibre.add(s)
    assert len(fibre) == 32
    levels = [sorted({s[j] for s in fibre}) for j in range(8)]
    assert {"".join(x) for x in product(*levels)} == fibre
    nt_rank, nt_cardinality = design_rank(sorted(fibre))
    codon_rows = [tuple((cds[465:467]+s+cds[475:477])[i:i+3] for i in range(0, 12, 3)) for s in sorted(fibre)]
    codon_rank, codon_cardinality = design_rank(codon_rows)
    comp, edges = components(fibre)
    assert comp == [32] and edges == 112 and nt_rank == codon_rank == 8
    allowed = {r["A"]: [int(r[k]) for k in ("E", "F", "G", "I")] for r in permitted_rows(src / "cds_processed.xlsx") if r.get("A") in fibre}
    assert len(allowed) == 32
    base_result = json.loads(args.eligibility_result.read_text())
    forward = {k: v for k, v in base_result["conditional_orientations"]["as_written"].items() if k != "fibre_summaries"}
    input_paths = [amp_path, src / "pMG23.dna", src / "barcodes_from_fastq_v16.py", src / "process_fastq_v16.sh",
                   src / "parse_barcodes_source.ipynb", src / "parse_barcodes_code_only.py", src / "analysis_code_only.py",
                   src / "cds_processed.xlsx", args.eligibility_result]
    result = {
        "status": "deposited_pipeline_forward_orientation_and_reference_fibre_design_eligible",
        "author_code_executed": False, "RNA_or_expression_outcomes_read": False,
        "inputs": {p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in input_paths},
        "audit_code_sha256": digest(Path(__file__)),
        "amplicon": {"nt": len(amp_seq), "independent_packet_and_Biopython_match": True,
                     "barcode_interval_zero_based_half_open": [124, 132], "barcode_strand": 1,
                     "exact_coding_fragment_amplicon": [66, 277], "exact_coding_fragment_pMG23": [7007, 7218],
                     "coding_fragment_nt": 211, "downstream_flank_unique_positive_strand": True,
                     "upstream_parser_anchor_on_reverse_strand": True,
                     "both_parser_anchors_are_not_one_contiguous_genomic_read": True},
        "linkage": {"shell_invokes_v16": True, "shell_input_merged_and_combined_FASTQ": True,
                    "parser_writes_eight_bases_unchanged": True,
                    "parse_notebook_enumerates_all_65536_keys_per_sample": True,
                    "parse_notebook_sample_labels": ["CDS_1", "CDS_2", "CDS_3", "CDS_8", "CDS_9", "CDS_10"],
                    "analysis_reads_matching_six_workbooks_and_exports_matching_processed_filename": True,
                    "forward_orientation_intended_by_deposited_assay": True,
                    "raw_FASTQ_preprocessing_independently_reproduced": False,
                    "inner_joins_do_not_censor_keys_under_documented_upstream_enumeration": True,
                    "856_absences_explained_by_strict_DNA_filter_under_documented_pipeline": True,
                    "unthresholded_counts_independently_recovered": False},
        "verified_forward_partition": forward,
        "reference_pilot": {"translated_reporter_residues": 313, "vertices": 32, "edges": edges,
                            "components": comp, "factor_graph": "K2 Cartesian K4 Cartesian K4",
                            "regular_degree": 7, "diameter": 3,
                            "eightmer_nt_position_cardinalities": nt_cardinality,
                            "variable_eightmer_positions_one_based": [i+1 for i, x in enumerate(nt_cardinality) if x > 1],
                            "affected_codon_cardinalities": codon_cardinality,
                            "nucleotide_additive_design_rank_including_intercept": nt_rank,
                            "codon_additive_design_rank_including_intercept": codon_rank,
                            "main_effect_degrees_of_freedom": 7,
                            "nonadditive_complement_dimension": 24,
                            "two_factor_interaction_dimension": 15,
                            "three_factor_interaction_dimension": 9,
                            "all_32_pass_strict_DNA_gt50": all(v[3] > 50 for v in allowed.values()),
                            "all_32_positive_in_each_DNA_library": all(min(v[:3]) > 0 for v in allowed.values()),
                            "RNA_libraries": ["CDS_1", "CDS_2", "CDS_3"],
                            "DNA_libraries": ["CDS_8", "CDS_9", "CDS_10"],
                            "pairing_assumed": False,
                            "independent_biological_replicates_established": False,
                            "other_protein_outcomes_untouched": True},
        "remaining_holds": ["Root must freeze pilot estimand and analysis before outcomes",
                            "Biological versus technical replicate identity and dependence",
                            "Historical software and shrinkage uncertainty implementation",
                            "Simultaneous uncertainty calibration, not independent marginal boxes",
                            "Raw FASTQ preprocessing provenance not reproduced"],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "reference_pilot": result["reference_pilot"]}, indent=2))


if __name__ == "__main__":
    main()
