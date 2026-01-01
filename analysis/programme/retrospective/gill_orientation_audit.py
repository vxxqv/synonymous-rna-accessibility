"""Read-only historical parser/construct orientation audit. No author execution."""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path

from Bio import SeqIO


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    parser = args.source_dir / "barcodes_from_fastq_v16.py"
    expected = "88a84bfdeb15fa1fe00d6d268fb95052bd868299e6f09e90cdd373e6b2f3838c"
    assert digest(parser) == expected
    tree = ast.parse(parser.read_text())
    barcode_slices, anchors, unchanged_writes = [], [], 0
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "barcode_sequence" for t in node.targets):
            sub = node.value
            assert isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name) and sub.value.id == "seq"
            barcode_slices.append([ast.literal_eval(sub.slice.lower), ast.literal_eval(sub.slice.upper)])
        if isinstance(node, ast.Compare) and len(node.ops) == 1 and isinstance(node.ops[0], ast.Eq):
            sub = node.left
            if isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name) and sub.value.id == "seq":
                literal = ast.literal_eval(node.comparators[0])
                anchors.append((ast.literal_eval(sub.slice.lower), ast.literal_eval(sub.slice.upper), literal))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "write" and isinstance(node.func.value, ast.Name) and node.func.value.id == "mapped_outFile":
            value = node.args[0]
            assert isinstance(value, ast.BinOp) and isinstance(value.op, ast.Add)
            assert isinstance(value.left, ast.Name) and value.left.id == "barcode_sequence"
            assert ast.literal_eval(value.right) == "\n"
            unchanged_writes += 1
    assert barcode_slices == [[135, 143]] and unchanged_writes == 1
    assert [(a, b) for a, b, _ in anchors] == [(53, 68), (143, 153)]
    rec = SeqIO.read(args.source_dir / "pMG23.dna", "snapgene")
    seq = str(rec.seq).upper()
    n = len(seq)
    matches, anchor_hits = [], []
    for strand, oriented in [(1, seq), (-1, str(rec.seq.reverse_complement()).upper())]:
        circular = oriented + oriented
        for a, b, literal in anchors:
            positions = [i for i in range(n) if circular[i:i+len(literal)] == literal]
            anchor_hits.append({"parser_anchor_slice": [a, b], "record_strand": strand,
                                "anchor_positions_in_oriented_record": positions,
                                "implied_barcode_start_in_oriented_record": [(i+135-a) % n for i in positions]})
        for start in range(n):
            if all(circular[start+a:start+b] == literal for a, b, literal in anchors):
                lo, hi = start+135, start+143
                interval = [lo % n, hi % n] if strand == 1 else [(n-hi) % n, (n-lo) % n]
                matches.append({"read_strand_relative_to_record": strand,
                                "barcode_record_interval_zero_based_half_open": interval,
                                "read_window_start_in_oriented_record": start})
    resolved = len(matches) == 1 and matches[0]["read_strand_relative_to_record"] == 1 and matches[0]["barcode_record_interval_zero_based_half_open"] == [7065, 7073]
    out = {"status": "orientation_authenticated_as_written" if resolved else "orientation_unresolved_deposited_v16_parser_does_not_match_construct",
           "method": "Static AST inspection, then exact two-flank matching at specified relative offsets on both circular construct strands",
           "author_code_executed": False, "outcome_values_read": False,
           "parser_sha256": digest(parser), "construct_sha256": digest(args.source_dir / "pMG23.dna"),
           "audit_code_sha256": digest(Path(__file__)),
           "parser_barcode_slice": [135, 143], "parser_flank_slices": [[a, b] for a, b, _ in anchors],
           "unchanged_mapped_barcode_write_verified": True,
           "unique_matches_both_strands_and_all_offsets": matches,
           "individual_anchor_hits": anchor_hits,
           "interpretation": "Use as_written partition from eligibility_v1" if resolved else "The parser writes read-strand bases unchanged, but its exact flank locations do not authenticate the coding-site construct. Do not select an orientation from outcomes. No fallback offsets or mismatch tolerance were fitted.",
           "remaining_gates": ([] if resolved else ["historical CDS parser/amplicon identity and orientation"]) + ["replicate dependence", "historical statistical software version", "joint uncertainty calibration"]}
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
