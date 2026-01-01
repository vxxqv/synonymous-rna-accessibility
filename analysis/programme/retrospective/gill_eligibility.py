"""Outcome-blind, conditional sequence/coverage eligibility of an archived library.

This reads only XLSX columns A/E/F/G/I. It neither executes author code nor reads
RNA counts, effect estimates, errors, tests, or notebook outputs. No network use.
Sequence orientations remain conditional until the historical parser is verified.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
import re
import sys
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

BASES = "ACGT"
CODONS = ["".join(x) for x in product("TCAG", repeat=3)]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = dict(zip(CODONS, AA))
NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
ALLOWED = {"A", "E", "F", "G", "I"}
EXPECTED = {
    "A": "Sequence", "E": "CDS_8 Barcode Count", "F": "CDS_9 Barcode Count",
    "G": "CDS_10 Barcode Count", "I": "Combined 154-CDS-gDNA Barcode Counts",
}


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def manual_translate(dna):
    assert len(dna) % 3 == 0
    # Standard sense translation; a fixed GTG initiation codon is handled below.
    return "".join(CODE[dna[i:i+3]] for i in range(0, len(dna), 3))


def snapgene_raw(path):
    raw = Path(path).read_bytes()
    i, dna, features = 0, None, []
    while i < len(raw):
        assert i + 5 <= len(raw)
        packet_type = raw[i]
        size = int.from_bytes(raw[i+1:i+5], "big")
        payload = raw[i+5:i+5+size]
        assert len(payload) == size
        i += size + 5
        if packet_type == 0:
            assert dna is None
            dna = payload[1:].decode("ascii").upper()
        if packet_type == 10:
            features.extend(ET.fromstring(payload).iter("Feature"))
    assert i == len(raw) and dna is not None
    return dna, features


def construct(path):
    dna, raw_features = snapgene_raw(path)
    rec = SeqIO.read(path, "snapgene")
    assert dna == str(rec.seq).upper()
    cds_feature = next(f for f in rec.features if f.type == "CDS" and
                       any("aprE" in x for x in f.qualifiers.get("label", [])))
    raw = next(f for f in raw_features if f.attrib.get("type") == "CDS" and
               "aprE" in f.attrib.get("name", ""))
    ranges = [tuple(int(n) for n in s.attrib["range"].split("-")) for s in raw.findall("Segment")]
    strand = -1 if raw.attrib.get("directionality") == "2" else 1
    assert [(a-1, b) for a, b in ranges] == [(int(p.start), int(p.end)) for p in cds_feature.location.parts]
    assert strand == cds_feature.location.strand
    manual_cds = "".join(dna[a-1:b] for a, b in ranges)
    if strand == -1:
        manual_cds = manual_cds.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    cds = str(cds_feature.extract(rec.seq)).upper()
    assert manual_cds == cds
    prot = manual_translate(cds)
    assert prot == str(Seq(cds).translate(table=1))
    # Table 11 recognizes GTG at initiation. This does not affect fibre equality.
    biological_protein = str(Seq(cds).translate(table=11, cds=True))
    assert biological_protein == "M" + prot[1:-1]
    return rec, raw_features, cds, {
        "file": path.name, "raw_sha256": digest(path), "record_nt": len(dna),
        "cds_nt": len(cds), "cds_strand": strand,
        "cds_intervals_zero_based_half_open": [[a-1, b] for a, b in ranges],
        "cds_sha256": sha256(cds.encode()).hexdigest(),
        "protein_residues_excluding_stop": len(prot)-1,
        "sense_translation_sha256": sha256(prot.encode()).hexdigest(),
        "initiator_recognized_table11": True, "terminal_stop": prot[-1] == "*",
        "internal_stops": prot[:-1].count("*"), "translation_routes_agree": True,
    }


def permitted_rows(path):
    """Stream worksheet XML without accessing values of nonpermitted cells."""
    with ZipFile(path) as z:
        # Shared strings are decoded only for an allowed cell's index.
        shared = ET.fromstring(z.read("xl/sharedStrings.xml")) if "xl/sharedStrings.xml" in z.namelist() else None
        strings = list(shared) if shared is not None else []
        sheets = sorted(n for n in z.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n))
        assert len(sheets) == 1
        with z.open(sheets[0]) as stream:
            for _, row in ET.iterparse(stream, events=("end",)):
                if row.tag != "{" + NS["x"] + "}row":
                    continue
                values = {}
                for cell in row:
                    col = re.sub(r"\d", "", cell.attrib.get("r", ""))
                    if col not in ALLOWED:
                        continue
                    typ = cell.attrib.get("t")
                    if typ == "inlineStr":
                        value = "".join(cell.find("x:is", NS).itertext())
                    else:
                        node = cell.find("x:v", NS)
                        value = node.text if node is not None else None
                        if typ == "s" and value is not None:
                            value = "".join(strings[int(value)].itertext())
                    values[col] = value
                if int(row.attrib["r"]) == 1:
                    assert values == EXPECTED, values
                else:
                    yield values
                row.clear()


def components(vertices):
    remaining = set(vertices)
    sizes, edges = [], 0
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        stack, size = [seed], 0
        while stack:
            s = stack.pop()
            size += 1
            for i, original in enumerate(s):
                for b in BASES:
                    if b == original:
                        continue
                    nxt = s[:i] + b + s[i+1:]
                    if nxt in vertices:
                        edges += 1
                    if nxt in remaining:
                        remaining.remove(nxt)
                        stack.append(nxt)
        sizes.append(size)
    assert edges % 2 == 0
    return sorted(sizes), edges // 2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    src = args.source_dir
    rec, features, cds, current = construct(src / "pMG23.dna")
    _, _, native_cds, native = construct(src / "pMG15_integrated.dna")
    variable = next(f for f in features if f.attrib.get("name") == "CDS Random 8mer")
    a, b = (int(n) for n in variable.find("Segment").attrib["range"].split("-"))
    offset = (a-1) - current["cds_intervals_zero_based_half_open"][0][0]
    assert b-a+1 == 8 and current["cds_strand"] == 1 and offset == 467
    assert cds[:939] == native_cds[:939]
    row_errors, observed = Counter(), {}
    total = 0
    for row in permitted_rows(src / "cds_processed.xlsx"):
        total += 1
        s = row.get("A", "")
        if not isinstance(s, str) or len(s) != 8 or set(s) - set(BASES):
            row_errors["invalid_sequence"] += 1
            continue
        if s in observed:
            row_errors["duplicate_sequence"] += 1
            continue
        try:
            numbers = [int(row[x]) for x in ("E", "F", "G", "I")]
        except (ValueError, TypeError, KeyError):
            row_errors["invalid_permitted_DNA_count"] += 1
            continue
        if min(numbers) < 0 or sum(numbers[:3]) != numbers[3]:
            row_errors["DNA_count_sum_or_sign_mismatch"] += 1
            continue
        if numbers[3] <= 50:
            row_errors["not_strictly_above_DNA_threshold"] += 1
            continue
        observed[s] = numbers
    all_strings = ["".join(x) for x in product(BASES, repeat=8)]
    groups, key_of, stops = defaultdict(set), {}, set()
    for s in all_strings:
        altered = cds[:offset] + s + cds[offset+8:]
        manual = manual_translate(altered)
        assert manual == str(Seq(altered).translate(table=1))
        if "*" in manual[:-1]:
            stops.add(s)
            continue
        key = sha256(manual.encode()).hexdigest()
        groups[key].add(s)
        key_of[s] = key
    reference_key = sha256(manual_translate(cds).encode()).hexdigest()
    theoretical_histogram = Counter(len(v) for v in groups.values())
    output = {
        "status": "conditional_sequence_coverage_eligible_orientation_and_uncertainty_gates_open",
        "exposure": {"workbook_value_columns_read": sorted(ALLOWED),
                     "RNA_or_expression_outcome_values_read": False,
                     "external_author_code_executed": False, "network_requests_by_this_script": 0},
        "runtime": {"python": platform.python_version(), "biopython": Bio.__version__},
        "inputs": {p.name: {"bytes": p.stat().st_size, "sha256": digest(p)} for p in
                   [src / "pMG23.dna", src / "pMG15_integrated.dna", src / "cds_processed.xlsx",
                    src / "analysis_code_only.py", src / "audit_documentation_manifest.json"]},
        "audit_code_sha256": digest(Path(__file__)),
        "constructs": [current, native],
        "construct_relationship": {"identical_coding_prefix_nt": 939,
                                   "premature_stop_reporter_supported_by_published_methods": True,
                                   "library_variable_offset_zero_based": offset,
                                   "variable_codon_offsets": [2, 0, 1, 2, 0, 1, 2, 0],
                                   "affected_codon_numbers_one_based": [156, 157, 158, 159]},
        "coverage": {"rows": total, "accepted_unique_sequences": len(observed),
                     "missing_of_65536": 65536-len(observed), "exclusions": dict(row_errors),
                     "strict_DNA_gt50_verified": True,
                     "all_three_DNA_counts_positive": sum(all(x > 0 for x in n[:3]) for n in observed.values()),
                     "source_prethreshold_inner_join_on_all_six_libraries": True,
                     "missingness_not_attributable_to_DNA_threshold_alone": True},
        "theory": {"enumerated_strings": 65536, "translation_route_disagreements": 0,
                   "premature_stop_strings_excluded": len(stops),
                   "full_length_sense_strings": len(key_of), "same_protein_fibres": len(groups),
                   "fibre_size_histogram": dict(sorted(theoretical_histogram.items())),
                   "reference_protein_fibre_size": len(groups[reference_key])},
        "conditional_orientations": {},
        "unclosed_gates": ["historical source parser and amplicon files inaccessible HTTP403; no retry",
                           "workbook Sequence transcription orientation not independently authenticated",
                           "biological versus technical replication and library dependence not established",
                           "historical PyDESeq2 version and shrinkage uncertainty interpretation not pinned",
                           "marginal shrinkage uncertainty is not a simultaneous independent rectangular box"],
    }
    for orientation in ("as_written", "reverse_complement"):
        seen = set(observed) if orientation == "as_written" else {
            s.translate(str.maketrans("ACGT", "TGCA"))[::-1] for s in observed}
        summaries = []
        for key, vertices in sorted(groups.items()):
            present = vertices & seen
            expected_components, expected_edges = components(vertices)
            observed_components, observed_edges = components(present)
            summaries.append({"protein_sha256": key, "expected_vertices": len(vertices),
                              "observed_vertices": len(present), "complete": present == vertices,
                              "expected_component_sizes": expected_components,
                              "observed_component_sizes": observed_components,
                              "expected_edges": expected_edges, "observed_edges": observed_edges,
                              "reference_protein": key == reference_key})
        complete = [g for g in summaries if g["complete"]]
        output["conditional_orientations"][orientation] = {
            "complete_fibres": len(complete), "incomplete_fibres": len(groups)-len(complete),
            "nontrivial_complete_fibres": sum(g["expected_vertices"] > 1 for g in complete),
            "vertices_in_complete_fibres": sum(g["expected_vertices"] for g in complete),
            "edges_in_complete_fibres": sum(g["expected_edges"] for g in complete),
            "complete_disconnected_fibres": sum(len(g["expected_component_sizes"]) > 1 for g in complete),
            "observed_premature_stop_strings": len(seen & stops),
            "reference_protein": next(g for g in summaries if g["reference_protein"]),
            "fibre_summaries": summaries,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    compact = {k: v for k, v in output.items() if k in ("status", "coverage", "theory")}
    compact["orientations"] = {k: {x: y for x, y in v.items() if x != "fibre_summaries"}
                               for k, v in output["conditional_orientations"].items()}
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
