"""Recover the frozen E. coli primary inputs without changing analysis results.

The six full-transcript mirror inputs are checksum-pinned. CDS coordinates are
reconstructed from the bundled U00096.2 annotation; all annotated CDS strands
contribute to overlap exclusions. Output coordinates are 1-based inclusive.
"""
from __future__ import annotations

import argparse
from array import array
from collections import defaultdict
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import tempfile
import urllib.request

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "6f29547d70a2b72ceccf28218c790c45406946e8"
REFERENCE_SHA256 = "935e1c7464e8fc6315b4c777254b690ec4bec92b1cfdc55fbda928a611ea18ed"
CONDITIONS = ("incell", "cellfree", "kasugamycin")
CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_records():
    records = json.loads((ROOT / "data/provenance/source_manifest.json").read_text())
    needed = {f"{c}.{suffix}" for c in CONDITIONS for suffix in ("dbn", "react")}
    selected = {record["url"].rsplit("/", 1)[1]: record for record in records
                if record["url"].rsplit("/", 1)[1] in needed}
    require(set(selected) == needed, "Manifest does not identify all six primary source files")
    require(all(f"/{COMMIT}/" in r["url"] for r in selected.values()), "Unpinned source URL")
    return selected


def recover_sources(output, cached_raw=None):
    raw = output / "weeks194_orig"
    raw.mkdir(parents=True, exist_ok=True)
    verified = []
    for name, record in sorted(source_records().items()):
        destination = raw / name
        if cached_raw is not None:
            cached = cached_raw / name
            if not cached.exists():
                cached = cached_raw / "weeks194_orig" / name
            require(cached.is_file(), f"Cached source missing: {name}")
            require(digest(cached) == record["sha256"], f"Cached checksum mismatch: {name}")
            if cached.resolve() != destination.resolve():
                shutil.copyfile(cached, destination)
        elif not destination.exists() or digest(destination) != record["sha256"]:
            request = urllib.request.Request(record["url"], headers={"User-Agent": "native-input-reproduction"})
            content = urllib.request.urlopen(request, timeout=120).read()
            require(hashlib.sha256(content).hexdigest() == record["sha256"],
                    f"Downloaded checksum mismatch: {name}")
            destination.write_bytes(content)
        require(digest(destination) == record["sha256"], f"Source checksum mismatch: {name}")
        verified.append({"file": f"weeks194_orig/{name}", "url": record["url"],
                         "sha256": record["sha256"], "bytes": destination.stat().st_size})
    return verified


def read_dbn(path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    require(len(lines) % 3 == 0, f"Malformed DBN triples: {path.name}")
    records = {}
    for offset in range(0, len(lines), 3):
        header, sequence, structure = lines[offset:offset + 3]
        require(header.startswith(">") and header[1:] not in records, "Invalid or duplicate transcript header")
        require(set(sequence) <= set("ACGU") and len(sequence) == len(structure), "Invalid DBN sequence")
        stack, pairs = [], []
        for i, symbol in enumerate(structure, 1):
            if symbol == "(":
                stack.append(i)
            elif symbol == ")":
                require(bool(stack), "Unbalanced closing parenthesis")
                pairs.append((stack.pop(), i))
            else:
                require(symbol == ".", "Unsupported DBN symbol")
        require(not stack, "Unbalanced opening parenthesis")
        require(all(sequence[i - 1] + sequence[j - 1] in CANONICAL for i, j in pairs),
                "Noncanonical decoded DBN pair")
        records[header[1:]] = {"sequence": sequence, "pairs": pairs}
    return records


def read_profiles(path):
    result, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            require(name not in result, "Duplicate profile transcript")
            result[name] = []
        elif line.strip():
            require(name is not None, "Profile values before transcript header")
            index, value = line.split()
            require(int(index) == len(result[name]) + 1, "Nonconsecutive profile index")
            number = None if value == "NA" else float(value)
            require(number is None or math.isfinite(number), "Nonfinite measured profile value")
            result[name].append(number)
    return result


def exact_hits(genome, dna):
    hits = []
    for strand, sequence in (("+", dna), ("-", dna.translate(str.maketrans("ACGT", "TGCA"))[::-1])):
        start = genome.find(sequence)
        while start >= 0:
            hits.append((strand, start + 1, start + len(sequence)))
            start = genome.find(sequence, start + 1)
    return hits


def write_csv(path, rows, fields=None):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reconstruct(output, reference):
    require(digest(reference) == REFERENCE_SHA256, "Bundled U00096.2 GenBank checksum differs")
    genome_record = SeqIO.read(reference, "genbank")
    require(genome_record.id == "U00096.2" and len(genome_record) == 4639675, "Unexpected reference")
    genome = str(genome_record.seq)
    structures = {c: read_dbn(output / "weeks194_orig" / f"{c}.dbn") for c in CONDITIONS}
    profiles = {c: read_profiles(output / "weeks194_orig" / f"{c}.react") for c in CONDITIONS}
    baseline = structures["incell"]
    require(len(baseline) == 194, "Expected 194 transcript records")
    format_audit = []
    for condition in CONDITIONS:
        require(structures[condition].keys() == baseline.keys() == profiles[condition].keys(),
                "Condition/source transcript identifiers differ")
        for name, record in structures[condition].items():
            require(record["sequence"] == baseline[name]["sequence"], "Condition sequences differ")
            require(len(record["sequence"]) == len(profiles[condition][name]), "Profile length differs")
        format_audit.append({"condition": condition, "n_transcripts": 194,
                             "n_nucleotides": sum(len(r["sequence"]) for r in baseline.values()),
                             "n_pairs": sum(len(r["pairs"]) for r in structures[condition].values()),
                             "reactivity_missing": sum(v is None for vs in profiles[condition].values() for v in vs)})

    all_cds = [f for f in genome_record.features if f.type == "CDS"]
    simple = [f for f in all_cds if len(f.location.parts) == 1
              and isinstance(f.location.start, ExactPosition) and isinstance(f.location.end, ExactPosition)]
    require(len(all_cds) == 4321 and len(simple) == 4266, "Reference CDS feature count changed")
    # Combine positions by locus so repeated positions within joined features
    # cannot count as an overlap with the same CDS. Strand is deliberately not
    # used here: antisense coding obligations also exclude variable positions.
    positions_by_locus = defaultdict(set)
    for feature in all_cds:
        locus = feature.qualifiers["locus_tag"][0]
        positions_by_locus[locus].update(feature.location)
    multiplicity = array("H", [0]) * len(genome)
    for positions in positions_by_locus.values():
        for position in positions:
            multiplicity[position] += 1

    transcripts, coding = [], []
    for name, record in baseline.items():
        hits = exact_hits(genome, record["sequence"].replace("U", "T"))
        require(len(hits) == 1, f"Transcript lacks a unique exact full-length mapping: {name}")
        strand, start, end = hits[0]
        transcripts.append({"transcript_id": name, "length": len(record["sequence"]), "exact_hits": hits,
                            "chromosome": "U00096.2", "strand": strand, "start": start, "end": end,
                            "mapping_confidence": "unique_exact_full_length"})
        for feature in simple:
            fstart, fend = int(feature.location.start) + 1, int(feature.location.end)
            fstrand = "+" if feature.location.strand == 1 else "-"
            if fstrand != strand or fstart < start or fend > end:
                continue
            q = feature.qualifiers
            tstart = fstart - start + 1 if strand == "+" else end - fend + 1
            tend = fend - start + 1 if strand == "+" else end - fstart + 1
            cds = record["sequence"][tstart - 1:tend]
            native_cds = feature.extract(genome_record.seq)
            require(cds == str(native_cds).replace("T", "U"), "Transcript/CDS sequence mismatch")
            codon_start = int(q.get("codon_start", ["1"])[0])
            translated = "M" + str(native_cds[3:-3].translate(table=11))
            coding.append(dict(
                transcript_id=name, species="Escherichia coli", strain="K-12 MG1655", reference="U00096.2",
                annotation_date="2013-02-26", start=fstart, end=fend, strand=strand,
                gene=q.get("gene", [""])[0], locus_tag=q["locus_tag"][0], protein_id=q.get("protein_id", [""])[0],
                codon_start=codon_start, transl_table=q.get("transl_table", ["11"])[0],
                transcript_start=tstart, transcript_end=tend, complete_cds=True,
                mapping_method="unique exact full transcript sequence followed by same-strand fully contained GenBank CDS",
                coding_phase=(1 - codon_start) % 3, cds_length=len(cds), length_mod3=len(cds) % 3,
                start_codon=cds[:3], stop_codon=cds[-3:],
                overlap_positions_with_other_CDS=sum(multiplicity[p] > 1 for p in feature.location),
                protein_translation_verified=translated == q.get("translation", [""])[0]))
    require(len(coding) == 409, "Expected 409 complete same-strand CDS occurrences")
    require(all(r["length_mod3"] == 0 and r["stop_codon"] in {"UAA", "UAG", "UGA"} for r in coding),
            "Unexpected incomplete or non-stop-terminated mapped CDS")
    frozen = list(csv.DictReader((ROOT / "data/provenance/cds_mapping.csv").open()))
    rendered = [{key: str(value) for key, value in row.items()} for row in coding]
    require(rendered == frozen, "Reconstructed CDS ledger differs from frozen provenance")
    write_csv(output / "cds_mapping.csv", coding)
    transcript_by_id = {row["transcript_id"]: row for row in transcripts}
    fields = ["transcript_id", "transcript_pos", "genomic_pos", "reference", "gene", "locus_tag",
              "protein_id", "strand", "codon_index", "codon_phase", "overlaps_other_CDS"]
    position_count = overlap_count = 0
    with (output / "coding_position_map.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in coding:
            transcript = transcript_by_id[row["transcript_id"]]
            for tpos in range(row["transcript_start"], row["transcript_end"] + 1):
                gpos = transcript["start"] + tpos - 1 if row["strand"] == "+" else transcript["end"] - tpos + 1
                overlap = int(multiplicity[gpos - 1] > 1)
                writer.writerow(dict(transcript_id=row["transcript_id"], transcript_pos=tpos, genomic_pos=gpos,
                                     reference="U00096.2", gene=row["gene"], locus_tag=row["locus_tag"],
                                     protein_id=row["protein_id"], strand=row["strand"],
                                     codon_index=(tpos - row["transcript_start"]) // 3 + 1,
                                     codon_phase=(tpos - row["transcript_start"]) % 3 + 1,
                                     overlaps_other_CDS=overlap))
                position_count += 1
                overlap_count += overlap
    (output / "transcript_mapping.json").write_text(json.dumps(transcripts, indent=2))
    (output / "format_audit.json").write_text(json.dumps(format_audit, indent=2))
    return {"unique_exact_transcripts": len(transcripts), "CDS_occurrences": len(coding),
            "unique_loci": len({r["locus_tag"] for r in coding}), "CDS_position_rows": position_count,
            "overlap_position_rows": overlap_count, "mapped_CDS_with_overlap": sum(
                r["overlap_positions_with_other_CDS"] > 0 for r in coding),
            "nonstandard_translation_loci": sorted({r["locus_tag"] for r in coding
                                                    if not r["protein_translation_verified"]}),
            "frozen_CDS_ledger_identical": True, "format_audit": format_audit}


def verify_selection(output):
    """Redirect select() output into a temporary root; never touch frozen data."""
    spec = importlib.util.spec_from_file_location("native_selection_reproduction", ROOT / "analysis/native_models.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory(prefix="native-selection-verification-") as temporary:
        module.ROOT = Path(temporary)
        (module.ROOT / "data").mkdir()
        (module.ROOT / "results").mkdir()
        module.select(output)
        reproduced = json.loads((module.ROOT / "data/native_windows.json").read_text())
        frozen = json.loads((ROOT / "data/native_windows.json").read_text())
        require(reproduced == frozen, "Reproduced blocks differ from frozen blocks")
        audit = json.loads((module.ROOT / "results/native_selection.json").read_text())
        frozen_audit = json.loads((ROOT / "results/native_selection.json").read_text())
        require(audit == frozen_audit, "Reproduced selection/exclusion audit differs from frozen audit")
        require(len(reproduced) == 404, "Expected 404 frozen native blocks")
        return {"selected_blocks": len(reproduced), "unique_transcripts": len({r["transcript_id"] for r in reproduced}),
                "complete_block_objects_identical": True, "selection_audit_identical": True,
                "frozen_files_written": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Directory for recovered raw inputs and complete maps")
    parser.add_argument("--cached-raw", type=Path, help="Use verified existing six source files without network requests")
    parser.add_argument("--reference", type=Path, default=ROOT / "data/pairmap/U00096.2.gb")
    parser.add_argument("--verify-frozen", action="store_true", help="Reproduce the exact 404-block frozen selection in temporary files")
    args = parser.parse_args()
    output = args.output.resolve()
    require(not output.is_relative_to(ROOT / "data") and not output.is_relative_to(ROOT / "results"),
            "Choose a separate input-cache directory outside frozen data and results")
    output.mkdir(parents=True, exist_ok=True)
    sources = recover_sources(output, args.cached_raw)
    summary = reconstruct(output, args.reference)
    if args.verify_frozen:
        summary["selection_verification"] = verify_selection(output)
    summary["sources"] = sources
    summary["reference"] = {"accession": "U00096.2", "sha256": REFERENCE_SHA256}
    (output / "preparation_manifest.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({key: value for key, value in summary.items() if key not in {"sources", "format_audit"}}))


if __name__ == "__main__":
    main()
