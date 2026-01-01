"""Verify cached Git metadata and log the source README exposure incident."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "measured_synonymy_sources"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    output = HERE / "provenance_and_exposure.json"
    if output.exists():
        raise ValueError("Refusing overwrite")
    inventory = json.loads((SOURCE / "scaling_law_tree.json").read_text())
    commit = json.loads(json.loads((SOURCE / "author_pin_0.json").read_text())["structuredContent"]["content"])
    body = json.loads((SOURCE / "author_pin_1.json").read_text())["structuredContent"]["content"]
    raw = body.encode()
    blob = hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()
    entry = [x for x in inventory["tree"] if x["path"].endswith(
        "001_Phylogeny_library_in_HEK293/005_remove_sequencing_errors.sh")]
    if len(entry) != 1 or entry[0]["sha"] != blob:
        raise ValueError("Source Git blob mismatch")
    if (SOURCE / "source_author_pin_1.txt").read_text() != body:
        raise ValueError("Saved script content differs from source response")
    entries = inventory["tree"]
    trees = {"":commit["commit"]["tree"]["sha"]}
    trees.update({x["path"]:x["sha"] for x in entries if x["type"] == "tree"})
    for directory, expected in trees.items():
        prefix = directory+"/" if directory else ""
        direct = [x for x in entries if x["path"].startswith(prefix)
                  and "/" not in x["path"][len(prefix):]]
        direct.sort(key=lambda x:x["path"][len(prefix):]+("/" if x["type"]=="tree" else ""))
        binary = b"".join((x["mode"].lstrip("0")+" "+x["path"][len(prefix):]+"\0").encode()
                          +bytes.fromhex(x["sha"]) for x in direct)
        actual = hashlib.sha1(f"tree {len(binary)}\0".encode()+binary).hexdigest()
        if actual != expected:
            raise ValueError(f"Recursive Git tree mismatch: {directory}")
    known = set(json.loads((SOURCE / "declared_human_reference_class_sequences.json").read_text())["sequences"])
    exposed = []
    readme = SOURCE / "source_author_source_0.txt"
    for number, line in enumerate(readme.read_text().splitlines(),1):
        match = re.match(r"^\|\s+([ACGT]{63})\s+\|", line)
        if match:
            exposed.append({"source_line":number,"sequence":match[1],
                            "native_reference_class_member":match[1] in known,
                            "exposure":"Occurrence example displayed in source README; no assay effect inferred"})
    if len(exposed) != 5:
        raise ValueError("Unexpected source exposure inventory")
    result = {"status":"PASS_CACHED_GIT_OBJECT_CHAIN",
              "source_commit":commit["sha"],"canonical_tree":commit["commit"]["tree"]["sha"],
              "verified_tree_objects":len(trees),"script_git_blob_sha1":blob,
              "source_script_sha256":sha(SOURCE / "source_author_pin_1.txt"),
              "source_url":json.loads((SOURCE / "author_pin_1.json").read_text())["structuredContent"]["url"],
              "tree_inventory_root_sha_semantics":"Inventory top sha is the requested commit; independently reconstructed canonical tree matches commit object",
              "source_README_sha256":sha(readme),"exposure_incident_records":exposed,
              "native_reference_class_exposed_ids":sum(x["native_reference_class_member"] for x in exposed),
              "count_files_opened":False,"counts_used_in_sequence_calculation":False,
              "numeric_occurrence_values_reextracted":False,
              "qualification":"Cached source provenance only, no fresh network authentication or measured genotype coverage"}
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"trees":len(trees),
                      "incident_native_class_members":result["native_reference_class_exposed_ids"],
                      "sha256":sha(output)}))


if __name__ == "__main__":
    main()
