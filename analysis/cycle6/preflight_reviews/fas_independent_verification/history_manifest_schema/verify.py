"""Independent exact audit of a published FAS design, without measured counts.

Uses a literal standard-code table, a restricted GenBank parser, codon-level
factor enumeration, and all-pairs nucleotide adjacency. No author module import.
"""
from pathlib import Path
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import platform
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
SOURCE = BASE / "work/research_cycle6/measured_synonymy_sources"
NATIVE_SHA = "e8964c4db7d19acd5f6a1075106ec67087ad4bc5a8dc3a986f7eccb0a7f31ca7"
CODE = dict(zip(("".join(x) for x in itertools.product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_pattern(pattern):
    domains = []
    i = 0
    while i < len(pattern):
        if pattern[i] == "[":
            stop = pattern.index("]", i)
            token = pattern[i + 1:stop]
            i = stop + 1
        else:
            token = pattern[i]
            i += 1
        require(token and set(token) <= set("ACGT") and len(set(token)) == len(token),
                "Unsupported pattern token")
        domains.append(tuple(token))
    return domains


def translate(cds):
    require(len(cds) % 3 == 0, "Incomplete CDS")
    return "".join(CODE[cds[i:i+3]] for i in range(0, len(cds), 3))


def hamming(a, b):
    require(len(a) == len(b), "Unequal sequence lengths")
    return sum(x != y for x, y in zip(a, b))


def main():
    result_path = HERE / "result.json"
    require(not result_path.exists(), "Refusing overwrite")
    started = datetime.now(timezone.utc).isoformat()
    timer = time.perf_counter()
    inputs = [SOURCE / x for x in ["native_refseq.gb", "source_author_pin_1.txt",
              "source_author_source_2.txt", "sequence_contract_check.json",
              "declared_human_reference_class_sequences.json", "scaling_law_tree.json",
              "author_pin_1.json", "author_pin_0.json"]]
    before = {str(p.relative_to(BASE)): sha(p) for p in inputs}
    require(sha(inputs[0]) == NATIVE_SHA, "Native reference hash drift")
    native_text = inputs[0].read_text(encoding="utf-8")
    sequence = re.sub("[^acgt]", "", native_text.split("ORIGIN", 1)[1].split("//", 1)[0]).upper()
    version = re.search(r"^VERSION\s+(\S+)", native_text, re.M).group(1)
    cds_matches = list(re.finditer(r"^     CDS             (\d+)\.\.(\d+)\n", native_text, re.M))
    require(len(cds_matches) == 1, "Need one explicit linear forward CDS")
    match = cds_matches[0]
    cds_start, cds_end = int(match[1])-1, int(match[2])
    next_feature = re.search(r"^     \S", native_text[match.end():], re.M)
    cds_block = native_text[match.start():match.end()+next_feature.start()]
    require('/codon_start=1' in cds_block and '/transl_except=' not in cds_block,
            "Unsupported coding convention")
    transl_table = re.search(r"/transl_table=(\d+)", cds_block)
    require(transl_table is None or transl_table[1] == "1", "Not standard nuclear code")
    protein_id = re.search(r'/protein_id="([^"]+)"', cds_block)[1]
    annotated = re.sub(r"\s+", "", re.search(r'/translation="([^"]+)"', cds_block)[1])
    cds = sequence[cds_start:cds_end]
    protein = translate(cds)
    require(protein[-1] == "*" and protein[:-1] == annotated, "Reference translation mismatch")

    shell = inputs[1].read_text(encoding="utf-8")
    patterns = re.findall(r"^\s*grep\s+(\S+)\s+", shell, re.M)
    require(len(patterns) == 2 and patterns[0] == patterns[1], "Different input/output designs")
    pattern = patterns[0]
    domains = parse_pattern(pattern)
    exon_features = [(int(a)-1, int(b)) for a,b in
                     re.findall(r"^     exon            (\d+)\.\.(\d+)$", native_text, re.M)]
    candidates = [(a,b) for a,b in exon_features if re.fullmatch(pattern, sequence[a:b])]
    require(len(candidates) == 1, "Design does not uniquely select a native annotated exon")
    exon_start, exon_end = candidates[0]
    wt = sequence[exon_start:exon_end]
    offset = exon_start - cds_start
    require(0 <= offset and offset + len(wt) <= len(cds), "Exon outside CDS")
    require(len(wt) == len(domains), "Pattern length mismatch")
    variable = [i for i,d in enumerate(domains) if len(d) > 1]
    codon_indices = sorted({(offset+i)//3 for i in variable})
    factors = []
    for index in codon_indices:
        local_domains = []
        for j in range(3):
            exon_index = index*3+j-offset
            local_domains.append(domains[exon_index] if 0 <= exon_index < len(wt)
                                 else (cds[index*3+j],))
        codons = sorted("".join(c) for c in itertools.product(*local_domains))
        amino_counts = Counter(CODE[c] for c in codons)
        native_codon = cds[index*3:index*3+3]
        keep = [c for c in codons if CODE[c] == CODE[native_codon]]
        local_edges = [[i,j] for i,j in itertools.combinations(range(len(keep)), 2)
                       if hamming(keep[i], keep[j]) == 1]
        factors.append({"cds_codon_1based": index+1, "native_codon": native_codon,
                        "amino_acid": CODE[native_codon], "all_codons": codons,
                        "amino_multiplicities": dict(sorted(amino_counts.items())),
                        "human_codons": keep, "local_edges": local_edges,
                        "variable_exon_positions_1based":
                            [i+1 for i in variable if (offset+i)//3 == index]})

    # Orthogonal class-size calculation by local amino-acid multiplicities.
    class_sizes = [1]
    for factor in factors:
        class_sizes = [a*b for a,b in itertools.product(
            class_sizes, factor["amino_multiplicities"].values())]
    expected_design_count = math.prod(len(d) for d in domains)
    require(sum(class_sizes) == expected_design_count, "Factor count does not cover whole design")

    groups = defaultdict(list)
    for choices in itertools.product(*(f["all_codons"] for f in factors)):
        chars = list(cds)
        for index, codon in zip(codon_indices, choices):
            chars[index*3:index*3+3] = codon
        full = "".join(chars)
        exon = full[offset:offset+len(wt)]
        require(re.fullmatch(pattern, exon) is not None, "Codon product escapes design")
        groups[translate(full)].append(exon)
    visited = [s for g in groups.values() for s in g]
    require(len(visited) == len(set(visited)) == expected_design_count, "Incomplete or duplicated enumeration")
    require(Counter(map(len, groups.values())) == Counter(class_sizes), "Class multiplicity disagreement")
    human = sorted(groups[protein])
    edges = [[i,j] for i,j in itertools.combinations(range(len(human)), 2)
             if hamming(human[i], human[j]) == 1]
    coordinates = []
    for exon in human:
        full = cds[:offset] + exon + cds[offset+len(wt):]
        coordinates.append(tuple(f["human_codons"].index(full[index*3:index*3+3])
                                 for index, f in zip(codon_indices, factors)))
    require(len(set(coordinates)) == len(human), "Noninjective factor mapping")
    require(set(coordinates) == set(itertools.product(*(range(len(f["human_codons"])) for f in factors))),
            "Fibre does not fill product")
    graph_checks = 0
    edge_set = {tuple(e) for e in edges}
    for i,j in itertools.combinations(range(len(human)), 2):
        differing = [k for k,(a,b) in enumerate(zip(coordinates[i], coordinates[j])) if a != b]
        product_edge = False
        if len(differing) == 1:
            k = differing[0]
            local_pair = sorted([coordinates[i][k], coordinates[j][k]])
            product_edge = local_pair in factors[k]["local_edges"]
        require(product_edge == ((i,j) in edge_set), "Product adjacency mismatch")
        graph_checks += 1
    active = [f for f in factors if len(f["human_codons"]) > 1]
    require(sorted((len(f["human_codons"]),len(f["local_edges"])) for f in active)
            == [(2,1)]*3 + [(3,2)]*2, "Wrong nontrivial factors")
    for f in active:
        if len(f["human_codons"]) == 3:
            degrees = Counter(i for e in f["local_edges"] for i in e)
            require(sorted(degrees.values()) == [1,1,2], "Three-vertex factor is not P3")
    adj = {i:set() for i in range(len(human))}
    for i,j in edges: adj[i].add(j);adj[j].add(i)
    distance_checks = 0
    diameter = 0
    for start in adj:
        dist = {start:0}; queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in adj[i]:
                if j not in dist:dist[j]=dist[i]+1;queue.append(j)
        require(len(dist) == len(human), "Disconnected human fibre")
        for j,d in dist.items():
            require(d == hamming(human[start],human[j]), "Nongeodesic fibre")
            distance_checks += 1
            diameter=max(diameter,d)
    ancestor = re.search(r"^WTSEQ=([ACGT]+)$", inputs[2].read_text(encoding="utf-8"), re.M)[1]
    require(re.fullmatch(pattern, ancestor) is not None, "Declared ancestor outside pattern")
    ancestor_protein = translate(cds[:offset]+ancestor+cds[offset+len(wt):])
    observed = {"declared_variants": expected_design_count, "protein_classes": len(groups),
                "protein_class_size_histogram": {str(k):v for k,v in sorted(Counter(class_sizes).items())},
                "human_reference_class_vertices":len(human), "human_reference_class_edges":len(edges),
                "native_exon_start0":exon_start,"native_cds_start0":cds_start,
                "native_exon_start_cds_offset0":offset,"native_exon_start_codon_offset0":offset%3,
                "native_accession":version,"native_protein_id":protein_id,
                "ancestor_same_protein_as_human":ancestor_protein==protein,
                "human_class_hamming_distance_from_native_histogram":
                    {str(k):v for k,v in sorted(Counter(hamming(s,wt) for s in human).items())}}
    author = json.loads(inputs[3].read_text(encoding="utf-8"))
    for key,value in observed.items(): require(author[key] == value, f"Reported claim mismatch: {key}")
    declared = json.loads(inputs[4].read_text(encoding="utf-8"))
    require(declared["sequences"] == human, "Human sequence list differs")
    require(expected_design_count == 3072 and len(groups)==80 and len(human)==72 and len(edges)==204,
            "Headline count mismatch")
    result = {"status":"PASS_FINITE_SEQUENCE_AND_GRAPH_CLAIMS", "observed":observed,
              "native_exon_1based_in_transcript":[exon_start+1,exon_end],
              "native_CDS_1based_in_transcript":[cds_start+1,cds_end],
              "exon_length":len(wt),"variable_positions":len(variable),
              "codon_factorization":factors,"nontrivial_graph_factors":"P3 x P3 x K2 x K2 x K2",
              "diameter":diameter,"degree_histogram":dict(sorted(Counter(map(len,adj.values())).items())),
              "all_pairs_adjacency_checks":graph_checks,"ordered_distance_checks":distance_checks,
              "protein_includes_native_termination":True,"premature_stop_protein_classes":sum('*' in p[:-1] for p in groups),
              "source_count_files_opened":False,"measured_coverage_authenticated":False,
              "minigene_protein_translation_authenticated":False,
              "exposure_incident":"Five illustrative README occurrence values were displayed during source review; not used by this computation. Separate ledger retains identifiers.",
              "checks_against_report":list(observed),"independent_algorithm":"Codon-block product and literal code; direct all-pairs Hamming graph; no author or Biopython import"}
    result_path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    after = {str(p.relative_to(BASE)):sha(p) for p in inputs}
    require(before == after, "Inputs changed during computation")
    repo = BASE / "outputs/synonymous-rna-accessibility"
    commit = subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],text=True).strip()
    dirty = bool(subprocess.check_output(["git","-C",str(repo),"status","--porcelain"],text=True).strip())
    manifest = {"schema_version":2,"claim_id":"FAS_declared_native_reference_fibre_independent_audit",
                "repository":{"commit":commit,"dirty":dirty},"command":[sys.executable,str(Path(__file__).resolve())],
                "environment":{"software":[{"name":"CPython","version":platform.python_version()}],"hardware":platform.machine()},
                "mathematics":{"assertion_tested":"Exact declared design enumeration, protein classes and native-reference fibre Cartesian graph",
                               "coefficient_domain":"Exact finite strings and integer counts","conventions":"Standard nuclear code; native complete included-isoform CDS; nucleotide Hamming-one edges",
                               "inputs":list(before),"bounds":{"designs":3072,"human_vertices":72},
                               "non_claims":["Measured coverage","Splicing outcomes","Minigene translation or protein output","Independent biological replication"]},
                "randomness":{"used":False,"generator":"","seed":None},
                "run":{"started_at":started,"runtime_seconds":time.perf_counter()-timer,"exit_status":0,"status":"completed",
                       "resource_limits":{"memory_bytes":None,"cpu_seconds":None,"max_cores":1,"max_threads":1}},
                "outputs":[{"path":str(result_path.relative_to(BASE)),"sha256":sha(result_path),"kind":"result"}],
                "checks":["Source stable before/after","Full native annotation translation","Codon-product class counts","2556 direct adjacency equivalences","5184 exact distances","All reported human vertices identical"],
                "result":"Implementation and finite assertion verified in stated design",
                "residual_risks":["Experimental minigene identity and measured coverage not authenticated", "Five illustrative README counts previously displayed, excluded from computation"],
                "input_artifacts":[{"path":p,"sha256":h,"sha256_after":after[p]} for p,h in before.items()],
                "declared_results":[{"path":str(result_path.relative_to(BASE)),"description":"Independent exact finite graph verification"}],
                "verifier_sha256":sha(Path(__file__).resolve())}
    (HERE / "manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"observed":observed,"factor_details":active,
                      "result_sha256":sha(result_path),"runtime_seconds":manifest["run"]["runtime_seconds"]}))


if __name__ == "__main__":
    main()
