"""Portable independent calculations copied unchanged from the pre-outcome reviewer.
create_data is synthetic fixture construction only. No author model imports.
"""
import hashlib
import itertools
import numpy as np
from Bio.Data import CodonTable

CODONS = sorted(set(map("".join, itertools.product("ACGT", repeat=3))) - {"TAA", "TAG", "TGA"})
AAS = list("ACDEFGHIKLMNPQRSTVWY")
PENALTIES = [0.1, 1.0, 10.0, 100.0, 1000.0]
MODELS = ["no_change", "baseline_linear", "baseline_core", "amino_acid_core", "codon_core", "association", "codon_noSE", "association_noSE", "association_nonlinear"]

def fold(group, which):
    return int(hashlib.sha256(f"cnot3-focus-v1|{which}|{group}".encode()).hexdigest(), 16) % 5


def standardize(a, b):
    mean, sd = np.mean(a, axis=0), np.std(a, axis=0, ddof=0)
    constant = np.ptp(a, axis=0) == 0
    sd[constant] = 1
    z, w = (a-mean)/sd, (b-mean)/sd
    z[:, constant] = 0
    w[:, constant] = 0
    return z, w


def transformed(train, assess, model):
    baseline = np.log([r["HL_WT"] for r in train])
    mean = np.mean(baseline)
    knots = np.quantile(baseline, [0.1, 0.5, 0.9], method="linear")
    medians = {p: np.median([r[p+"_GC"] for r in train if r[p+"_length_nt"] > 0]) for p in ("UTR5", "UTR3")}
    def base(rows):
        t = np.log([r["HL_WT"] for r in rows])
        columns = [t, (t-mean)**2, (t-mean)**3] + [np.maximum(t-k, 0) for k in knots]
        columns += [np.log([r["CDS_length_nt"] for r in rows]), np.array([r["CDS_GC"] for r in rows])]
        for p in ("UTR5", "UTR3"):
            columns += [np.log1p([r[p+"_length_nt"] for r in rows]),
                        np.array([r[p+"_GC"] if r[p+"_length_nt"] else medians[p] for r in rows])]
        columns += [np.array([float(r[p+"_length_nt"] == 0) for r in rows]) for p in ("UTR5", "UTR3")]
        if model not in ("codon_noSE", "association_noSE"):
            columns += [np.log([r["SE_IP"] for r in rows])]
        if model == "amino_acid_core":
            columns += [np.array([r["aa_"+aa] for r in rows]) for aa in AAS]
        if model in ("codon_core", "association", "codon_noSE", "association_noSE", "association_nonlinear"):
            columns += [np.array([r["codon_"+c] for r in rows]) for c in CODONS]
        return np.array(columns).T
    a, b = base(train), base(assess)
    if model.startswith("association"):
        z, w = standardize(a, b)
        z, w = np.c_[np.ones(len(z)), z], np.c_[np.ones(len(w)), w]
        raw_a = np.array([r["LOG2_FC_IP"] for r in train])
        raw_b = np.array([r["LOG2_FC_IP"] for r in assess])
        ca, cb = [raw_a], [raw_b]
        if model == "association_nonlinear":
            for knot in np.quantile(raw_a, [0.1, 0.5, 0.9], method="linear"):
                ca.append(np.maximum(raw_a-knot, 0)); cb.append(np.maximum(raw_b-knot, 0))
        for ra, rb in zip(ca, cb):
            sx, tx = standardize(ra[:, None], rb[:, None])
            beta = np.linalg.pinv(z, rcond=1e-12) @ sx[:, 0]
            residual, assessment_residual = sx[:, 0]-z@beta, tx[:, 0]-w@beta
            if np.linalg.norm(residual) <= 1e-10*np.sqrt(len(train)):
                residual, assessment_residual = np.zeros(len(train)), np.zeros(len(assess))
            else:
                rmean, rsd = residual.mean(), residual.std(ddof=0)
                residual, assessment_residual = (residual-rmean)/rsd, (assessment_residual-rmean)/rsd
            a, b = np.c_[a, residual], np.c_[b, assessment_residual]
    return a, b


def gram_prediction(a, y, b, lam):
    z, w = standardize(a, b)
    beta = np.linalg.solve(z.T@z + lam*np.eye(z.shape[1]), z.T@(y-y.mean()))
    return y.mean()+w@beta


def create_data():
    rng = np.random.default_rng(912731)
    translation = CodonTable.unambiguous_dna_by_id[1].forward_table
    rows, targets = [], {}
    for g in range(160):
        for j in range(2):
            group, gene = f"SYN_GROUP_{g:04d}", f"SYN_{g:04d}_{j}"
            n = int(rng.integers(100, 450)); counts = rng.multinomial(n, rng.dirichlet(np.ones(61)*3))
            fraction = counts/n
            aa_fraction = {a: sum(f for c, f in zip(CODONS, fraction) if translation[c] == a) for a in AAS}
            t, ip, se = rng.normal(2.3, 0.55), rng.normal(), np.exp(rng.normal(-1.5, 0.4))
            r = {"GeneID": gene, "group_ID": group, "HL_WT": float(np.exp(t)), "LOG2_FC_IP": float(ip), "SE_IP": float(se),
                 "CDS_length_nt": (n+1)*3, "CDS_GC": float(sum(counts[i]*sum(ch in "GC" for ch in c) for i,c in enumerate(CODONS))/((n+1)*3))}
            for p in ("UTR5", "UTR3"):
                length = 0 if rng.random() < 0.16 else int(rng.integers(15, 700))
                r[p+"_length_nt"] = length; r[p+"_GC"] = None if not length else float(rng.uniform(0.25, 0.75))
            r.update({"codon_"+c: float(v) for c,v in zip(CODONS, fraction)})
            r.update({"aa_"+a: float(v) for a,v in aa_fraction.items()})
            targets[gene] = float(0.7*t+0.9*ip+3*fraction[0]+rng.normal(0,0.07))
            rows.append(r)
    train = [r for r in rows if fold(r["group_ID"], "outer") != 0]
    test = [r for r in rows if fold(r["group_ID"], "outer") == 0]
    return train, test, targets


