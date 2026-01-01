"""Portable independent verifier. No author model imports or XLSX reader."""
import sys
if not __debug__:
    raise RuntimeError("Optimized Python is refused before input access")
sys.dont_write_bytecode = True

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import re
import statistics

import Bio
import numpy as np

HERE = Path(__file__).resolve().parent
PROTOCOL_SHA = "219dda125d7a8167b2a45fa67499ac4dc24118e6edcbdee83b540f8b77e27f22"
AUTHOR_SHA = "be75c90e4076f29952553fa2a22529c2640db0691cc8237e5078c1d1e2cc25d7"
CORE_SHA = "a5255167a8182d36e107db804d42b87a488e3dc22361cc75c57772a7729c8786"
PRED_ATOL = 1e-9
PRED_RTOL = 1e-10
SCORE_ATOL = 1e-12
SCORE_RTOL = 1e-10
VERSION = "cnot3-operational-exploratory-v2"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    def unique(pairs):
        d = {}
        for key,value in pairs:
            if key in d:
                raise ValueError("Duplicate JSON key: "+key)
            d[key] = value
        return d
    def bad(value):
        raise ValueError("Nonfinite JSON constant: "+value)
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique, parse_constant=bad)


class Checks:
    def __init__(self):
        self.count = 0
        self.values = 0
        self.maximum_abs = 0.0
        self.failures = []

    def exact(self, name, actual, expected):
        self.count += 1
        if actual != expected:
            self.failures.append({"check":name,"kind":"exact_mismatch","actual":actual,"expected":expected})

    def number(self, name, actual, expected, prediction=False):
        self.count += 1
        if expected is None or isinstance(expected, bool):
            if actual is not expected:
                self.failures.append({"check":name,"kind":"null_or_boolean_mismatch","actual":actual,"expected":expected})
            return
        try:
            a,b = np.asarray(actual,dtype=float),np.asarray(expected,dtype=float)
            if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
                raise ValueError("Shape/nonfinite mismatch")
            d = np.abs(a-b)
            self.values += a.size
            self.maximum_abs = max(self.maximum_abs,float(np.max(d)) if d.size else 0.0)
            atol,rtol = (PRED_ATOL,PRED_RTOL) if prediction else (SCORE_ATOL,SCORE_RTOL)
            if np.any(d > atol+rtol*np.maximum(np.abs(a),np.abs(b))):
                self.failures.append({"check":name,"kind":"numerical_mismatch","max_abs":float(np.max(d)),"atol":atol,"rtol":rtol})
        except (TypeError,ValueError,OverflowError) as error:
            self.failures.append({"check":name,"kind":"malformed_value","error":str(error)})


def rank_rho(x,y):
    if len(x)<2 or len(x)!=len(y):
        return None
    # Independent scalar rank definition, including exact ties.
    def ranks(a):
        sorted_values = sorted(a)
        from bisect import bisect_left,bisect_right
        return [(bisect_left(sorted_values,v)+bisect_right(sorted_values,v)+1)/2 for v in a]
    a,b = ranks(x),ranks(y)
    ma,mb = statistics.mean(a),statistics.mean(b)
    d = math.sqrt(math.fsum((v-ma)**2 for v in a)*math.fsum((v-mb)**2 for v in b))
    return math.fsum((u-ma)*(v-mb) for u,v in zip(a,b))/d if d else None


def validate_rows(rows,core):
    fields = {"GeneID","group_ID","HL_WT","LOG2_FC_IP","SE_IP","CDS_length_nt","CDS_GC","UTR5_length_nt","UTR5_GC","UTR3_length_nt","UTR3_GC"}
    fields |= {"codon_"+c for c in core.CODONS}|{"aa_"+a for a in core.AAS}
    if not isinstance(rows,list) or not rows:
        raise ValueError("Predictor list must be nonempty")
    for row in rows:
        if not isinstance(row,dict) or set(row)!=fields:
            raise ValueError("Exact predictor whitelist violated")
        if not all(isinstance(row[k],str) and row[k] for k in ("GeneID","group_ID")):
            raise ValueError("Malformed predictor identity")
        for k in fields-{"GeneID","group_ID"}:
            value=row[k]
            if k in ("UTR5_GC","UTR3_GC") and row[k[:4]+"_length_nt"]==0:
                if value is not None:
                    raise ValueError("Zero-length UTR GC must be null")
                continue
            if type(value) not in (int,float) or not math.isfinite(value):
                raise ValueError("Nonfinite/non-numeric predictor: "+k)
        if any(row[k]<=0 for k in ("HL_WT","SE_IP","CDS_length_nt")):
            raise ValueError("Nonpositive required predictor")
    ids=[r["GeneID"] for r in rows]
    if len(ids)!=len(set(ids)):
        raise ValueError("Duplicate predictor GeneID")
    return sorted(rows,key=lambda r:r["GeneID"])


def verify_predictions(core,train,test,y,artifact,checks):
    checks.exact("version",artifact["version"],VERSION)
    checks.exact("model_order",artifact["model_order"],core.MODELS)
    checks.exact("exact_model_keys",sorted(artifact["models"]),sorted(core.MODELS))
    checks.exact("training_membership",artifact["training_ids"],[r["GeneID"] for r in train])
    checks.exact("holdout_membership",artifact["holdout_ids"],[r["GeneID"] for r in test])
    checks.number("training_targets",artifact["training_logKO"],y)
    folds=np.array([core.fold(r["group_ID"],"inner") for r in train])
    checks.exact("inner_fold_assignments",artifact["training_fold_ids"],folds.tolist())
    x=np.log([r["HL_WT"] for r in train]); h=np.log([r["HL_WT"] for r in test])
    checks.number("baseline_quartiles",artifact["baseline_quartiles"],np.quantile(x,[.25,.5,.75],method="linear"))
    checks.number("h_ref",artifact["h_ref"],np.median([r["HL_WT"] for r in train]))
    for prefix,rows in [("training",train),("holdout",test)]:
        context=[{k:r[k] for k in ["GeneID","group_ID","HL_WT","LOG2_FC_IP"]} for r in rows]
        checks.exact(prefix+"_context",artifact[prefix+"_context"],context)
    for model in core.MODELS:
        actual=artifact["models"][model]
        if model=="no_change":
            pred=h
        elif model=="baseline_linear":
            z,w=core.standardize(x[:,None],h[:,None])
            coefficient=np.linalg.pinv(np.c_[np.ones(len(x)),z],rcond=1e-12)@y
            pred=np.c_[np.ones(len(h)),w]@coefficient
        else:
            oof={lam:np.empty(len(train)) for lam in core.PENALTIES}
            for f in range(5):
                mask=folds!=f
                tr=[r for i,r in enumerate(train) if mask[i]]
                va=[r for i,r in enumerate(train) if not mask[i]]
                if not tr or not va:
                    raise ValueError("Empty inner fitting/assessment fold")
                if any(not any(r[p+"_length_nt"]>0 for r in tr) for p in ("UTR5","UTR3")):
                    raise ValueError("No training nonempty UTR")
                a,b=core.transformed(tr,va,model)
                record=actual["fold_fits"][f]
                checks.exact(model+f"_fold{f}_training",record["training_ids"],[r["GeneID"] for r in tr])
                checks.exact(model+f"_fold{f}_assessment",record["assessment_ids"],[r["GeneID"] for r in va])
                for lam in core.PENALTIES:
                    oof[lam][~mask]=core.gram_prediction(a,y[mask],b,lam)
            losses={p:float(np.mean(abs(v-y))) for p,v in oof.items()}
            chosen=max(p for p in core.PENALTIES if losses[p]<=min(losses.values())+1e-12)
            checks.exact(model+"_selected_penalty",actual["selected_lambda"],chosen)
            for lam in core.PENALTIES:
                checks.number(model+"_OOF_"+str(lam),actual["oof_predictions"][str(lam)],oof[lam],prediction=True)
                checks.number(model+"_CV_loss_"+str(lam),actual["pooled_oof_mae"][str(lam)],losses[lam])
            a,b=core.transformed(train,test,model)
            pred=core.gram_prediction(a,y,b,chosen)
        checks.number(model+"_holdout_predictions",actual["holdout_prediction"],pred,prediction=True)


def losses(pred,y,groups):
    signed=np.asarray(pred)-y; absolute=np.abs(signed)
    return {"MAE":float(np.mean(absolute)),"RMSE":float(np.sqrt(np.mean(signed*signed))),
            "median_absolute_error":float(np.median(absolute)),"p90_absolute_error":float(np.quantile(absolute,.9,method="linear")),
            "mean_signed_error":float(np.mean(signed)),"spearman":rank_rho(list(pred),list(y)),
            "equal_component_MAE":statistics.mean(float(np.mean(absolute[np.array(groups)==g])) for g in sorted(set(groups)))}


def compare_entry(checks,name,actual,base,proposed):
    checks.number(name+"_absolute",actual["absolute_MAE_reduction"],base-proposed)
    checks.number(name+"_relative",actual["relative_MAE_reduction"],(base-proposed)/base if base>0 else None)
    if "meets_5_percent" in actual:
        checks.exact(name+"_criterion",actual["meets_5_percent"],bool(base>0 and proposed<=.95*base))


def verify_kinetic(checks,name,actual,context,y,h_ref):
    b=np.array([r["HL_WT"] for r in context]);ip=np.array([r["LOG2_FC_IP"] for r in context])
    arrays={"baseline_HL":b,"log2_lifetime_change":(y-np.log(b))/math.log(2),
            "dimensionless_rate_contrast":math.log(2)*h_ref*(1/b-1/np.exp(y)),"LOG2_FC_IP":ip}
    checks.exact(name+"_n",actual["n"],len(y));checks.number(name+"_h_ref",actual["h_ref"],h_ref)
    checks.exact(name+"_row_identifiers",[r["GeneID"] for r in actual["rows"]],[r["GeneID"] for r in context])
    for key,v in arrays.items():
        frozen_values=np.array([r[key] for r in actual["rows"]])
        checks.number(name+key+"_rows",frozen_values,v)
        s=actual["summaries"][key]
        for quantile,label in [(0.5,"median"),(.1,"p10"),(.9,"p90")]:
            checks.number(name+key+label,s[label],float(np.quantile(v,quantile,method="linear")))
        # Discrete signs describe the archived floating values. Numerically
        # equivalent reciprocal expressions can round differently at zero.
        for sign,n in [("negative",sum(frozen_values<0)),("zero",sum(frozen_values==0)),("positive",sum(frozen_values>0))]:
            checks.exact(name+key+sign,s[sign],int(n))
    for a,b in __import__("itertools").combinations(arrays,2):
        checks.number(name+a+b,actual["spearman"][a+"__"+b],rank_rho(list(arrays[a]),list(arrays[b])))


def verify_evaluation(core,artifact,evaluation,target,checks):
    context=artifact["holdout_context"];ids=artifact["holdout_ids"]
    if set(target)!=set(ids):raise ValueError("Heldout target membership mismatch")
    y=np.array([target[g] for g in ids]);groups=[r["group_ID"] for r in context]
    if not np.all(np.isfinite(y)):
        raise ValueError("Nonfinite heldout target")
    checks.exact("evaluation_version",evaluation["version"],VERSION)
    checks.exact("evaluation_status",evaluation["status"],artifact["status"])
    checks.exact("evaluation_model_keys",sorted(evaluation["metrics"]),sorted(core.MODELS))
    pred={m:np.array(artifact["models"][m]["holdout_prediction"]) for m in core.MODELS}
    metric={m:losses(v,y,groups) for m,v in pred.items()}
    for m in core.MODELS:
        for k,v in metric[m].items():checks.number(m+"_metric_"+k,evaluation["metrics"][m][k],v)
    baseline=core.MODELS[:5];proposed=metric["association"]["MAE"]
    checks.exact("primary_success",evaluation["primary_success"],all(metric[m]["MAE"]>0 and proposed<=.95*metric[m]["MAE"] for m in baseline))
    for m in baseline:compare_entry(checks,m,evaluation["primary_comparisons"][m],metric[m]["MAE"],proposed)
    for m,b in [("association_noSE","codon_noSE"),("association_nonlinear","codon_core")]:
        checks.exact(m+"_sensitivity_comparator",evaluation["sensitivities"][m]["baseline"],b)
        compare_entry(checks,m,evaluation["sensitivities"][m],metric[b]["MAE"],metric[m]["MAE"])
    bins=np.searchsorted(np.array(artifact["baseline_quartiles"]),np.log([r["HL_WT"] for r in context]),side="left")
    for f in range(4):
        mask=bins==f;entry=evaluation["baseline_strata"][f]
        checks.exact("strata"+str(f)+"_n",entry["n"],int(sum(mask)))
        q=artifact["baseline_quartiles"]
        checks.number("strata"+str(f)+"_lower",entry["lower"],None if f==0 else q[f-1])
        checks.number("strata"+str(f)+"_upper",entry["upper"],None if f==3 else q[f])
        checks.exact("strata"+str(f)+"_interval",entry["convention"],"lower open, upper closed")
        for m in core.MODELS:checks.number(m+"_stratum"+str(f),entry["MAE"][m],float(np.mean(abs(pred[m][mask]-y[mask]))) if any(mask) else None)
    checks.exact("evaluation_row_ids",[r["GeneID"] for r in evaluation["rows"]],ids)
    for i,row in enumerate(evaluation["rows"]):
        checks.number(f"row{i}_target",row["log_HL_KO"],y[i]);checks.exact(f"row{i}_bin",row["baseline_interval"],int(bins[i]))
        checks.exact(f"row{i}_group",row["group_ID"],groups[i])
        for m in core.MODELS:
            for key,v in [("predicted",pred[m][i]),("signed_error",pred[m][i]-y[i]),("absolute_error",abs(pred[m][i]-y[i]))]:
                checks.number(f"row{i}_{m}_{key}",row["models"][m][key],float(v))
    deletions=evaluation["component_deletions"];group_ids=sorted(set(groups))
    checks.exact("deleted_group_ids",[d["deleted_group"] for d in deletions],group_ids)
    gain=[];absolute=[];success=[];positive_all=[]
    for d,g in zip(deletions,group_ids):
        keep=np.array(groups)!=g
        errors={m:float(np.mean(abs(v[keep]-y[keep]))) for m,v in pred.items()}
        p=errors["association"]
        checks.exact("deletion"+g+"_n",d["remaining_genes"],int(sum(keep)))
        for m in core.MODELS[:6]:checks.number("deletion"+g+m+"_MAE",d["MAE"][m],errors[m])
        for m in baseline:compare_entry(checks,"deletion"+g+m,d["comparisons"][m],errors[m],p)
        ok=all(errors[m]>0 and p<=.95*errors[m] for m in baseline);success.append(ok)
        checks.exact("deletion"+g+"_success",d["primary_success"],ok)
        b=errors["codon_core"];absolute.append(b-p);gain.append((b-p)/b if b>0 else None)
        positive_all.append(all(errors[m]-p>0 for m in baseline))
    inf=evaluation["influence"]
    checks.number("influence_absolute_range",inf["codon_absolute_gain_range"],[min(absolute),max(absolute)])
    checks.number("influence_relative_range",inf["codon_relative_gain_range"],[min(gain),max(gain)] if all(v is not None for v in gain) else None)
    checks.exact("influence_undefined",inf["undefined_relative_deletions"],sum(v is None for v in gain))
    checks.exact("influence_codon_positive",inf["positive_codon_gain_every_deletion"],all(v>0 for v in absolute))
    checks.exact("influence_all_positive",inf["positive_gain_every_baseline_every_deletion"],all(positive_all))
    checks.exact("influence_all_5percent",inf["five_percent_every_deletion"],all(success))
    verify_kinetic(checks,"holdout",evaluation["kinetic"]["holdout"],context,y,artifact["h_ref"])
    verify_kinetic(checks,"whole",evaluation["kinetic"]["entire_eligible"],artifact["training_context"]+context,np.r_[artifact["training_logKO"],y],artifact["h_ref"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ["training-features","training-targets","holdout-features","design-receipt","predictions","protocol","author-code","authorization","out"]:
        parser.add_argument("--"+key,required=True)
    for key in ["evaluation","holdout-targets","commit-receipt"]:parser.add_argument("--"+key)
    args=parser.parse_args();out=Path(args.out)
    if out.exists():raise FileExistsError("Verifier output must be fresh")
    # Authorization metadata and hashes precede any scientific JSON decoding.
    authorization=load(args.authorization)
    status=authorization.get("status")
    if status not in ("synthetic_only","root_authorized_actual_output_verification"):
        raise ValueError("Separate verifier authorization required")
    if authorization.get("verifier_sha256")!=sha(__file__) or authorization.get("reference_core_sha256")!=CORE_SHA or sha(HERE/"reference_core.py")!=CORE_SHA:
        raise ValueError("Verifier/core hash mismatch")
    runtime={"implementation":platform.python_implementation(),"python":platform.python_version(),"numpy":np.__version__,"Biopython":Bio.__version__}
    if runtime!={"implementation":"CPython","python":"3.12.14","numpy":"2.5.4","Biopython":"1.88"}:
        raise RuntimeError("Frozen verifier runtime differs")
    paths={k:getattr(args,k) for k in ["training_features","training_targets","holdout_features","design_receipt","predictions","protocol","author_code"]}
    scoring=any(getattr(args,k) for k in ["evaluation","holdout_targets","commit_receipt"])
    if scoring:
        if not all(getattr(args,k) for k in ["evaluation","holdout_targets","commit_receipt"]):raise ValueError("All scoring paths required together")
        paths.update({k:getattr(args,k) for k in ["evaluation","holdout_targets","commit_receipt"]})
    before={k:sha(p) for k,p in paths.items()}
    if authorization.get("inputs")!=before or before["protocol"]!=PROTOCOL_SHA or before["author_code"]!=AUTHOR_SHA:
        raise ValueError("Frozen input/protocol/author code mismatch before decoding")
    spec=importlib.util.spec_from_file_location("independent_reference_core",HERE/"reference_core.py")
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    checks=Checks();out.mkdir(parents=True)
    try:
        raw_train,raw_test=load(args.training_features),load(args.holdout_features)
        train,test=validate_rows(raw_train,core),validate_rows(raw_test,core)
        targets=load(args.training_targets);design=load(args.design_receipt);artifact=load(args.predictions)
        if set(targets)!={r["GeneID"] for r in train}:raise ValueError("Training target membership mismatch")
        if status=="synthetic_only":
            if artifact.get("status")!="synthetic_only" or any(not r["GeneID"].startswith("SYN") for r in train+test):raise ValueError("Synthetic authorization cannot decode biological fixtures")
        elif artifact.get("status")!="root_authorized_operational_prediction":raise ValueError("Actual prediction authorization mismatch")
        checks.exact("design_artifact",artifact["design"]["source_receipt"],design)
        if {r["GeneID"] for r in train}&{r["GeneID"] for r in test} or {r["group_ID"] for r in train}&{r["group_ID"] for r in test}:raise ValueError("Train/holdout overlap")
        if any(core.fold(r["group_ID"],"outer")==0 for r in train) or any(core.fold(r["group_ID"],"outer")!=0 for r in test):raise ValueError("Frozen outer partition differs")
        y=np.array([targets[r["GeneID"]] for r in train],dtype=float)
        if not np.all(np.isfinite(y)):raise ValueError("Invalid training log targets")
        for key in ["training_features","training_targets","holdout_features","design_receipt","protocol"]:
            checks.exact("prediction_provenance_"+key,artifact["file_provenance"]["inputs"][key],before[key])
        checks.exact("prediction_author_code",artifact["file_provenance"]["implementation_sha256"],AUTHOR_SHA)
        verify_predictions(core,train,test,y,artifact,checks)
        if scoring:
            commit=load(args.commit_receipt)
            if commit.get("predictions_sha256")!=before["predictions"] or commit.get("holdout_targets_sha256")!=before["holdout_targets"]:raise ValueError("Commit receipt file mismatch")
            if status!="synthetic_only" and (commit.get("status")!="root_verified_operational_prediction_commit" or commit.get("verified") is not True or not re.fullmatch("[0-9a-f]{40}",str(commit.get("commit","")))):raise ValueError("Root-verified prediction commit required")
            evaluation=load(args.evaluation)
            for k,role in [("predictions_sha256","predictions"),("holdout_targets_sha256","holdout_targets"),("commit_receipt_sha256","commit_receipt")]:
                checks.exact("evaluation_provenance_"+k,evaluation["provenance"][k],before[role])
            checks.exact("evaluation_author_code",evaluation["provenance"]["implementation_sha256"],AUTHOR_SHA)
            verify_evaluation(core,artifact,evaluation,load(args.holdout_targets),checks)
        after={k:sha(p) for k,p in paths.items()}
        checks.exact("inputs_immutable",after,before)
        result={"status":"passed" if not checks.failures else "failed","authorization_status":status,"input_sha256":before,
                "verifier_sha256":sha(__file__),"reference_core_sha256":CORE_SHA,"runtime":runtime,
                "checks":checks.count,"numerical_values_compared":checks.values,"maximum_abs_difference":checks.maximum_abs,"failures":checks.failures,
                "prediction_vectors":9,"OOF_vectors":35,"ridge_penalty_selections":7,"evaluation_checked":scoring,
                "tolerances":{"prediction_atol":PRED_ATOL,"prediction_rtol":PRED_RTOL,"score_atol":SCORE_ATOL,"score_rtol":SCORE_RTOL},
                "limits":"Author code hashed, never imported. No source eligibility/protein-group recomputation. Historical Git not reauthenticated. Named algorithms only; not biological validity or causal calibration."}
        (out/"verification.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print(json.dumps({k:result[k] for k in ["status","checks","numerical_values_compared","maximum_abs_difference","evaluation_checked"]}))
        if checks.failures:raise RuntimeError("Independent verification discrepancies retained")
    except Exception as error:
        (out/"failure.json").write_text(json.dumps({"error_type":type(error).__name__,"message":str(error),"checks_before_failure":checks.count,"discrepancies":checks.failures},indent=2)+"\n",encoding="utf-8")
        raise


if __name__=="__main__":
    main()
