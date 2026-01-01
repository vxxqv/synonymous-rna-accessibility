# AUTHOR NOTEBOOK CELL 0
!which python
!python --version

# AUTHOR NOTEBOOK CELL 1
import pandas as pd
import seaborn as sns
import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

# AUTHOR NOTEBOOK CELL 2
basepath = "/home/commander-data/Storage/users/mave/search_sequences_300_var"
datapath = '/'.join([basepath, "data/input/fastqs/X205SC23111498-Z01-F020/01.RawData"])
sspath = '/'.join([basepath, "data/input/search_sequences"])
plotpath = '/'.join([basepath, "data/output/plots"])

# AUTHOR NOTEBOOK CELL 3
# read search sequences
search_sequences = []

for filepath in Path(sspath).iterdir():
    if filepath.is_file():
        new_data = pd.read_csv(filepath, sep = " ", header = None)
        new_data["modality"] = filepath.name.split('_')[2].split('.')[0]
        search_sequences.append(new_data)

search_sequences = pd.concat(search_sequences)
search_sequences.columns = ["combined_var_id", "seq", "modality"]
search_sequences

# AUTHOR NOTEBOOK CELL 4
# deduplicate search sequences
search_sequences["_is_old"] = search_sequences["combined_var_id"].str.endswith("_old")

tmp = (
    search_sequences.sort_values(["combined_var_id", "_is_old"])  # False comes before True
      .drop_duplicates(subset=["seq", "modality"], keep="first")
      .drop(columns="_is_old")
)
search_sequences = tmp.copy()

# AUTHOR NOTEBOOK CELL 5


# AUTHOR NOTEBOOK CELL 6
# read in counts
result = []

for subdir in Path(datapath).iterdir():
    if subdir.is_dir() and "copy" not in str(subdir) and "OLD" not in str(subdir) and "_after_" not in str(subdir) and "_before_" not in str(subdir):
        matches = list(subdir.glob("*_results_total.tsv"))
        if matches:
            if len(matches) > 1:
                print("multiple")
            for filepath in matches:
                new_data = pd.read_csv(filepath, sep = "\t", comment = "#")
                new_data["file"] = filepath.name
                result.append(new_data)
        else:
            print("NOTHING FOUND: " + str(subdir))
result = pd.concat(result)

result["modality"] = [x.split('_')[0] for x in result["file"]]
result["exon"] = [x.split('_')[1] for x in result["file"]]
result["cell_line"] = [x.split('_')[2] for x in result["file"]]
result["replicate"] = [x.split('_')[3] for x in result["file"]]
result = result.drop_duplicates()

result["fraction"] = (result["counts_complete"] +1) / (result["total_reads"] +1)
result["RPM"] = result["fraction"] * 1000000

result = result.merge(search_sequences, on = ["seq", "modality"], how = "left")
result

# AUTHOR NOTEBOOK CELL 7
# pivot such that each combined variant has the rpm dna and rpm rna in one line
wide = (
    result.pivot_table(
        index=["combined_var_id", "exon", "cell_line", "replicate"],
        columns="modality",
        values="RPM"
    ).reset_index().rename_axis(None, axis=1).rename(columns={"DNA": "RPM_DNA", "RNA": "RPM_RNA"})
)
wide

# AUTHOR NOTEBOOK CELL 8
# annotate library information
library = pd.read_csv("/home/commander-data/Storage/users/mave/cPRIME/multi_exon_screen/Overview_edit_names.tsv", sep = "\t")

wide_annotated = wide.merge(library, on = "combined_var_id")
wide_annotated["clinvar_simple"] = wide_annotated["clinvar"].map({
    "unknown": "unknown",
    "Uncertain_significance": "VUS",
    "Likely_benign": "benign",
    "Conflicting_classifications_of_pathogenicity": "conflicting",
    "Pathogenic": "pathogenic",
    "Likely_pathogenic": "pathogenic",
    "Pathogenic/Likely_pathogenic": "pathogenic",
    "Benign/Likely_benign": "benign"
})
wide_annotated["cprime_score"] = (wide_annotated["RPM_RNA"]) / (wide_annotated["RPM_DNA"])

# AUTHOR NOTEBOOK CELL 9
combined_spliceai_scores_path = "/home/commander-data/Storage/users/mave/search_sequences/data/library_design/combined_variants_spliceai.vcf"

combined_splice_scores = []

with open(combined_spliceai_scores_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or line == "":
            continue

        parts = line.split('\t')
        info = parts[7]
        info_parts = info.split(';')

        spliceai_gain_score = None
        spliceai_gain_pos = None
        spliceai_loss_score = None
        spliceai_loss_pos = None

        for info_part in info_parts:
            if info_part.startswith('SpliceAI='):
                if ',' in info_part:
                    print("MULTIPLE!")
                spliceai_parts = info_part[9:].split('|')
                #ALLELE|NAME|STRAND|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL|DS_AG_REF|DS_AG_ALT|DS_AL_REF|DS_AL_ALT|DS_DG_REF|DS_DG_ALT|DS_DL_REF|DS_DL_ALT

                acceptor_gain_score = float(spliceai_parts[3])
                acceptor_loss_score = float(spliceai_parts[4])
                donor_gain_score = float(spliceai_parts[5])
                donor_loss_score = float(spliceai_parts[6])
                acceptor_gain_pos = float(spliceai_parts[7])
                acceptor_loss_pos = float(spliceai_parts[8])
                donor_gain_pos = float(spliceai_parts[9])
                donor_loss_pos = float(spliceai_parts[10])

                spliceai_gain_score = donor_gain_score
                spliceai_gain_pos = donor_gain_pos
                if acceptor_gain_score > donor_gain_score:
                    spliceai_gain_score = acceptor_gain_score
                    spliceai_gain_pos = acceptor_gain_pos

                spliceai_loss_score = donor_loss_score
                spliceai_loss_pos = donor_loss_pos
                if acceptor_loss_score > donor_loss_score:
                    spliceai_loss_score = acceptor_loss_score
                    spliceai_loss_pos = acceptor_loss_pos
            if info_part.startswith('var_id='):
                var_id = info_part[7:]
            if info_part.startswith('marker_var_id='):
                marker_var_id = info_part[14:]

        combined_splice_scores.append({
            "combined_var_id": '_'.join([var_id, marker_var_id]),
            "combined_spliceai_gain_score": spliceai_gain_score,
            "combined_spliceai_gain_pos": spliceai_gain_pos,
            "combined_spliceai_loss_score": spliceai_loss_score,
            "combined_spliceai_loss_pos": spliceai_loss_pos,
            "combined_spliceai": '|'.join(spliceai_parts)
        })


combined_splice_scores = pd.DataFrame(combined_splice_scores)
combined_splice_scores = combined_splice_scores.drop_duplicates()

wide_annotated = wide_annotated.merge(combined_splice_scores, on = "combined_var_id", how = "left")
wide_annotated

# AUTHOR NOTEBOOK CELL 10
wide_annotated.to_csv("/home/commander-data/Storage/users/mave/search_sequences_300_var/data/output/wide_annotated_new.tsv", sep = "\t", index = False)

savedat = wide_annotated.copy()
savedat = savedat[[
    "combined_var_id", "var_id", "hgvsc", "hgvsp", "exon_number", "syn_var", "cell_line", 
    "RPM_DNA", "RPM_RNA", "cprime_score", 
    "pangolin_gain_score_combined", "pangolin_gain_pos_combined", "pangolin_loss_score_combined", "pangolin_loss_pos_combined",
    "pangolin_gain", "pangolin_gain_pos", "pangolin_loss", "pangolin_loss_pos", # only intended edit
    "spliceai_gain_score", "spliceai_gain_pos_y", "spliceai_loss_score", "spliceai_loss_pos_y", # only intended edit
    "combined_spliceai_gain_score", "combined_spliceai_gain_pos", "combined_spliceai_loss_score", "combined_spliceai_loss_pos",
    "consequence", "clinvar", "clinvar_simple", "pangolin_gain_combined_bin", "pangolin_gain_bin", "source",
    "Target_sequence", "Spacer_sequence", "PAM_sequence", "Extension_sequence", "Strand", "PBS_length", "PBS_GC_content", "RTT_length", "RTT_GC_content", "RHA_length", "longest_T_stretch", "complete_peg_rna"
]]
savedat.to_csv("/home/commander-data/Storage/users/mave/search_sequences_300_var/data/output/rpm_new.tsv", sep = "\t", index = False)

# AUTHOR NOTEBOOK CELL 11


# AUTHOR NOTEBOOK CELL 12
# take mean of replicates
wide_mean = wide_annotated[["combined_var_id", "cell_line", "RPM_DNA", "RPM_RNA"]].groupby(["combined_var_id", "cell_line"]).mean().reset_index()
wide_mean["cprime_score"] = (wide_mean["RPM_RNA"]) / (wide_mean["RPM_DNA"])
wide_mean = wide_mean.merge(library, on = "combined_var_id", how = "left")
# fill missing clinvar
missing_values = pd.DataFrame([
    {"combined_var_id": "chr17-31233005-T-G_old", "clinvar": "Pathogenic"},
    {"combined_var_id": "chr17-31233016-A-T_old", "clinvar": "unknown"},
    {"combined_var_id": "chr17-31233072-A-G_old", "clinvar": "Conflicting_classifications_of_pathogenicity"},
    {"combined_var_id": "chr17-31233123-G-A_old", "clinvar": "Likely_benign"},
    {"combined_var_id": "chr17-31233126-A-G_old", "clinvar": "Conflicting_classifications_of_pathogenicity"},
    {"combined_var_id": "chr17-31233128-T-A_old", "clinvar": "Pathogenic"},
    {"combined_var_id": "chr17-31233133-G-T_old", "clinvar": "Pathogenic"},
    {"combined_var_id": "chr17-31233156-T-G_old", "clinvar": "Likely_pathogenic"},
    {"combined_var_id": "chr17-31233157-C-T_old", "clinvar": "Pathogenic"},
])
clinvar_map = missing_values.set_index('combined_var_id')['clinvar']
wide_mean['clinvar'] = wide_mean['clinvar'].fillna(
    wide_mean['combined_var_id'].map(clinvar_map)
)
wide_mean["clinvar_simple"] = wide_mean["clinvar"].map({
    "unknown": "unknown",
    "Uncertain_significance": "VUS",
    "Likely_benign": "benign",
    "Conflicting_classifications_of_pathogenicity": "conflicting",
    "Pathogenic": "pathogenic",
    "Likely_pathogenic": "pathogenic",
    "Pathogenic/Likely_pathogenic": "pathogenic",
    "Benign/Likely_benign": "benign"
})
wide_mean

# AUTHOR NOTEBOOK CELL 13
wide_mean.to_csv("/home/commander-data/Storage/users/mave/search_sequences_300_var/data/output/wide_mean_new.tsv", sep = "\t", index = False)

savedat = wide_mean.copy()
savedat = savedat[[
    "combined_var_id", "var_id", "hgvsc", "hgvsp", "exon_number", "syn_var", "cell_line", 
    "RPM_DNA", "RPM_RNA", "cprime_score", 
    "pangolin_gain_score_combined", "pangolin_gain_pos_combined", "pangolin_loss_score_combined", "pangolin_loss_pos_combined",
    "pangolin_gain", "pangolin_gain_pos", "pangolin_loss", "pangolin_loss_pos",
    "spliceai_gain_score", "spliceai_gain_pos_y", "spliceai_loss_score", "spliceai_loss_pos_y",
    "consequence", "clinvar", "clinvar_simple", "pangolin_gain_combined_bin", "pangolin_gain_bin", "source",
    "Target_sequence", "Spacer_sequence", "PAM_sequence", "Extension_sequence", "Strand", "PBS_length", "PBS_GC_content", "RTT_length", "RTT_GC_content", "RHA_length", "longest_T_stretch", "complete_peg_rna"
]]
savedat.to_csv("/home/commander-data/Storage/users/mave/search_sequences_300_var/data/output/rpm_mean_new.tsv", sep = "\t", index = False)

# AUTHOR NOTEBOOK CELL 14


# AUTHOR NOTEBOOK CELL 15
# compute heatmap of bad and good marker variants
for cell_line_oi in ["KBM7", "K562"]:
    plotdat = wide_mean.copy()
    plotdat = plotdat[plotdat["cell_line"] == cell_line_oi]
    plotdat = plotdat[plotdat["RPM_DNA"] >= 5]
    plotdat = plotdat[plotdat["source"].isin(["ss_lib", "bad_marker_lib"])]
    plotdat["cprime_score"] = plotdat["cprime_score"].clip(upper=2)#np.log2(plotdat["cprime_score"])
    plotdat["var_id"] = [x.split('_')[0] for x in plotdat["combined_var_id"]]
    plotdat["syn_var_id"] = [x.split('_')[1] for x in plotdat["combined_var_id"]]
    new_plotdat = []
    for group in plotdat.groupby("var_id"):
        group = group[1]
        new_data = []
        good_markers = group[group["source"] != "bad_marker_lib"].copy()
        if good_markers.shape[0] > 0:
            good_markers["syn_var_index"] = range(good_markers.shape[0])
            new_data.append(good_markers)
        bad_markers = group[group["source"] == "bad_marker_lib"].copy()
        if bad_markers.shape[0] > 0:
            bad_markers["syn_var_index"] =  [(x + 1) * -1 for x in range(bad_markers.shape[0])]
            new_data.append(bad_markers)
        new_plotdat.append(pd.concat(new_data))
    new_plotdat = pd.concat(new_plotdat)
    new_plotdat["syn_var_index"] = new_plotdat["syn_var_index"].astype(str)


    # plot the figure
    heatmap_df = new_plotdat.drop_duplicates().pivot(index="var_id", columns="syn_var_index", values="cprime_score")
    heatmap_df = heatmap_df[(~heatmap_df.isna()).sum(axis = 1) > 1] # remove rows with only one marker left

    new_plotdat["pangolin_annot"] = new_plotdat["pangolin_gain_score_combined"].astype(str)+ " / " + new_plotdat["pangolin_gain"].astype(str) + " / " + new_plotdat["RPM_DNA"].round(2).astype(str)
    pangolin_df = new_plotdat.pivot(index = "var_id", columns = "syn_var_index", values = "pangolin_annot")
    pangolin_df = pangolin_df.reindex_like(heatmap_df)

    plt.figure(figsize=(7, 5))
    sns.heatmap(
        heatmap_df,
        cmap = "coolwarm",
        center = 1,
        annot = pangolin_df.round(2),
        fmt = "",
        annot_kws = {"size": 8}
    )
    plt.title(cell_line_oi + " RNA/DNA with Pangolin scores (0, 1: good markers, -1 bad marker)")
    #plt.tight_layout()
    
    plt.savefig(
        '/'.join([plotpath, cell_line_oi + "_good_vs_bad_markers.pdf"]),
        format = "pdf",
        bbox_inches = "tight"
    )

    plt.show()

    source_data = new_plotdat.drop_duplicates()[["var_id", "syn_var_index", "cprime_score", "pangolin_gain_score_combined"]]
    source_data["syn_var_index"] =  source_data["syn_var_index"].replace({
        "-1": "confounding_marker",
        "0": "non_confounding_marker_1",
        "1": "non_confounding_marker_2"
    })
    source_data = source_data.pivot(index="var_id", columns="syn_var_index", values=["cprime_score", "pangolin_gain_score_combined"])
    source_data = source_data[(~source_data.isna()).sum(axis = 1) > 2]
    source_data.columns = ['_'.join(x) for x in source_data.columns]
    source_data.to_excel('/'.join([plotpath, cell_line_oi + "_good_vs_bad_markers.xlsx"]))
    print('/'.join([plotpath, cell_line_oi + "_good_vs_bad_markers.xlsx"]))


# AUTHOR NOTEBOOK CELL 16


# AUTHOR NOTEBOOK CELL 17
cell_lines = ["K562", "KBM7"]

for cell_line_oi in cell_lines:
        rocdat = wide_mean.copy()
        rocdat = rocdat[rocdat["cell_line"].isin([cell_line_oi]) & (rocdat["RPM_DNA"] >= 5) & (rocdat["source"].isin(["ss_lib", "stop_syn_gained_lib"])) & (rocdat["consequence"] != "stop_gained")].copy()
        print(rocdat["source"].value_counts())
        rocdat = rocdat[rocdat["clinvar_simple"].isin(["benign", "pathogenic"])] # differentiate between benign and pathogenic variants
        
        rocdat = rocdat[["var_id", "cprime_score", "clinvar_simple"]].groupby("var_id", as_index=False).agg(cprime_score=("cprime_score", "mean"), clinvar_simple=("clinvar_simple", "first")) # take mean per intended edit
        rocdat["label"] = rocdat["clinvar_simple"] == "pathogenic"

        print(rocdat.shape)
        print(rocdat["clinvar_simple"].value_counts())

        y_true = (rocdat["clinvar_simple"] == "pathogenic").astype(int)
        y_score = -rocdat["cprime_score"] # flip score because lower = more pathogenic

        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)

        # plot
        plt.figure()
        plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
        plt.plot([0, 1], [0, 1], linestyle="--")

        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve " + cell_line_oi)
        plt.legend()

        plt.savefig(
            '/'.join([plotpath, cell_line_oi + "_roc.pdf"]),
            format = "pdf",
            bbox_inches="tight"
        )

        plt.show()

        ## do some diagnostics

        # LOO
        aucs = []
        for i in range(len(rocdat)):
            df_loo = rocdat.drop(rocdat.index[i]).copy()
            auc_val = roc_auc_score(df_loo["clinvar_simple"], -df_loo["cprime_score"])
            aucs.append(auc_val)

        print("LOO AUC range:", min(aucs), max(aucs))


        # calculate confidence interval
        y_true = (rocdat["clinvar_simple"] == "pathogenic").astype(int)
        y_score = -rocdat["cprime_score"] # flip direction

        # observed AUC
        observed_auc = roc_auc_score(y_true, y_score)

        # bootstrap
        n_bootstraps = 1000
        bootstrapped_aucs = []
        rng = np.random.default_rng(42)
        for _ in range(n_bootstraps):
            indices = rng.choice(len(rocdat), len(rocdat), replace=True)

            if len(np.unique(y_true.iloc[indices])) < 2:
                continue  # skip invalid samples
            
            auc = roc_auc_score(y_true.iloc[indices], y_score.iloc[indices])
            bootstrapped_aucs.append(auc)

        bootstrapped_aucs = np.array(bootstrapped_aucs)
        ci_lower = np.percentile(bootstrapped_aucs, 2.5)
        ci_upper = np.percentile(bootstrapped_aucs, 97.5)

        print("AUC:", observed_auc)
        print("95% CI:", ci_lower, "-", ci_upper)

        # calculate p value
        n_permutations = 1000
        permuted_aucs = []
        for _ in range(n_permutations):
            shuffled = rng.permutation(y_true)
            auc = roc_auc_score(shuffled, y_score)
            permuted_aucs.append(auc)

        permuted_aucs = np.array(permuted_aucs)
        p_value = np.mean(permuted_aucs >= observed_auc)

        print("p-value:", p_value)

        rocdat.to_excel('/'.join([plotpath, cell_line_oi + "_roc.xlsx"]))



# AUTHOR NOTEBOOK CELL 18


# AUTHOR NOTEBOOK CELL 19

for column_oi in ["clinvar_simple"]:
    for cell_line in ["K562", "KBM7"]:
        plotdat = wide_mean[wide_mean["cell_line"].isin([cell_line]) & (wide_mean["RPM_DNA"] >= 5) & (wide_mean["source"].isin(["ss_lib", "stop_syn_gained_lib"])) & (wide_mean["consequence"] != "stop_gained")].copy()
        print(plotdat["source"].value_counts())

        print(len(plotdat))

        print(len(plotdat["var_id"].unique()))

        plotdat = plotdat[["var_id", "cprime_score", "clinvar_simple"]].groupby("var_id", as_index=False).agg(cprime_score=("cprime_score", "mean"),clinvar_simple=("clinvar_simple", "first")) # take mean per intended edit
        plotdat['cprime_score_clipped'] = plotdat['cprime_score'].clip(lower=0, upper=2) # clip to 0 - 2 -> nicer for visualization
        print(len(plotdat["var_id"].unique()))

        if column_oi == "clinvar_simple":
            plotdat = plotdat[~plotdat["clinvar_simple"].isin(["unknown"])]
        categories = sorted(plotdat[column_oi].dropna().unique())

        print(len(plotdat["var_id"].unique()))

        color_map = {
            "pathogenic": "#E88A8B",
            "benign": "#8FD18F",
            "conflicting": "#ECB576",
            "unknown": "#A5A5A5",
            "VUS": "#8FAFE3",
        }
        colors = [color_map[c] for c in categories]
        data = [plotdat.loc[plotdat[column_oi] == cat, "cprime_score_clipped"] for cat in categories]

        display_labels = {
            "pathogenic": "(likely) pathogenic",
            "benign": "(likely) benign",
            "conflicting": "conflicting",
            "VUS": "VUS",
            "unknown": "unknown"
        }
        legend_labels = [display_labels[c] for c in categories]

        plt.figure(figsize=(10, 6))
        plt.hist(
            data,
            bins=15,
            stacked=True,
            label=legend_labels,
            color=colors,
            edgecolor="black",
        )
        plt.xlabel("cprime_score_clipped")
        plt.ylabel("Count")
        plt.title(cell_line + " cPRIME score vs " + column_oi)
        plt.legend(title="ClinVar pathogenicity")
        plt.tight_layout()

        plt.savefig(
            '/'.join([plotpath, cell_line + "_cprime_score_vs_pangolin.pdf"]),
            format = "pdf",
            bbox_inches="tight"
        )

        plt.show()

        plotdat.to_excel('/'.join([plotpath, cell_line + "_cprime_score_vs_pangolin.xlsx"]))
