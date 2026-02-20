import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score

plt.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "Arial"
})

# Load data
df = pd.read_csv(
    "/Users/angela/Documents/hall_jiang_lab/ecSIR_Manuscript/roc_table2.tsv",
    sep="\t"
)

COPRO_CUTOFF = 6          # log10(area_coprostanol + 1) threshold
y_true = (np.log10(df['area_coprostanol'] + 1) >= COPRO_CUTOFF).astype(int)

species_cols = [
    "ismA_only_species_mapped",
    "spiR_and_ismA_species_mapped",
    "spiR_only_species_mapped",
    "spiR_or_ismA_species_mapped",
    "reads_mapped_spiR_species",
    "reads_mapped_ismA_species",
]

gene_cols = [
    "reads_mapped_spiR_gene",
    "reads_mapped_ismA_gene",
    "spiR_or_ismA_gene_mapped",
]

# Convert to CPM
# CPM = (counts / total_reads) * 1e6
for col in species_cols + gene_cols:
    df[f"{col}_cpm"] = df[col] / df["total_reads"] * 1e6

# Replace the original column lists with the new CPM ones
species_cols = [f"{col}_cpm" for col in species_cols]
gene_cols    = [f"{col}_cpm" for col in gene_cols]

# Curve colors, legend tables
STYLE = {
    # ------------- species -------------
    "reads_mapped_spiR_species_cpm" : ("spiR"            , "#387fb9"),
    "reads_mapped_ismA_species_cpm": ("ismA"            , "#f57f20"),
    "spiR_or_ismA_species_mapped_cpm":("spiR or ismA"    , "#4eaf49"),
    "spiR_and_ismA_species_mapped_cpm":("spiR and ismA"  , "#2d4196"),
    "ismA_only_species_mapped_cpm" : ("SpiR– ismA+"     , "#c63c40"),
    "spiR_only_species_mapped_cpm"  : ("SpiR+ ismA–"     , "#a65627"),

    # ------------- gene ---------------
    "reads_mapped_spiR_gene_cpm"    : ("spiR"            , "#387fb9"),
    "reads_mapped_ismA_gene_cpm"   : ("ismA"            , "#f57f20"),
    "spiR_or_ismA_gene_mapped_cpm"  : ("spiR or ismA"    , "#4eaf49"),
}


# Plot one ROC curve function
def add_roc(ax, y_true_series, score_series, label, colour, **plot_kw):
    mask = y_true_series.notna() & score_series.notna()
    y_cln, s_cln = y_true_series[mask].astype(int), score_series[mask]

    if len(np.unique(y_cln)) < 2:
        return  # need both classes present

    fpr, tpr, _ = roc_curve(y_cln, s_cln)
    auc = roc_auc_score(y_cln, s_cln)
    ax.plot(fpr, tpr, label=f"{label} (AUC={auc:.4f})",
            color=colour, lw=2, **plot_kw)


def add_pr(ax, y_true_series, score_series, label, colour, **plot_kw):
    mask = y_true_series.notna() & score_series.notna()
    y_cln, s_cln = y_true_series[mask].astype(int), score_series[mask]

    if len(np.unique(y_cln)) < 2:
        return

    prec, rec, _ = precision_recall_curve(y_cln, s_cln)
    ap = average_precision_score(y_cln, s_cln)
    ax.plot(rec, prec, label=f"{label} (AP={ap:.4f})",
            color=colour, lw=2, **plot_kw)



# Draw panels
fig, (ax_species, ax_gene) = plt.subplots(
    1, 2, figsize=(12, 6), sharex=True, sharey=True
)

# Species
for col in species_cols:
    lbl, col_hex = STYLE[col]
    add_roc(ax_species, y_true, df[col], lbl, col_hex)
ax_species.set_title("ROC: Species level")
ax_species.set_ylabel("True-Positive Rate")
ax_species.set_xlabel("False-Positive Rate")
ax_species.plot([0, 1], [0, 1], "k--", lw=1)
ax_species.grid(alpha=0.3)
ax_species.legend(fontsize=8)

# Gene
for col in gene_cols:
    lbl, col_hex = STYLE[col]
    add_roc(ax_gene, y_true, df[col], lbl, col_hex)
ax_gene.set_title("ROC: Gene level")
ax_gene.set_xlabel("False-Positive Rate")
ax_gene.plot([0, 1], [0, 1], "k--", lw=1)
ax_gene.grid(alpha=0.3)
ax_gene.legend(fontsize=8)

# save figure
fig.tight_layout()
fig.savefig("roc_species_gene.svg", format="svg")
plt.show()


# PR curve
fig_pr, (ax_pr_species, ax_pr_gene) = plt.subplots(
    1, 2, figsize=(12, 6), sharex=True, sharey=True
)

# Species level
for col in species_cols:
    lbl, col_hex = STYLE[col]
    add_pr(ax_pr_species, y_true, df[col], lbl, col_hex)
ax_pr_species.set_title("Precision-Recall: Species level")
ax_pr_species.set_xlabel("Recall")
ax_pr_species.set_ylabel("Precision")
ax_pr_species.grid(alpha=0.3)
ax_pr_species.legend(fontsize=8)

# Gene level
for col in gene_cols:
    lbl, col_hex = STYLE[col]
    add_pr(ax_pr_gene, y_true, df[col], lbl, col_hex)
ax_pr_gene.set_title("Precision-Recall: Gene level")
ax_pr_gene.set_xlabel("Recall")
ax_pr_gene.grid(alpha=0.3)
ax_pr_gene.legend(fontsize=8)

fig_pr.tight_layout()
fig_pr.savefig("pr_species_gene.svg", format="svg")
plt.show()

