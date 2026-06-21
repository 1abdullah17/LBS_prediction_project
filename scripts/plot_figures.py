import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

# ==================================================
# INPUT
# ==================================================

INPUT_FILE = "results/cleaning_summary.xlsx"

OUTPUT_FIGURE = (
    "figures/alphafold_quality_assessment.tiff"
)

# ==================================================
# LOAD DATA
# ==================================================

df = pd.read_excel(INPUT_FILE)

print("Kinases:", len(df))

# ==================================================
# CREATE FIGURE
# ==================================================

fig = plt.figure(
    figsize=(18, 12)
)

gs = fig.add_gridspec(
    2,
    1,
    height_ratios=[1, 1.4]
)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# ==================================================
# PANEL A
# Mean pLDDT
# ==================================================

mean_plddt = df["avg_plddt"]

ax1.axvspan(
    0,
    70,
    alpha=0.15,
    color="lightcoral"
)

ax1.axvspan(
    70,
    90,
    alpha=0.15,
    color="khaki"
)

ax1.axvspan(
    90,
    100,
    alpha=0.15,
    color="lightgreen"
)

ax1.hist(
    mean_plddt,
    bins=15,
    edgecolor="black"
)

ax1.axvline(
    mean_plddt.mean(),
    linestyle="--",
    linewidth=2,
    label=f"Mean = {mean_plddt.mean():.1f}"
)

ax1.set_title(
    "A. Mean AlphaFold Confidence (pLDDT)",
    fontsize=18,
    weight="bold"
)

ax1.set_xlabel(
    "Mean pLDDT",
    fontsize=14,
    weight="bold"
)

ax1.set_ylabel(
    "Number of Kinases",
    fontsize=14,
    weight="bold"
)

ymax = ax1.get_ylim()[1]

ax1.text(
    65,
    ymax * 0.85,
    "Low\n(<70)",
    ha="center",
    fontsize=16,
    weight="bold"
)

ax1.text(
    85,
    ymax * 0.85,
    "Moderate\n(70–90)",
    ha="center",
    fontsize=16,
    weight="bold"
)

ax1.text(
    95,
    ymax * 0.85,
    "High\n(>90)",
    ha="center",
    fontsize=16,
    weight="bold"
)

stats_text = (
    f"Mean = {mean_plddt.mean():.1f}\n"
    f"Median = {mean_plddt.median():.1f}\n"
    f"Min = {mean_plddt.min():.1f}\n"
    f"Max = {mean_plddt.max():.1f}"
)

ax1.text(
    0.02,
    0.95,
    stats_text,
    transform=ax1.transAxes,
    fontsize=15,
    va="top",
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        alpha=1
    )
)

ax1.legend()

# ==================================================
# PANEL C
# Residue retention
# ==================================================

sorted_df = df.sort_values(
    "percent_removed"
).reset_index(drop=True)

x = np.arange(len(sorted_df))

ax2.bar(
    x,
    sorted_df["kept_residues"],
    label="Kept"
)

ax2.bar(
    x,
    sorted_df["removed_residues"],
    bottom=sorted_df["kept_residues"],
    label="Removed"
)

ax2.set_title(
    "B. Residue Retention After Cleaning",
    fontsize=18,
    weight="bold"
)

ax2.set_ylabel(
    "Number of Residues",
    fontsize=14,
    weight="bold"
)

ax2.set_xlabel(
    "Kinases (sorted by % residues removed)",
    fontsize=14,
    weight="bold"
)

if "af_id" in sorted_df.columns:

    labels = sorted_df["af_id"]

elif "kinase" in sorted_df.columns:

    labels = sorted_df["kinase"]

else:

    labels = range(len(sorted_df))

ax2.set_xticks(x)

ax2.set_xticklabels(
    labels,
    rotation=90,
    fontsize=8
)

ax2.legend()

# ==================================================
# SAVE
# ==================================================

os.makedirs(
    "figures",
    exist_ok=True
)

plt.tight_layout()

plt.savefig(
    OUTPUT_FIGURE,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(OUTPUT_FIGURE)