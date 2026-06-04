import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from scipy.stats import pearsonr


# =====================================================
# INPUT FILES
# =====================================================

POCKET_FILE = "results/alphafold_enhanced_master_table.xlsx"

CONFIDENCE_FILE = "data/cleaning_summary.xlsx"

OUTPUT_DIR = "results/week3_plddt_analysis"

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# =====================================================
# LOAD DATA
# =====================================================

pocket_df = pd.read_excel(
    POCKET_FILE
)

confidence_df = pd.read_excel(
    CONFIDENCE_FILE
)


# =====================================================
# MERGE DATASETS
# =====================================================

merged_df = pocket_df.merge(
    confidence_df,
    on="kinase",
    how="inner"
)

print(
    f"Merged rows: {len(merged_df)}"
)


# =====================================================
# CORRELATIONS
# =====================================================

results = []


def calculate_correlation(
        x,
        y,
        x_name,
        y_name
):

    r, p = pearsonr(x, y)

    results.append({

        "variable_x": x_name,
        "variable_y": y_name,
        "pearson_r": round(r, 4),
        "p_value": round(p, 6)

    })

    print(
        f"{x_name} vs {y_name}"
    )

    print(
        f"r = {r:.3f}"
    )

    print(
        f"p = {p:.5f}\n"
    )


calculate_correlation(

    merged_df["avg_plddt"],

    merged_df["score"],

    "avg_plddt",

    "score"

)

calculate_correlation(

    merged_df["avg_plddt"],

    merged_df["probability"],

    "avg_plddt",

    "probability"

)

calculate_correlation(

    merged_df["low_conf_percent"],

    merged_df["score"],

    "low_conf_percent",

    "score"

)

calculate_correlation(

    merged_df["low_conf_percent"],

    merged_df["probability"],

    "low_conf_percent",

    "probability"

)


# =====================================================
# SAVE CORRELATION TABLE
# =====================================================

correlation_df = pd.DataFrame(
    results
)

correlation_df.to_excel(

    os.path.join(
        OUTPUT_DIR,
        "correlation_statistics.xlsx"
    ),

    index=False
)


# =====================================================
# FIGURE 1
# SCORE VS AVG PLDDT
# =====================================================

plt.figure(figsize=(6, 5))

plt.scatter(

    merged_df["avg_plddt"],

    merged_df["score"]

)

plt.xlabel(
    "Average pLDDT"
)

plt.ylabel(
    "P2Rank Score"
)

plt.title(
    "Pocket Score vs Average pLDDT"
)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "score_vs_avg_plddt.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# FIGURE 2
# PROBABILITY VS AVG PLDDT
# =====================================================

plt.figure(figsize=(6, 5))

plt.scatter(

    merged_df["avg_plddt"],

    merged_df["probability"]

)

plt.xlabel(
    "Average pLDDT"
)

plt.ylabel(
    "P2Rank Probability"
)

plt.title(
    "Probability vs Average pLDDT"
)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "probability_vs_avg_plddt.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# FIGURE 3
# SCORE VS LOW-CONFIDENCE %
# =====================================================

plt.figure(figsize=(6, 5))

plt.scatter(

    merged_df["low_conf_percent"],

    merged_df["score"]

)

plt.xlabel(
    "Low-confidence Residues (%)"
)

plt.ylabel(
    "P2Rank Score"
)

plt.title(
    "Pocket Score vs Low-confidence Regions"
)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "score_vs_low_confidence.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# FIGURE 4
# BOXPLOT
# =====================================================

merged_df["confidence_group"] = merged_df[
    "low_conf_percent"
].apply(

    lambda x:
    "High confidence"
    if x < 20
    else "Low confidence"
)

high_conf_scores = merged_df.loc[
    merged_df["confidence_group"]
    == "High confidence",
    "score"
]

low_conf_scores = merged_df.loc[
    merged_df["confidence_group"]
    == "Low confidence",
    "score"
]

plt.figure(figsize=(6, 5))

plt.boxplot(

    [
        high_conf_scores,
        low_conf_scores
    ],

    labels=[
        "High confidence",
        "Low confidence"
    ]
)

plt.ylabel(
    "P2Rank Score"
)

plt.title(
    "Pocket Scores by Confidence Group"
)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "boxplot_confidence_groups.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# SUMMARY TABLE
# =====================================================

summary = merged_df.groupby(
    "confidence_group"
).agg({

    "score": ["mean", "std"],

    "probability": ["mean", "std"],

    "pocket_size": ["mean", "std"]

})

summary.to_excel(

    os.path.join(
        OUTPUT_DIR,
        "confidence_group_summary.xlsx"
    )
)

print(
    "\nAnalysis complete."
)

print(
    f"Results saved to: {OUTPUT_DIR}"
)