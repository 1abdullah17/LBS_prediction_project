import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# =====================================================
# ARGUMENTS
# =====================================================

parser = argparse.ArgumentParser(
    description="Create P2Rank master table and plots"
)

parser.add_argument(
    "--input_dir",
    default="results/alphafold",
    help="Directory containing P2Rank output files"
)

parser.add_argument(
    "--output",
    default="results/alphafold_master_table.xlsx",
    help="Output Excel file"
)

parser.add_argument(
    "--plots_dir",
    default="results/plots",
    help="Directory for plots"
)

args = parser.parse_args()

INPUT_DIR = Path(args.input_dir)
OUTPUT_FILE = Path(args.output)
PLOTS_DIR = Path(args.plots_dir)

PLOTS_DIR.mkdir(parents=True, exist_ok=True)


# =====================================================
# RESIDUE CLASSES
# =====================================================

HYDROPHOBIC = {
    "ALA", "VAL", "ILE", "LEU",
    "MET", "PHE", "TRP", "PRO"
}

CHARGED = {
    "ASP", "GLU", "LYS",
    "ARG", "HIS"
}


# =====================================================
# DESCRIPTOR FUNCTION
# =====================================================

def calculate_descriptors(residue_ids_string, residue_lookup):

    if pd.isna(residue_ids_string):
        return 0, 0.0, 0.0

    residue_names = []

    tokens = str(residue_ids_string).split()

    for token in tokens:

        try:
            residue_number = int(
                token.split("_")[1]
            )

            if residue_number in residue_lookup:

                residue_names.append(
                    residue_lookup[residue_number]
                )

        except Exception:
            pass

    total = len(residue_names)

    if total == 0:
        return 0, 0.0, 0.0

    hydrophobic_count = sum(
        residue in HYDROPHOBIC
        for residue in residue_names
    )

    charged_count = sum(
        residue in CHARGED
        for residue in residue_names
    )

    return (
        total,
        hydrophobic_count / total,
        charged_count / total
    )


# =====================================================
# FIND FILES
# =====================================================

prediction_files = sorted(
    INPUT_DIR.glob("*_predictions.csv")
)

print(f"\nFound {len(prediction_files)} prediction files\n")

all_data = []


# =====================================================
# PROCESS FILES
# =====================================================

for pred_file in prediction_files:

    kinase_name = (
        pred_file.name
        .replace("_clean.pdb_predictions.csv", "")
        .replace("_predictions.csv", "")
    )

    residue_file = Path(
        str(pred_file).replace(
            "_predictions.csv",
            "_residues.csv"
        )
    )

    print(f"Processing {kinase_name}")

    try:

        pred_df = pd.read_csv(pred_file)
        pred_df.columns = pred_df.columns.str.strip()

        if not residue_file.exists():

            print(
                f"WARNING: Missing residue file: {residue_file.name}"
            )

            continue

        res_df = pd.read_csv(residue_file)

        # -----------------------------
        # CLEAN COLUMN NAMES
        # -----------------------------
        res_df.columns = res_df.columns.str.strip()

        # -----------------------------
        # CLEAN VALUES
        # THIS FIXES THE ZERO PROBLEM
        # -----------------------------
        res_df["residue_label"] = pd.to_numeric(
            res_df["residue_label"],
            errors="coerce"
        )

        res_df["residue_name"] = (
            res_df["residue_name"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        # -----------------------------
        # LOOKUP DICTIONARY
        # -----------------------------
        residue_lookup = dict(
            zip(
                res_df["residue_label"],
                res_df["residue_name"]
            )
        )

        # -----------------------------
        # POCKETS
        # -----------------------------
        for _, row in pred_df.iterrows():

            residue_ids = row.get(
                "residue_ids",
                ""
            )

            (
                num_residues,
                hydrophobic_fraction,
                charged_fraction
            ) = calculate_descriptors(
                residue_ids,
                residue_lookup
            )

            all_data.append({

                "kinase":
                    kinase_name,

                "pocket_rank":
                    row.get("rank"),

                "score":
                    row.get("score"),

                "probability":
                    row.get("probability"),

                "num_residues":
                    num_residues,

                "hydrophobic_fraction":
                    round(
                        hydrophobic_fraction,
                        3
                    ),

                "charged_fraction":
                    round(
                        charged_fraction,
                        3
                    ),

                "center_x":
                    row.get("center_x"),

                "center_y":
                    row.get("center_y"),

                "center_z":
                    row.get("center_z"),

                "surface_atoms":
                    row.get("surf_atoms"),

                "pocket_size":
                    row.get("sas_points"),

                "residue_ids":
                    residue_ids

            })

    except Exception as e:

        print(
            f"ERROR processing {pred_file.name}"
        )

        print(e)


# =====================================================
# MASTER TABLE
# =====================================================

master_df = pd.DataFrame(all_data)

if master_df.empty:

    print("\nNo data found.")
    quit()

master_df = master_df.sort_values(
    ["kinase", "pocket_rank"]
)

master_df.to_excel(
    OUTPUT_FILE,
    index=False
)

print(
    f"\nMaster table saved:\n{OUTPUT_FILE}"
)

print(
    f"Total pockets: {len(master_df)}"
)


# =====================================================
# BAR CHART
# NUMBER OF POCKETS
# =====================================================

plt.figure(figsize=(12, 6))

(
    master_df.groupby("kinase")
    .size()
    .sort_values(ascending=False)
    .plot(kind="bar")
)

plt.title("Number of Predicted Pockets per Kinase")
plt.ylabel("Pocket Count")
plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "bar_pockets_per_kinase.png",
    dpi=300
)

plt.close()


# =====================================================
# BAR CHART
# MEAN SCORE
# =====================================================

plt.figure(figsize=(12, 6))

(
    master_df.groupby("kinase")["score"]
    .mean()
    .sort_values(ascending=False)
    .plot(kind="bar")
)

plt.title("Mean P2Rank Score per Kinase")
plt.ylabel("Mean Score")
plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "bar_mean_score_per_kinase.png",
    dpi=300
)

plt.close()


# =====================================================
# SCATTER
# SCORE VS PROBABILITY
# =====================================================

plt.figure(figsize=(7, 6))

plt.scatter(
    master_df["score"],
    master_df["probability"],
    alpha=0.7
)

plt.xlabel("Score")
plt.ylabel("Probability")
plt.title("P2Rank Score vs Probability")

plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "scatter_score_probability.png",
    dpi=300
)

plt.close()


# =====================================================
# SCATTER
# POCKET SIZE VS SCORE
# =====================================================

plt.figure(figsize=(7, 6))

plt.scatter(
    master_df["pocket_size"],
    master_df["score"],
    alpha=0.7
)

plt.xlabel("Pocket Size (SAS points)")
plt.ylabel("Score")
plt.title("Pocket Size vs Score")

plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "scatter_pocket_size_score.png",
    dpi=300
)

plt.close()


# =====================================================
# SCATTER
# HYDROPHOBIC FRACTION VS SCORE
# =====================================================

plt.figure(figsize=(7, 6))

plt.scatter(
    master_df["hydrophobic_fraction"],
    master_df["score"],
    alpha=0.7
)

plt.xlabel("Hydrophobic Fraction")
plt.ylabel("Score")
plt.title("Hydrophobic Fraction vs Score")

plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "scatter_hydrophobic_fraction_score.png",
    dpi=300
)

plt.close()


# =====================================================
# SCATTER
# CHARGED FRACTION VS SCORE
# =====================================================

plt.figure(figsize=(7, 6))

plt.scatter(
    master_df["charged_fraction"],
    master_df["score"],
    alpha=0.7
)

plt.xlabel("Charged Fraction")
plt.ylabel("Score")
plt.title("Charged Fraction vs Score")

plt.tight_layout()

plt.savefig(
    PLOTS_DIR /
    "scatter_charged_fraction_score.png",
    dpi=300
)

plt.close()

print("\nPlots saved to:")
print(PLOTS_DIR)

print("\nDone.")