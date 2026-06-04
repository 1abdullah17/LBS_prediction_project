import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser


# =====================================================
# COMMAND LINE ARGUMENTS
# =====================================================

parser = argparse.ArgumentParser(
    description="Analyze AlphaFold P2Rank pocket predictions"
)

parser.add_argument(
    "--prediction_dir",
    default="results/alphafold",
    help="Directory containing P2Rank prediction files"
)

parser.add_argument(
    "--pdb_dir",
    default="data/alphafold_clean",
    help="Directory containing cleaned AlphaFold structures"
)

parser.add_argument(
    "--output",
    default="results/alphafold_enhanced_master_table.xlsx",
    help="Output Excel file"
)

args = parser.parse_args()


# =====================================================
# RESIDUE GROUPS
# =====================================================

HYDROPHOBIC = {
    "ALA", "VAL", "ILE", "LEU",
    "MET", "PHE", "TRP", "PRO"
}

CHARGED = {
    "ASP", "GLU",
    "LYS", "ARG", "HIS"
}


# =====================================================
# LOAD RESIDUE MAP FROM PDB
# =====================================================

def load_residue_map(pdb_file):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        pdb_file
    )

    residue_map = {}

    for model in structure:

        for chain in model:

            for residue in chain:

                if residue.id[0] != " ":
                    continue

                residue_number = residue.id[1]

                residue_name = residue.get_resname()

                residue_map[residue_number] = residue_name

    return residue_map


# =====================================================
# CALCULATE DESCRIPTORS
# =====================================================

def calculate_descriptors(
        residue_string,
        residue_map
):

    residue_entries = str(
        residue_string
    ).split()

    residue_names = []

    for item in residue_entries:

        try:

            residue_number = int(
                item.split("_")[1]
            )

            if residue_number in residue_map:

                residue_names.append(
                    residue_map[residue_number]
                )

        except Exception:
            continue

    total = len(residue_names)

    if total == 0:

        return (
            0,
            0,
            0
        )

    hydrophobic_count = sum(
        r in HYDROPHOBIC
        for r in residue_names
    )

    charged_count = sum(
        r in CHARGED
        for r in residue_names
    )

    hydrophobic_fraction = (
        hydrophobic_count / total
    )

    charged_fraction = (
        charged_count / total
    )

    return (
        total,
        hydrophobic_fraction,
        charged_fraction
    )


# =====================================================
# FIND PREDICTION FILES
# =====================================================

prediction_files = [

    f

    for f in os.listdir(
        args.prediction_dir
    )

    if f.endswith(
        "_predictions.csv"
    )
]

print(
    f"Found {len(prediction_files)} prediction files"
)


# =====================================================
# PROCESS FILES
# =====================================================

all_data = []

for file in prediction_files:

    prediction_path = os.path.join(
        args.prediction_dir,
        file
    )

    kinase_name = file.replace(
        "_clean.pdb_predictions.csv",
        ""
    )

    pdb_file = os.path.join(
        args.pdb_dir,
        kinase_name + "_clean.pdb"
    )

    if not os.path.exists(pdb_file):

        print(
            f"Missing structure: {kinase_name}"
        )

        continue

    try:

        residue_map = load_residue_map(
            pdb_file
        )

        df = pd.read_csv(
            prediction_path
        )

        df.columns = (
            df.columns
            .str.strip()
        )

        for _, row in df.iterrows():

            (
                num_residues,
                hydrophobic_fraction,
                charged_fraction

            ) = calculate_descriptors(

                row.get(
                    "residue_ids",
                    ""
                ),

                residue_map
            )

            all_data.append({

                "kinase":
                    kinase_name,

                "pocket_rank":
                    row.get(
                        "rank"
                    ),

                "score":
                    row.get(
                        "score"
                    ),

                "probability":
                    row.get(
                        "probability"
                    ),

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

                "pocket_size":
                    row.get(
                        "sas_points"
                    ),

                "center_x":
                    row.get(
                        "center_x"
                    ),

                "center_y":
                    row.get(
                        "center_y"
                    ),

                "center_z":
                    row.get(
                        "center_z"
                    )
            })

    except Exception as e:

        print(
            f"Error processing {file}: {e}"
        )


# =====================================================
# CREATE MASTER TABLE
# =====================================================

master_df = pd.DataFrame(
    all_data
)

master_df = master_df.sort_values(

    by=[
        "kinase",
        "pocket_rank"
    ]
)

master_df.to_excel(
    args.output,
    index=False
)

print(
    f"\nMaster table saved:\n{args.output}"
)


# =====================================================
# CREATE FIGURES DIRECTORY
# =====================================================

figure_dir = "results/week3_figures"

os.makedirs(
    figure_dir,
    exist_ok=True
)


# =====================================================
# BAR CHART
# =====================================================

top_pockets = master_df[
    master_df["pocket_rank"] == 1
]

plt.figure(
    figsize=(10, 5)
)

plt.bar(

    top_pockets["kinase"],

    top_pockets[
        "hydrophobic_fraction"
    ]
)

plt.xticks(
    rotation=90
)

plt.ylabel(
    "Hydrophobic Fraction"
)

plt.tight_layout()

plt.savefig(

    os.path.join(

        figure_dir,

        "hydrophobic_fraction_barplot.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# SCATTER 1
# =====================================================

plt.figure()

plt.scatter(

    master_df["score"],

    master_df["probability"]
)

plt.xlabel(
    "P2Rank Score"
)

plt.ylabel(
    "Probability"
)

plt.tight_layout()

plt.savefig(

    os.path.join(

        figure_dir,

        "score_vs_probability.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# SCATTER 2
# =====================================================

plt.figure()

plt.scatter(

    master_df["pocket_size"],

    master_df["score"]
)

plt.xlabel(
    "Pocket Size (SAS Points)"
)

plt.ylabel(
    "P2Rank Score"
)

plt.tight_layout()

plt.savefig(

    os.path.join(

        figure_dir,

        "size_vs_score.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# SCATTER 3
# =====================================================

plt.figure()

plt.scatter(

    master_df[
        "hydrophobic_fraction"
    ],

    master_df["score"]
)

plt.xlabel(
    "Hydrophobic Fraction"
)

plt.ylabel(
    "P2Rank Score"
)

plt.tight_layout()

plt.savefig(

    os.path.join(

        figure_dir,

        "hydrophobicity_vs_score.png"
    ),

    dpi=300
)

plt.close()


# =====================================================
# SCATTER 4
# =====================================================

plt.figure()

plt.scatter(

    master_df[
        "charged_fraction"
    ],

    master_df["score"]
)

plt.xlabel(
    "Charged Fraction"
)

plt.ylabel(
    "P2Rank Score"
)

plt.tight_layout()

plt.savefig(

    os.path.join(

        figure_dir,

        "charge_vs_score.png"
    ),

    dpi=300
)

plt.close()

print(
    "\nFigures saved to:",
    figure_dir
)