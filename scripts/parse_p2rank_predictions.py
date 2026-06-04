import os
import argparse
import pandas as pd


# =====================================================
# COMMAND LINE ARGUMENTS
# =====================================================

parser = argparse.ArgumentParser(
    description="Run P2Rank and create a master pocket table."
)

parser.add_argument(
    "--input_dir",
    default="results/alphafold",
    help="Directory containing P2Rank prediction files"
)

parser.add_argument(
    "--output",
    default="results/alphafold_master_table.xlsx",
    help="Output Excel file"
)

args = parser.parse_args()

PREDICTIONS_DIR = args.input_dir
OUTPUT_FILE = args.output


# =====================================================
# RUN P2RANK
# =====================================================

# Create .ds file containing all cleaned structures
os.system(
    'find data/alphafold_clean -name "*.pdb" > alphafold.ds'
)

# Run P2Rank
status = os.system(
    '~/LBS_prediction_project/p2rank_2.5/prank predict '
    'alphafold.ds '
    '-o results/alphafold '
    '-c alphafold'
)

if status == 0:
    print("P2Rank prediction completed")
else:
    print("P2Rank failed")
    exit()


# =====================================================
# FIND PREDICTION FILES
# =====================================================

files = [
    f
    for f in os.listdir(PREDICTIONS_DIR)
    if f.endswith("_predictions.csv")
]

print(f"\nFound {len(files)} prediction files")


# =====================================================
# STORE ALL POCKET DATA
# =====================================================

all_data = []


# =====================================================
# PROCESS EACH PREDICTION FILE
# =====================================================

for file in files:

    filepath = os.path.join(
        PREDICTIONS_DIR,
        file
    )

    try:

        # Read P2Rank output
        df = pd.read_csv(filepath)

        # Remove accidental spaces in column names
        df.columns = df.columns.str.strip()

        # Example:
        # AF-P00533-F1_clean.pdb_predictions.csv
        # -> AF-P00533-F1
        kinase_name = file.replace(
            "_clean.pdb_predictions.csv",
            ""
        )

        # Loop through pockets
        for _, row in df.iterrows():

            all_data.append({

                "kinase": kinase_name,

                # Pocket rank
                "pocket_rank": row.get(
                    "rank",
                    None
                ),

                # P2Rank score
                "score": row.get(
                    "score",
                    None
                ),

                # Probability already ranges from 0–1
                "probability": row.get(
                    "probability",
                    None
                ),

                # Pocket center coordinates
                "center_x": row.get(
                    "center_x",
                    None
                ),

                "center_y": row.get(
                    "center_y",
                    None
                ),

                "center_z": row.get(
                    "center_z",
                    None
                ),

                # Pocket size descriptors
                "num_surface_atoms": row.get(
                    "surf_atoms",
                    None
                ),

                "sas_points": row.get(
                    "sas_points",
                    None
                )
            })

    except Exception as e:

        print(
            f"Error reading {file}: {e}"
        )


# =====================================================
# CREATE MASTER DATAFRAME
# =====================================================

master_df = pd.DataFrame(all_data)


# Convert numeric columns
numeric_cols = [
    "pocket_rank",
    "score",
    "probability",
    "center_x",
    "center_y",
    "center_z",
    "num_surface_atoms",
    "sas_points"
]

for col in numeric_cols:

    if col in master_df.columns:

        master_df[col] = pd.to_numeric(
            master_df[col],
            errors="coerce"
        )


# Sort by kinase then pocket rank
master_df = master_df.sort_values(
    by=[
        "kinase",
        "pocket_rank"
    ]
)


# Round values for cleaner Excel output
master_df["score"] = master_df["score"].round(3)
master_df["probability"] = master_df["probability"].round(3)

master_df["center_x"] = master_df["center_x"].round(3)
master_df["center_y"] = master_df["center_y"].round(3)
master_df["center_z"] = master_df["center_z"].round(3)


# =====================================================
# SAVE MASTER TABLE
# =====================================================

master_df.to_excel(
    OUTPUT_FILE,
    index=False
)

print(
    f"\nMaster table saved to:\n{OUTPUT_FILE}"
)

print(
    f"Total predicted pockets: {len(master_df)}"
)