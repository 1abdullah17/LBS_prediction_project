import os
import argparse
import pandas as pd
from pathlib import Path

# =====================================================
# FIXED PROJECT ROOT (IMPORTANT)
# =====================================================

PREDICTIONS_DIR = "results/alphafold"
OUTPUT_FILE = "results/alphafold_master_table.xlsx"


# =====================================================
# COMMAND LINE ARGUMENTS
# =====================================================

parser = argparse.ArgumentParser(
    description="Run P2Rank and create a filtered master pocket table."
)

parser.add_argument(
    "--input_dir",
    default=str(PREDICTIONS_DIR),
    help="Directory containing P2Rank prediction files"
)

parser.add_argument(
    "--output",
    default=str(OUTPUT_FILE),
    help="Output Excel file"
)

args = parser.parse_args()

PREDICTIONS_DIR = Path(args.input_dir)
OUTPUT_FILE = Path(args.output)


# =====================================================
# RUN P2RANK
# =====================================================

os.system(
    f'find {"data/alphafold_clean"} -name "*.pdb" > alphafold.ds'
)

status = os.system(
    '~/LBS_prediction_project/p2rank_2.5/prank predict '
    'alphafold.ds '
    '-o results/alphafold '
    '-c alphafold'
)

if status != 0:
    print("P2Rank failed")
    exit()

print("P2Rank prediction completed")


# =====================================================
# FIND PREDICTION FILES (SAFE)
# =====================================================

files = sorted([
    f for f in os.listdir(PREDICTIONS_DIR)
    if f.endswith("_clean.pdb_predictions.csv")
])

print(f"\nFound {len(files)} prediction files")


# =====================================================
# STORE ALL DATA
# =====================================================

all_data = []

total_before_filter = 0
total_after_filter = 0


# =====================================================
# PROCESS FILES (ROBUST VERSION)
# =====================================================

for file in files:

    filepath = PREDICTIONS_DIR / file

    try:
        df = pd.read_csv(filepath)

        # FIX COLUMN SPACES
        df.columns = df.columns.str.strip()

        # ENSURE PROBABILITY IS NUMERIC
        df["probability"] = pd.to_numeric(df["probability"], errors="coerce")

        # DROP INVALID ROWS
        df = df.dropna(subset=["probability"])

        pockets_before = len(df)

        # FILTER
        df = df[df["probability"] >= 0.2]

        pockets_after = len(df)

        total_before_filter += pockets_before
        total_after_filter += pockets_after

        kinase_name = file.replace("_clean.pdb_predictions.csv", "")

        print(
            f"{kinase_name}: "
            f"{pockets_before - pockets_after} removed "
            f"({pockets_after} kept)"
        )

        # STORE DATA
        for _, row in df.iterrows():

            all_data.append({
                "kinase": kinase_name,
                "pocket_rank": row.get("rank"),
                "score": row.get("score"),
                "probability": row.get("probability"),
                "center_x": row.get("center_x"),
                "center_y": row.get("center_y"),
                "center_z": row.get("center_z"),
                "num_surface_atoms": row.get("surf_atoms"),
                "sas_points": row.get("sas_points")
            })

    except Exception as e:
        print(f"FAILED FILE: {file}")
        print(f"ERROR: {repr(e)}")


# =====================================================
# MASTER TABLE
# =====================================================

master_df = pd.DataFrame(all_data)

if not master_df.empty:

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
            master_df[col] = pd.to_numeric(master_df[col], errors="coerce")

    master_df = master_df.sort_values(
        by=["kinase", "pocket_rank"]
    )

    master_df["score"] = master_df["score"].round(3)
    master_df["probability"] = master_df["probability"].round(3)

    master_df.to_excel(OUTPUT_FILE, index=False)

    print("\nMASTER TABLE SAVED:", OUTPUT_FILE)

    print("Before filtering:", total_before_filter)
    print("After filtering:", total_after_filter)
    print("Removed:", total_before_filter - total_after_filter)

else:
    print("\nWARNING: Master table is empty!")