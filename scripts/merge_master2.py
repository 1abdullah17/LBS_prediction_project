import pandas as pd
import numpy as np

# =====================================================
# FILES
# =====================================================

AF_FILE = "results/alphafold_master_annotated.xlsx"

CRYSTAL_FILE = "results/crystal_master_annotated.xlsx"

OUTPUT_FILE = "results/combined_pocket_annotations.xlsx"

# =====================================================
# LOAD
# =====================================================

af = pd.read_excel(AF_FILE)

crystal = pd.read_excel(CRYSTAL_FILE)

# =====================================================
# ADD SOURCE
# =====================================================

af["Prediction_source"] = "AlphaFold"

crystal["Prediction_source"] = "Experimental"

# =====================================================
# RENAME CRYSTAL COLUMNS
# =====================================================

crystal = crystal.rename(
    columns={
        "Kinase_name": "kinase",
        "pocket_center_x": "center_x",
        "pocket_center_y": "center_y",
        "pocket_center_z": "center_z"
    }
)

# =====================================================
# ADD MISSING COLUMNS TO CRYSTAL
# =====================================================

if "num_surface_atoms" not in crystal.columns:
    crystal["num_surface_atoms"] = np.nan

if "sas_points" not in crystal.columns:
    crystal["sas_points"] = np.nan

# =====================================================
# ADD MISSING COLUMN TO AF
# =====================================================

if "PDB" not in af.columns:
    af["PDB"] = ""

# =====================================================
# COLUMN ORDER
# =====================================================

final_columns = [

    "Prediction_source",

    "kinase",
    "PDB",

    "pocket_rank",

    "score",
    "probability",

    "center_x",
    "center_y",
    "center_z",

    "num_surface_atoms",
    "sas_points",

    "annotation",

    "ligands",
    "ligand_type",
    "activity_label",
    "spatial_label",

    "n_ligands",

    "distance_to_pocket",

    "matched_pocket_score",
    "matched_pocket_probability"
]

af = af.reindex(columns=final_columns)

crystal = crystal.reindex(columns=final_columns)

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    [af, crystal],
    ignore_index=True
)

# =====================================================
# SAVE
# =====================================================

combined.to_excel(
    OUTPUT_FILE,
    index=False
)

# =====================================================
# SUMMARY
# =====================================================

print("\n================================")

print("Saved:")

print(OUTPUT_FILE)

print("\nRows:")

print(len(combined))

print("\nSource counts:\n")

print(
    combined["Prediction_source"]
    .value_counts()
)

print("\nAnnotation counts:\n")

print(
    pd.crosstab(
        combined["annotation"],
        combined["Prediction_source"]
    )
)

print("\nDone.")