import pandas as pd
import numpy as np

# =====================================================
# FILES
# =====================================================

CRYSTAL_MASTER = (
    "results/crystal_pockets_in_af_coordinates.xlsx"
)

KINCORE_FILE = (
    "data/kinases_list.xlsx"
)

OUTPUT_FILE = (
    "results/crystal_master_annotated.xlsx"
)

# =====================================================
# PARAMETERS
# =====================================================

PROBABILITY_THRESHOLD = 0.20

MAX_BINDING_DISTANCE = 8.0

# =====================================================
# LOAD
# =====================================================

crystal = pd.read_excel(
    CRYSTAL_MASTER
)

kincore = pd.read_excel(
    KINCORE_FILE
)

print(
    "Crystal rows:",
    len(crystal)
)

# =====================================================
# STANDARDIZE
# =====================================================

crystal["Kinase_name"] = (
    crystal["Kinase_name"]
    .astype(str)
    .str.strip()
    .str.upper()
)

kincore["Kinase_name"] = (
    kincore["Kinase_name"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# =====================================================
# FILTER LOW-CONFIDENCE POCKETS
# =====================================================

crystal = crystal[
    crystal["probability"]
    >= PROBABILITY_THRESHOLD
].copy()

print(
    "Rows after probability filter:",
    len(crystal)
)

# =====================================================
# KEEP REQUIRED KINCORE COLUMNS
# =====================================================

kincore = kincore[
    [
        "Kinase_name",
        "Binding_site",
        "ligand_type",
        "Activity_label",
        "Spatial_label"
    ]
].drop_duplicates()

# =====================================================
# MERGE METADATA
# =====================================================

crystal = crystal.merge(
    kincore,
    on="Kinase_name",
    how="left",
    suffixes=("", "_kincore")
)

# =====================================================
# BUILD MASTER POCKET TABLE
# =====================================================

pocket_table = crystal[
    [
        "Kinase_name",
        "PDB",
        "rank",
        "score",
        "probability",
        "pocket_center_x",
        "pocket_center_y",
        "pocket_center_z"
    ]
].drop_duplicates()

pocket_table = pocket_table.rename(
    columns={
        "rank": "pocket_rank"
    }
)

# =====================================================
# NEW COLUMNS
# =====================================================

pocket_table["annotation"] = "unassigned"

pocket_table["ligands"] = ""

pocket_table["ligand_type"] = ""

pocket_table["activity_label"] = ""

pocket_table["spatial_label"] = ""

pocket_table["n_ligands"] = 0

pocket_table["distance_to_pocket"] = np.nan

pocket_table["matched_pocket_score"] = np.nan

pocket_table["matched_pocket_probability"] = np.nan

# =====================================================
# FIND OCCUPIED POCKETS
# =====================================================

assignments = []

group_cols = [
    "Kinase_name",
    "PDB",
    "ligand_center_x",
    "ligand_center_y",
    "ligand_center_z"
]

for keys, group in crystal.groupby(
    group_cols
):

    kinase = keys[0]

    pdb_id = keys[1]

    group = group.sort_values(
        "distance_to_ligand_A"
    )

    best = group.iloc[0]

    assignments.append({

        "kinase":
            kinase,

        "pdb":
            pdb_id,

        "occupied_pocket":
            int(best["rank"]),

        "distance":
            float(
                best["distance_to_ligand_A"]
            ),

        "binding_site":
            best["Binding_site"],

        "ligand":
            best["ligand_type"],

        "ligand_type":
            best["ligand_type_kincore"],

        "activity_label":
            best["Activity_label"],

        "spatial_label":
            best["Spatial_label"],

        "score":
            best["score"],

        "probability":
            best["probability"]

    })

assignments = pd.DataFrame(
    assignments
)

print(
    "Ligand assignments:",
    len(assignments)
)

# =====================================================
# ANNOTATE POCKETS
# =====================================================

for idx, pocket in pocket_table.iterrows():

    kinase = pocket["Kinase_name"]

    rank = pocket["pocket_rank"]

    matches = assignments[
        (assignments["kinase"] == kinase)
        &
        (assignments["occupied_pocket"] == rank)
    ]

    if len(matches) == 0:
        continue

    matches = matches[
        matches["distance"]
        <= MAX_BINDING_DISTANCE
    ]

    if len(matches) == 0:
        continue

    unique_sites = set(
        matches["binding_site"]
        .dropna()
        .astype(str)
    )

    if "ATP_site/Allosteric_site" in unique_sites:

        annotation = (
            "ATP_site/Allosteric_site"
        )

    elif "Allosteric_site" in unique_sites:

        annotation = (
            "Allosteric_site"
        )

    elif "ATP_site" in unique_sites:

        annotation = (
            "ATP_site"
        )

    else:

        annotation = (
            "unassigned"
        )

    pocket_table.loc[
        idx,
        "annotation"
    ] = annotation

    pocket_table.loc[
        idx,
        "ligands"
    ] = ";".join(
        sorted(
            set(
                matches["ligand"]
                .astype(str)
            )
        )
    )

    pocket_table.loc[
        idx,
        "ligand_type"
    ] = ";".join(
        sorted(
            set(
                matches["ligand_type"]
                .dropna()
                .astype(str)
            )
        )
    )

    pocket_table.loc[
        idx,
        "activity_label"
    ] = ";".join(
        sorted(
            set(
                matches["activity_label"]
                .dropna()
                .astype(str)
            )
        )
    )

    pocket_table.loc[
        idx,
        "spatial_label"
    ] = ";".join(
        sorted(
            set(
                matches["spatial_label"]
                .dropna()
                .astype(str)
            )
        )
    )

    pocket_table.loc[
        idx,
        "n_ligands"
    ] = len(matches)

    pocket_table.loc[
        idx,
        "distance_to_pocket"
    ] = round(
        matches["distance"].mean(),
        3
    )

    pocket_table.loc[
        idx,
        "matched_pocket_score"
    ] = round(
        matches["score"].mean(),
        2
    )

    pocket_table.loc[
        idx,
        "matched_pocket_probability"
    ] = round(
        matches["probability"].mean(),
        3
    )

# =====================================================
# SAVE
# =====================================================

pocket_table = pocket_table.sort_values(
    [
        "Kinase_name",
        "pocket_rank"
    ]
)

pocket_table.to_excel(
    OUTPUT_FILE,
    index=False
)

# =====================================================
# SUMMARY
# =====================================================

print("\n================================")

print(
    "Saved:",
    OUTPUT_FILE
)

print(
    "\nPocket count:",
    len(pocket_table)
)

print(
    "\nAnnotation counts:\n"
)

print(
    pocket_table["annotation"]
    .value_counts(
        dropna=False
    )
)

print(
    "\nAnnotated pockets:\n"
)

print(
    pocket_table[
        pocket_table["annotation"]
        != "unassigned"
    ][[
        "Kinase_name",
        "PDB",
        "pocket_rank",
        "annotation",
        "distance_to_pocket"
    ]]
    .head(30)
)

print("\nDone.")