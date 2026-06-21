import pandas as pd
import numpy as np

# =====================================================
# FILES
# =====================================================

ALPHAFOLD_MASTER = "results/alphafold_master_table.xlsx"
CRYSTAL_ANNOTATION = "results/crystal_pocket_annotation.xlsx"
KINCORE_FILE = "data/kinases_list.xlsx"

OUTPUT_FILE = "results/alphafold_master_annotated.xlsx"

# =====================================================
# LOAD
# =====================================================

af = pd.read_excel(ALPHAFOLD_MASTER)

annot = pd.read_excel(CRYSTAL_ANNOTATION)

kincore = pd.read_excel(KINCORE_FILE)

# =====================================================
# STANDARDIZE KINASE NAMES
# =====================================================

kincore = kincore.rename(
    columns={
        "Kinase_name": "kinase"
    }
)

annot["kinase"] = (
    annot["kinase"]
    .astype(str)
    .str.strip()
    .str.upper()
)

kincore["kinase"] = (
    kincore["kinase"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# =====================================================
# KEEP ONLY REQUIRED KINCORE COLUMNS
# =====================================================

kincore = kincore[
    [
        "kinase",
        "Binding_site",
        "ligand_type",
        "Activity_label",
        "Spatial_label"
    ]
].drop_duplicates()

# =====================================================
# MERGE KINCORE METADATA
# =====================================================

annot = annot.merge(
    kincore,
    on="kinase",
    how="left"
)

print(
    "AlphaFold pockets:",
    len(af)
)

print(
    "Crystal annotations:",
    len(annot)
)

# =====================================================
# IGNORE NON-BINDING MOLECULES
# =====================================================

IGNORE = {
    "ACT",
    "FMT",
    "DTT",
    "BOG",
    "NO3",
    "NI"
}

# =====================================================
# NEW COLUMNS
# =====================================================

af["annotation"] = "unassigned"

af["ligands"] = ""

af["ligand_type"] = ""

af["activity_label"] = ""

af["spatial_label"] = ""

af["n_ligands"] = 0

af["distance_to_pocket"] = np.nan

af["matched_pocket_score"] = np.nan

af["matched_pocket_probability"] = np.nan

# =====================================================
# PROCESS EACH POCKET
# =====================================================

for idx, pocket in af.iterrows():

    af_id = pocket["kinase"]

    pocket_rank = pocket["pocket_rank"]

    matches = annot[
        (annot["af_id"] == af_id)
        &
        (annot["occupied_pocket"] == pocket_rank)
    ]

    if len(matches) == 0:
        continue

    ligands = []

    distances = []

    ligand_types = []

    activity_labels = []

    spatial_labels = []

    for _, row in matches.iterrows():

        ligand = str(
            row["ligand"]
        ).strip()

        if ligand in IGNORE:
            continue

        ligands.append(
            ligand
        )

        distances.append(
            row["distance_to_pocket"]
        )

        if pd.notna(
            row["ligand_type"]
        ):
            ligand_types.append(
                str(
                    row["ligand_type"]
                )
            )

        if pd.notna(
            row["Activity_label"]
        ):
            activity_labels.append(
                str(
                    row["Activity_label"]
                )
            )

        if pd.notna(
            row["Spatial_label"]
        ):
            spatial_labels.append(
                str(
                    row["Spatial_label"]
                )
            )

    if len(ligands) == 0:
        continue

    # =================================================
    # BINDING SITE CLASSIFICATION
    # =================================================

    unique_sites = set(
        matches["Binding_site"]
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

    # =================================================
    # STORE RESULTS
    # =================================================

    af.loc[
        idx,
        "annotation"
    ] = annotation

    af.loc[
        idx,
        "ligands"
    ] = ";".join(
        sorted(
            set(ligands)
        )
    )

    af.loc[
        idx,
        "ligand_type"
    ] = ";".join(
        sorted(
            set(ligand_types)
        )
    )

    af.loc[
        idx,
        "activity_label"
    ] = ";".join(
        sorted(
            set(activity_labels)
        )
    )

    af.loc[
        idx,
        "spatial_label"
    ] = ";".join(
        sorted(
            set(spatial_labels)
        )
    )

    af.loc[
        idx,
        "n_ligands"
    ] = len(
        set(ligands)
    )

    af.loc[
        idx,
        "distance_to_pocket"
    ] = round(
        np.mean(distances),
        3
    )

    af.loc[
        idx,
        "matched_pocket_score"
    ] = round(
        matches[
            "pocket_score"
        ].mean(),
        2
    )

    af.loc[
        idx,
        "matched_pocket_probability"
    ] = round(
        matches[
            "pocket_probability"
        ].mean(),
        3
    )

# =====================================================
# SAVE
# =====================================================

af.to_excel(
    OUTPUT_FILE,
    index=False
)

# =====================================================
# SUMMARY
# =====================================================

print("\n================================")

print("Saved:")

print(
    OUTPUT_FILE
)

print(
    "\nAnnotation counts:\n"
)

print(
    af["annotation"]
    .value_counts(
        dropna=False
    )
)

print(
    "\nAnnotated pockets:\n"
)

print(
    af[
        af["annotation"]
        != "unassigned"
    ][[
        "kinase",
        "pocket_rank",
        "annotation",
        "ligands",
        "ligand_type",
        "activity_label",
        "spatial_label",
        "distance_to_pocket",
        "matched_pocket_score",
        "matched_pocket_probability"
    ]]
    .head(30)
)

print("\nDone.")