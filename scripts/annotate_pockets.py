import pandas as pd
import numpy as np

LIGAND_FILE = "results/aligned_ligand_centroids.xlsx"
POCKET_FILE = "results/alphafold_master_table.xlsx"
OUTPUT_FILE = "results/crystal_pocket_annotation.xlsx"

# ==========================================
# LOAD
# ==========================================

ligands = pd.read_excel(LIGAND_FILE)
pockets = pd.read_excel(POCKET_FILE)

print("Ligands :", len(ligands))
print("Pockets :", len(pockets))

# ==========================================
# FILTER
# ==========================================

ligands = ligands[
    ligands["rmsd"] < 5
].copy()

print(
    "Good ligand alignments:",
    len(ligands)
)

# ==========================================
# CHECK IDS
# ==========================================

print("\nUnique AF IDs in ligands:")
print(ligands["af_id"].nunique())

print("\nUnique AF IDs in pockets:")
print(pockets["kinase"].nunique())

missing = sorted(
    set(ligands["af_id"]) -
    set(pockets["kinase"])
)

print("\nMissing AF IDs:")
print(missing)

# ==========================================
# ANNOTATION
# ==========================================

results = []

for _, lig in ligands.iterrows():

    af_id = lig["af_id"]

    kinase_pockets = pockets[
        pockets["kinase"] == af_id
    ]

    if len(kinase_pockets) == 0:

        print(
            f"[WARNING] No pockets found for {af_id}"
        )
        continue

    ligand_xyz = np.array([
        lig["x"],
        lig["y"],
        lig["z"]
    ])

    all_distances = []

    for _, pocket in kinase_pockets.iterrows():

        pocket_xyz = np.array([
            pocket["center_x"],
            pocket["center_y"],
            pocket["center_z"]
        ])

        dist = np.linalg.norm(
            ligand_xyz - pocket_xyz
        )

        all_distances.append(
            (
                dist,
                pocket["pocket_rank"],
                pocket["score"],
                pocket["probability"]
            )
        )

    all_distances.sort(
        key=lambda x: x[0]
    )

    best = all_distances[0]

    results.append({

        "kinase":
            lig["kinase"],

        "af_id":
            af_id,

        "ligand":
            lig["ligand"],

        "residue_id":
            lig["residue_id"],

        "rmsd":
            lig["rmsd"],

        "occupied_pocket":
            int(best[1]),

        "distance_to_pocket":
            round(best[0], 3),

        "pocket_score":
            best[2],

        "pocket_probability":
            best[3],

        "second_best_distance":
            round(
                all_distances[1][0], 3
            ) if len(all_distances) > 1 else np.nan,

        "third_best_distance":
            round(
                all_distances[2][0], 3
            ) if len(all_distances) > 2 else np.nan
    })

# ==========================================
# SAVE
# ==========================================

annot = pd.DataFrame(results)

if len(annot) == 0:

    print("\nNo annotations generated.")
    raise SystemExit

annot = annot.sort_values(
    [
        "kinase",
        "distance_to_pocket"
    ]
)

annot.to_excel(
    OUTPUT_FILE,
    index=False
)

print("\n=================================")
print("Saved:", OUTPUT_FILE)
print("Rows:", len(annot))

print("\nPreview:\n")
print(annot.head(20))

print("\nPocket usage:\n")
print(
    annot["occupied_pocket"]
    .value_counts()
    .sort_index()
)

print("\nDistance statistics:\n")

print(
    annot["distance_to_pocket"]
    .describe()
)

print("\nTop 20 closest assignments:\n")

print(
    annot.sort_values(
        "distance_to_pocket"
    )[
        [
            "kinase",
            "ligand",
            "occupied_pocket",
            "distance_to_pocket"
        ]
    ].head(20)
)

print("\nTop 20 worst assignments:\n")

print(
    annot.sort_values(
        "distance_to_pocket",
        ascending=False
    )[
        [
            "kinase",
            "ligand",
            "occupied_pocket",
            "distance_to_pocket"
        ]
    ].head(20)
)