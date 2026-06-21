import pandas as pd
import numpy as np

# =====================================================
# FILES
# =====================================================

AF_FILE = (
    "results/alphafold_master_annotated.xlsx"
)

CRYSTAL_FILE = (
    "results/crystal_pockets_in_af_coordinates.xlsx"
)

OUTPUT_FILE = (
    "results/pocket_correspondence.xlsx"
)

# =====================================================
# LOAD
# =====================================================

af = pd.read_excel(AF_FILE)

crystal = pd.read_excel(CRYSTAL_FILE)

# =====================================================
# KEEP ONLY GOOD POCKETS
# =====================================================

af = af[
    af["probability"] > 0.2
]

crystal = crystal[
    crystal["probability"] > 0.2
]

# =====================================================
# MATCH
# =====================================================

rows = []

for af_id in sorted(
    set(af["kinase"])
):

    af_pockets = af[
        af["kinase"]
        == af_id
    ]

    crystal_pockets = crystal[
        crystal["AlphaFold_id"]
        == af_id
    ]

    if len(crystal_pockets) == 0:
        continue

    for _, cp in crystal_pockets.iterrows():

        c_xyz = np.array([
            cp["af_x"],
            cp["af_y"],
            cp["af_z"]
        ])

        distances = []

        for _, ap in af_pockets.iterrows():

            a_xyz = np.array([
                ap["center_x"],
                ap["center_y"],
                ap["center_z"]
            ])

            d = np.linalg.norm(
                c_xyz - a_xyz
            )

            distances.append(
                (
                    d,
                    ap["pocket_rank"],
                    ap["score"],
                    ap["probability"]
                )
            )

        distances.sort(
            key=lambda x: x[0]
        )

        best = distances[0]

        rows.append({

            "AlphaFold_id":
                af_id,

            "Crystal_rank":
                cp["rank"],

            "AF_rank":
                best[1],

            "Distance":
                best[0],

            "Crystal_score":
                cp["score"],

            "AF_score":
                best[2],

            "Crystal_probability":
                cp["probability"],

            "AF_probability":
                best[3]
        })

# =====================================================
# SAVE
# =====================================================

results = pd.DataFrame(rows)

results["Pocket_status"] = np.where(
    results["Distance"] < 10,
    "Shared",
    "Crystal_only"
)

results.to_excel(
    OUTPUT_FILE,
    index=False
)

print(
    results["Pocket_status"]
    .value_counts()
)

# =====================================================
# SUMMARY
# =====================================================

print("\nDistance stats\n")

print(
    results["Distance"]
    .describe()
)

for cutoff in [
    3,
    5,
    8,
    10,
    15,
    20
]:

    n = (
        results["Distance"]
        < cutoff
    ).sum()

    print(
        f"{cutoff:>2} Å : {n}"
    )