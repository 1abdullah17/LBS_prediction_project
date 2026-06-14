import pandas as pd

# =========================================================
# FILES
# =========================================================

CORRESPONDENCE_FILE = "results/pocket_correspondence.xlsx"

AF_FILE = "results/alphafold_master_annotated.xlsx"

CRYSTAL_FILE = "results/crystal_master_annotated.xlsx"

KINASE_FILE = "data/kinases_list.xlsx"

OUTPUT_FILE = "results/annotation_agreement.xlsx"

DISTANCE_THRESHOLD = 10.0

# =========================================================
# LOAD
# =========================================================

corr = pd.read_excel(CORRESPONDENCE_FILE)

af = pd.read_excel(AF_FILE)

crystal = pd.read_excel(CRYSTAL_FILE)

kinases = pd.read_excel(KINASE_FILE)

print("\nLoaded")

print("Correspondence:", len(corr))
print("AF:", len(af))
print("Crystal:", len(crystal))

# =========================================================
# KEEP SHARED POCKETS ONLY
# =========================================================

corr = corr[
    corr["Distance"] <= DISTANCE_THRESHOLD
].copy()

print(
    "\nShared pockets (<="
    f"{DISTANCE_THRESHOLD} Å):",
    len(corr)
)

# =========================================================
# BUILD AF ID MAPPING
# =========================================================

mapping = kinases[
    [
        "Kinase_name",
        "AlphaFold_id"
    ]
].drop_duplicates()

# =========================================================
# ADD AF ID TO CRYSTAL TABLE
# =========================================================

crystal = crystal.merge(
    mapping,
    on="Kinase_name",
    how="left"
)

# =========================================================
# CHECK
# =========================================================

missing = crystal["AlphaFold_id"].isna().sum()

print(
    "\nCrystal rows without AF ID:",
    missing
)

# =========================================================
# MERGE CRYSTAL ANNOTATION
# =========================================================

crystal_subset = crystal[
    [
        "AlphaFold_id",
        "pocket_rank",
        "annotation",
        "ligands",
        "ligand_type"
    ]
].copy()

crystal_subset = crystal_subset.rename(
    columns={
        "annotation":
            "Crystal_annotation",
        "ligands":
            "Crystal_ligands",
        "ligand_type":
            "Crystal_ligand_type"
    }
)

corr = corr.merge(
    crystal_subset,
    left_on=[
        "AlphaFold_id",
        "Crystal_rank"
    ],
    right_on=[
        "AlphaFold_id",
        "pocket_rank"
    ],
    how="left"
)

# =========================================================
# MERGE AF ANNOTATION
# =========================================================

af_subset = af[
    [
        "kinase",
        "pocket_rank",
        "annotation",
        "ligands",
        "ligand_type"
    ]
].copy()

af_subset = af_subset.rename(
    columns={
        "annotation":
            "AF_annotation",
        "ligands":
            "AF_ligands",
        "ligand_type":
            "AF_ligand_type"
    }
)

corr = corr.merge(
    af_subset,
    left_on=[
        "AlphaFold_id",
        "AF_rank"
    ],
    right_on=[
        "kinase",
        "pocket_rank"
    ],
    how="left"
)

# =========================================================
# AGREEMENT
# =========================================================

corr["Annotation_match"] = (
    corr["Crystal_annotation"]
    ==
    corr["AF_annotation"]
)

# =========================================================
# SUMMARY
# =========================================================

n = len(corr)

same = corr["Annotation_match"].sum()

different = n - same

agreement = (
    100 * same / n
    if n > 0 else 0
)

print("\n================================")
print("ANNOTATION AGREEMENT")
print("================================\n")

print("Shared pockets:", n)

print("Same annotation:", same)

print("Different annotation:", different)

print(
    f"Agreement: {agreement:.1f}%"
)

# =========================================================
# CONFUSION MATRIX
# =========================================================

print("\n================================")
print("CONFUSION MATRIX")
print("================================\n")

cm = pd.crosstab(
    corr["Crystal_annotation"],
    corr["AF_annotation"]
)

print(cm)

# =========================================================
# SAVE
# =========================================================

corr.to_excel(
    OUTPUT_FILE,
    index=False
)

print("\nSaved:")
print(OUTPUT_FILE)

# =========================================================
# EXAMPLES OF DISAGREEMENT
# =========================================================

print("\n================================")
print("EXAMPLES OF MISMATCHES")
print("================================\n")

mismatch = corr[
    ~corr["Annotation_match"]
]

if len(mismatch):

    print(
        mismatch[
            [
                "AlphaFold_id",
                "Crystal_rank",
                "AF_rank",
                "Distance",
                "Crystal_annotation",
                "AF_annotation"
            ]
        ]
        .head(20)
    )

else:

    print("No mismatches.")