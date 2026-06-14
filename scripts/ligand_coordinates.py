import os
import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, Superimposer

# =====================================================
# PATHS
# =====================================================

KINASE_FILE = "data/kinase_list.xlsx"

PDB_DIR = "data/pdb_raw"
AF_DIR = "data/alphafold_clean"

OUTPUT_FILE = "results/aligned_ligand_centroids.xlsx"

parser = PDBParser(QUIET=True)

# =====================================================
# CA EXTRACTION
# =====================================================

def get_ca_by_resnum(chain):

    ca = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        if "CA" not in residue:
            continue

        ca[residue.id[1]] = residue["CA"]

    return ca


# =====================================================
# BUILD COMMON ATOMS
# =====================================================

def build_common_atoms(
    crystal_chain,
    af_chain
):

    crystal_ca = get_ca_by_resnum(
        crystal_chain
    )

    af_ca = get_ca_by_resnum(
        af_chain
    )

    common_residues = sorted(
        set(crystal_ca.keys())
        &
        set(af_ca.keys())
    )

    print(
        "PDB:",
        min(crystal_ca),
        max(crystal_ca)
    )

    print(
        "AF:",
        min(af_ca),
        max(af_ca)
    )

    print(
        "Common:",
        min(common_residues),
        max(common_residues)
    )

    print(
        "[INFO] Common residues:",
        len(common_residues)
    )

    crystal_atoms = [
        crystal_ca[r]
        for r in common_residues
    ]

    af_atoms = [
        af_ca[r]
        for r in common_residues
    ]

    return (
        crystal_atoms,
        af_atoms,
        len(common_residues),
        min(common_residues),
        max(common_residues)
    )


# =====================================================
# ALIGN STRUCTURES
# =====================================================

def get_transform(
    crystal_chain,
    af_chain
):

    (
        crystal_atoms,
        af_atoms,
        n_common,
        common_start,
        common_end
    ) = build_common_atoms(
        crystal_chain,
        af_chain
    )

    if n_common < 100:

        raise ValueError(
            f"Too few common residues ({n_common})"
        )

    sup = Superimposer()

    sup.set_atoms(
        crystal_atoms,
        af_atoms
    )

    rmsd = sup.rms

    print(
        "[INFO] RMSD:",
        round(rmsd, 3),
        "Å"
    )

    rotation, translation = (
        sup.rotran
    )

    return (
        rotation,
        translation,
        rmsd,
        n_common,
        common_start,
        common_end
    )


# =====================================================
# LIGAND EXTRACTION
# =====================================================

def extract_ligands(chain):

    ligands = []

    IGNORE = {

        "HOH",
        "WAT",
        "DOD",

        "SO4",
        "PO4",
        "CL",

        "NA",
        "K",

        "MG",
        "CA",
        "ZN",

        "PEG",
        "EDO",
        "GOL",
        "DMS",
        "MPD"
    }

    for residue in chain:

        hetflag = residue.id[0]

        if not str(
            hetflag
        ).startswith("H_"):
            continue

        ligand_name = (
            residue.resname.strip()
        )

        if ligand_name in IGNORE:
            continue

        coords = [
            atom.coord
            for atom in residue
        ]

        if len(coords) == 0:
            continue

        centroid = np.mean(
            coords,
            axis=0
        )

        ligands.append({

            "name":
                ligand_name,

            "chain":
                chain.id,

            "residue_id":
                residue.id[1],

            "centroid":
                centroid
        })

    return ligands


# =====================================================
# TRANSFORM POINT
# =====================================================

# =====================================================
# TRANSFORM POINT
# Crystal -> AlphaFold
# =====================================================

def transform_point(
    point,
    rotation,
    translation
):
    """
    Biopython Superimposer:

        sup.set_atoms(
            crystal_atoms,
            af_atoms
        )

    computes a transform that moves:

        AF -> Crystal

    Ligands are in Crystal coordinates,
    therefore we need the inverse:

        Crystal -> AF
    """

    rotation_inv = rotation.T

    translation_inv = (
        -translation @ rotation.T
    )

    return (
        np.dot(
            point,
            rotation_inv
        )
        + translation_inv
    )


# =====================================================
# MAIN
# =====================================================

kinases = pd.read_excel(
    KINASE_FILE
)

kinases.columns = (
    kinases.columns
    .str.strip()
)

rows = []

for _, row in kinases.iterrows():

    kinase = row["Kinase_name"]

    af_id = row["AlphaFold_id"]

    pdb_chain_info = str(
        row["Pdb_chain"]
    ).strip()

    pdb_id = pdb_chain_info[:4]

    chain_id = pdb_chain_info[-1]

    crystal_file = os.path.join(
        PDB_DIR,
        f"{pdb_id}.pdb"
    )

    af_file = os.path.join(
        AF_DIR,
        f"{af_id}_clean.pdb"
    )

    print("\n=================================")
    print("Kinase:", kinase)
    print("PDB:", pdb_id)
    print("Chain:", chain_id)
    print("AF:", af_id)

    if not os.path.exists(
        crystal_file
    ):

        print(
            "[MISSING PDB]"
        )

        continue

    if not os.path.exists(
        af_file
    ):

        print(
            "[MISSING AF]"
        )

        continue

    try:

        crystal = parser.get_structure(
            "crystal",
            crystal_file
        )

        af = parser.get_structure(
            "af",
            af_file
        )

        if chain_id not in crystal[0]:

            print(
                "[MISSING CHAIN]",
                chain_id
            )

            continue

        crystal_chain = crystal[0][chain_id]

        af_chain = list(
            af.get_chains()
        )[0]

        ligands = extract_ligands(
            crystal_chain
        )

        if len(ligands) == 0:

            print(
                "[NO LIGANDS]"
            )

            continue

        print(
            "[INFO] Ligands:",
            [
                x["name"]
                for x in ligands
            ]
        )

        (
            rotation,
            translation,
            rmsd,
            n_common,
            common_start,
            common_end

        ) = get_transform(
            crystal_chain,
            af_chain
        )

        if rmsd > 5:

            print(
                "[WARNING] High RMSD"
            )

        for ligand in ligands:

            transformed = (
                transform_point(
                    ligand["centroid"],
                    rotation,
                    translation
                )
            )

            rows.append({

                "kinase":
                    kinase,

                "af_id":
                    af_id,

                "ligand":
                    ligand["name"],

                "chain":
                    ligand["chain"],

                "residue_id":
                    ligand["residue_id"],

                "common_residues":
                    n_common,

                "common_start":
                    common_start,

                "common_end":
                    common_end,

                "rmsd":
                    rmsd,

                "rmsd_flag":
                    "HIGH"
                    if rmsd > 5
                    else "OK",

                "x":
                    transformed[0],

                "y":
                    transformed[1],

                "z":
                    transformed[2]
            })

    except Exception as e:

        print(
            "[ERROR]",
            kinase,
            str(e)
        )

# =====================================================
# SAVE
# =====================================================

df = pd.DataFrame(
    rows
)

os.makedirs(
    os.path.dirname(
        OUTPUT_FILE
    ),
    exist_ok=True
)

df.to_excel(
    OUTPUT_FILE,
    index=False
)

print("\n=================================")
print(
    "Saved:",
    OUTPUT_FILE
)

print(
    "Rows:",
    len(df)
)

if len(df) > 0:

    print("\nPreview:\n")

    print(
        df.head()
    )