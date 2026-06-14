import os
import numpy as np
import pandas as pd

from Bio.PDB import PDBParser
from Bio.PDB import Superimposer

# =====================================================
# FILES
# =====================================================

KINASE_FILE = "data/kinases_list.xlsx"

CRYSTAL_MASTER = (
    "results/master_pocket_table_no_duplicates.xlsx"
)

PDB_DIR = "data/pdb_raw"

AF_DIR = "data/alphafold_clean"

OUTPUT_FILE = (
    "results/crystal_pockets_in_af_coordinates.xlsx"
)

# =====================================================
# PDB
# =====================================================

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
# COMMON ATOMS
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

    common = sorted(
        set(crystal_ca.keys())
        &
        set(af_ca.keys())
    )

    crystal_atoms = [
        crystal_ca[r]
        for r in common
    ]

    af_atoms = [
        af_ca[r]
        for r in common
    ]

    return crystal_atoms, af_atoms

# =====================================================
# TRANSFORM
# =====================================================

def get_transform(
    crystal_chain,
    af_chain
):

    crystal_atoms, af_atoms = (
        build_common_atoms(
            crystal_chain,
            af_chain
        )
    )

    sup = Superimposer()

    sup.set_atoms(
        crystal_atoms,
        af_atoms
    )

    rotation, translation = (
        sup.rotran
    )

    return (
        rotation,
        translation,
        sup.rms
    )

# =====================================================
# CRYSTAL -> AF
# =====================================================

def transform_point(
    point,
    rotation,
    translation
):

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
# LOAD
# =====================================================

kinases = pd.read_excel(
    KINASE_FILE
)

master = pd.read_excel(
    CRYSTAL_MASTER
)

rows = []

# =====================================================
# LOOP
# =====================================================

for _, kinase_row in kinases.iterrows():

    kinase = kinase_row["Kinase_name"]

    af_id = kinase_row["AlphaFold_id"]

    pdb_id = str(
        kinase_row["Pdb_chain"]
    )[:4]

    chain_id = str(
        kinase_row["Pdb_chain"]
    )[-1]

    crystal_file = os.path.join(
        PDB_DIR,
        f"{pdb_id}.pdb"
    )

    af_file = os.path.join(
        AF_DIR,
        f"{af_id}_clean.pdb"
    )

    if not os.path.exists(
        crystal_file
    ):
        continue

    if not os.path.exists(
        af_file
    ):
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

        crystal_chain = (
            crystal[0][chain_id]
        )

        af_chain = list(
            af.get_chains()
        )[0]

        (
            rotation,
            translation,
            rmsd
        ) = get_transform(
            crystal_chain,
            af_chain
        )

        subset = master[
            master["AlphaFold_id"]
            == af_id
        ]

        for _, pocket in subset.iterrows():

            xyz = np.array([
                pocket["pocket_center_x"],
                pocket["pocket_center_y"],
                pocket["pocket_center_z"]
            ])

            xyz_af = transform_point(
                xyz,
                rotation,
                translation
            )

            row = pocket.to_dict()

            row["af_x"] = xyz_af[0]
            row["af_y"] = xyz_af[1]
            row["af_z"] = xyz_af[2]

            row["alignment_rmsd"] = rmsd

            rows.append(row)

    except Exception as e:

        print(
            af_id,
            str(e)
        )

# =====================================================
# SAVE
# =====================================================

out = pd.DataFrame(rows)

out.to_excel(
    OUTPUT_FILE,
    index=False
)

print(
    "Saved:",
    OUTPUT_FILE
)

print(
    "Rows:",
    len(out)
)