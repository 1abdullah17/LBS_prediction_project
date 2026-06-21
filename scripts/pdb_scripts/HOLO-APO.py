#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 14:34:50 2026

@author: sevval
"""

import os
from Bio.PDB import PDBParser

pdb_dir = "/Users/sevval/Desktop/kinase_project/data/pdb_clean"

parser = PDBParser(QUIET=True)

results = []

for file in os.listdir(pdb_dir):

    if not file.endswith(".pdb"):
        continue

    pdb_path = os.path.join(pdb_dir, file)

    structure = parser.get_structure("x", pdb_path)

    ligands = set()

    for residue in structure.get_residues():

        hetflag = residue.id[0]

        if hetflag != " ":

            resname = residue.get_resname()

            if resname not in ["HOH","WAT"]:
                ligands.add(resname)

    status = "HOLO" if len(ligands) > 0 else "APO"

    results.append({
        "PDB_file": file,
        "status": status,
        "ligands": ";".join(sorted(ligands))
    })

import pandas as pd

df = pd.DataFrame(results)

df.to_excel(
    "/Users/sevval/Desktop/apo_holo.xlsx",
    index=False
)

print(df["status"].value_counts())
print("Saved.")