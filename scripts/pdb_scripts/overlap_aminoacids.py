#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 13:57:50 2026

@author: sevval
"""

import os
import pandas as pd
from Bio.PDB import PDBParser

annotation_file = "/Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap_FINAL.xlsx"
pdb_dir = "/Users/sevval/Desktop/kinase_project/data/pdb_clean"
output_file = "/Users/sevval/Desktop/kinase_project/overlap_residue_aminoacids_top.xlsx"

df = pd.read_excel(annotation_file)
df = df[df["rank"] == 1]

parser = PDBParser(QUIET=True)

rows = []

for _, row in df.iterrows():
    pdb_id = str(row["PDB"]).upper().strip()
    pocket_name = row.get("pocket_name", "")
    overlap_residues = row.get("overlap_residues", "")

    if pd.isna(overlap_residues) or str(overlap_residues).strip() == "":
        continue

    pdb_file = os.path.join(pdb_dir, f"{pdb_id}_clean.pdb")

    if not os.path.exists(pdb_file):
        print("Missing PDB:", pdb_id)
        continue

    structure = parser.get_structure(pdb_id, pdb_file)

    residue_dict = {}

    for residue in structure.get_residues():
        chain_id = residue.get_parent().id
        resi = residue.id[1]
        resname = residue.get_resname().strip()

        key = f"{chain_id}_{resi}"
        residue_dict[key] = resname

    for res_id in str(overlap_residues).split():
        aa_name = residue_dict.get(res_id, "NOT_FOUND")

        rows.append({
            "PDB": pdb_id,
            "pocket_name": pocket_name,
            "overlap_residue_id": res_id,
            "amino_acid_3letter": aa_name,
            "ligand_type_from_excel": row.get("ligand_type_from_excel", ""),
            "ligand_site_annotation": row.get("ligand_site_annotation", ""),
            "pocket_overlap_annotation": row.get("pocket_overlap_annotation", "")
        })

out = pd.DataFrame(rows)

out.to_excel(output_file, index=False)

print("Overlap residue amino acid table created:")
print(output_file)
print("Rows:", len(out))