#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:31:54 2026

@author: sevval
"""

import os
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
import re


def parse_p2rank_residues(residue_string):
    if pd.isna(residue_string):
        return set()

    residues = str(residue_string).replace(",", " ").split()
    return set(residues)


from Bio.PDB import PDBParser, NeighborSearch
import numpy as np

def get_ligand_contact_residues(pdb_file, cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    protein_atoms = []
    ligand_atoms = []

    for atom in structure.get_atoms():
        residue = atom.get_parent()
        hetflag = residue.id[0]
        resname = residue.get_resname().strip()

        if hetflag == " ":
            protein_atoms.append(atom)
        elif resname not in ["HOH", "WAT"]:
            ligand_atoms.append(atom)

    ns = NeighborSearch(protein_atoms)

    contact_residues = set()
    ligand_names = set()

    for lig_atom in ligand_atoms:
        lig_res = lig_atom.get_parent()
        ligand_names.add(lig_res.get_resname().strip())

        nearby_atoms = ns.search(lig_atom.coord, cutoff)

        for atom in nearby_atoms:
            prot_res = atom.get_parent()
            chain_id = prot_res.get_parent().id
            resi = prot_res.id[1]
            contact_residues.add(f"{chain_id}_{resi}")

    return contact_residues, sorted(ligand_names)



def classify_ligand_type(text):
    text = str(text).lower()
    text = text.replace("/", " ").replace("-", " ").replace("_", " ")

    if "allosteric" in text:
        return "Allosteric"

    if "atplike" in text or "atp like" in text:
        return "ATP-like"

    if "type i" in text:
        return "Type I"

    if "type ii" in text:
        return "Type II"

    return "Unclassified"

def classify_overlap(ratio):
    if ratio >= 0.50:
        return "Strong ligand-site match"
    elif ratio >= 0.20:
        return "Partial ligand-site match"
    else:
        return "Weak/no ligand-site match"

def main():
    parser = argparse.ArgumentParser(
        description="Compare P2Rank pocket residues with ligand-contact residues."
    )

    parser.add_argument("--master_table", required=True)
    parser.add_argument("--pdb_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cutoff", type=float, default=5.0)

    args = parser.parse_args()

    df = pd.read_excel(args.master_table)

    rows = []

    for _, row in df.iterrows():
        pdb_id = str(row["PDB"]).upper().strip()

        pdb_file = os.path.join(args.pdb_dir, f"{pdb_id}_clean.pdb")

        if not os.path.exists(pdb_file):
            print(f"Missing PDB file: {pdb_id}")
            continue

        pocket_residues = parse_p2rank_residues(row["residue_ids"])

        ligand_contact_residues, ligand_names = get_ligand_contact_residues(
            pdb_file,
            cutoff=args.cutoff
        )

        overlap = pocket_residues.intersection(ligand_contact_residues)

        if len(pocket_residues) > 0:
            overlap_ratio_pocket = len(overlap) / len(pocket_residues)
        else:
            overlap_ratio_pocket = 0

        if len(ligand_contact_residues) > 0:
            overlap_ratio_ligand = len(overlap) / len(ligand_contact_residues)
        else:
            overlap_ratio_ligand = 0

        ligand_type = row.get("ligand_type_from_excel", "")
        annotation = classify_ligand_type(ligand_type)

        new_row = row.to_dict()
      

        new_row.update({
           "ligand_names_in_pdb": ";".join(ligand_names),
           "ligand_contact_residue_count": len(ligand_contact_residues),
           "ligand_contact_residues": " ".join(sorted(ligand_contact_residues)),
           "pocket_residue_count": len(pocket_residues),
           "overlap_residue_count": len(overlap),
           "overlap_residues": " ".join(sorted(overlap)),
           "overlap_ratio_vs_pocket": overlap_ratio_pocket,
           "overlap_ratio_vs_ligand": overlap_ratio_ligand,
           "ligand_site_annotation": annotation,
           "pocket_overlap_annotation": classify_overlap(overlap_ratio_ligand)
})

        rows.append(new_row)

    out = pd.DataFrame(rows)

    out.to_excel(args.output, index=False)

    print("Annotation table created:")
    print(args.output)
    print("Rows:", len(out))


if __name__ == "__main__":
    main()
    
     
    
#%run /Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap.py --master_table "/Users/sevval/Desktop/kinase_project/master_pocket_table_no_duplicates.xlsx" --pdb_dir "/Users/sevval/Desktop/kinase_project/data/pdb_clean" --output "/Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap.xlsx"


import pandas as pd

df = pd.read_excel(
"/Users/sevval/Desktop/kinase_project/master_pocket_table_no_duplicates.xlsx"
)

print(df.columns.tolist())


#SONRdan overlap kısmı eklendi

import pandas as pd

file = "/Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap.xlsx"

df = pd.read_excel(file)

def classify_ligand_type(text):
    text = str(text).lower()
    text = text.replace("/", " ").replace("-", " ").replace("_", " ")

    if "allosteric" in text:
        return "Allosteric"
    elif "atplike" in text or "atp like" in text:
        return "ATP-like"
    elif "type ii" in text:
        return "Type II"
    elif "type i" in text:
        return "Type I"
    else:
        return "Unclassified"

def classify_overlap(ratio):
    try:
        ratio = float(ratio)
    except:
        return "No overlap data"

    if ratio >= 0.50:
        return "Strong ligand-site match"
    elif ratio >= 0.20:
        return "Partial ligand-site match"
    else:
        return "Weak/no ligand-site match"

df["ligand_site_annotation"] = df["ligand_type_from_excel"].apply(classify_ligand_type)

df["pocket_overlap_annotation"] = df["overlap_ratio_vs_ligand"].apply(classify_overlap)

out = "/Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap_FINAL.xlsx"

df.to_excel(out, index=False)

print("Kaydedildi:", out)
print(df[["ligand_type_from_excel", "ligand_site_annotation", "overlap_ratio_vs_ligand", "pocket_overlap_annotation"]].head(20))








