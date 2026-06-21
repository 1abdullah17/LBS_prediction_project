#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 18:55:01 2026

@author: sevval
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Dosya yolları
# =========================

annotation_file = "/Users/sevval/Desktop/kinase_project/annotate_pocket_ligand_overlap_FINAL.xlsx"
aminoacid_file = "/Users/sevval/Desktop/kinase_project/overlap_residue_aminoacids_top.xlsx"

output_excel = "/Users/sevval/Desktop/kinase_project/pocket_descriptor_analysis.xlsx"
fig_dir = "/Users/sevval/Desktop/kinase_project/pocket_descriptor_figures"

os.makedirs(fig_dir, exist_ok=True)

# =========================
# Verileri oku
# =========================

df = pd.read_excel(annotation_file)
aa_df = pd.read_excel(aminoacid_file)

# sadece top pocket / rank 1
df = df[df["rank"] == 1].copy()

# =========================
# Aminoasit grupları
# =========================

hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
charged = {"ASP", "GLU", "LYS", "ARG", "HIS"}

# =========================
# Her PDB-pocket için descriptor hesapla
# =========================

descriptor_rows = []

for _, row in df.iterrows():
    pdb = row["PDB"]
    pocket = row["pocket_name"]

    sub = aa_df[
        (aa_df["PDB"] == pdb) &
        (aa_df["pocket_name"] == pocket)
    ]

    aa_list = sub["amino_acid_3letter"].dropna().tolist()

    total_overlap_residues = len(aa_list)

    if total_overlap_residues > 0:
        hydrophobic_count = sum(aa in hydrophobic for aa in aa_list)
        charged_count = sum(aa in charged for aa in aa_list)

        hydrophobic_fraction = hydrophobic_count / total_overlap_residues
        charged_fraction = charged_count / total_overlap_residues
    else:
        hydrophobic_count = 0
        charged_count = 0
        hydrophobic_fraction = 0
        charged_fraction = 0

    descriptor_rows.append({
        "Kinase_name": row.get("Kinase_name", ""),
        "PDB": pdb,
        "pocket_name": pocket,
        "rank": row.get("rank", ""),
        "score": row.get("score", ""),
        "probability": row.get("probability", ""),
        "pocket_size_sas_points": row.get("sas_points", ""),
        "surf_atoms": row.get("surf_atoms", ""),
        "pocket_center_x": row.get("pocket_center_x", ""),
        "pocket_center_y": row.get("pocket_center_y", ""),
        "pocket_center_z": row.get("pocket_center_z", ""),
        "pocket_residue_count": row.get("pocket_residue_count", ""),
        "ligand_contact_residue_count": row.get("ligand_contact_residue_count", ""),
        "overlap_residue_count": row.get("overlap_residue_count", ""),
        "overlap_ratio_vs_ligand": row.get("overlap_ratio_vs_ligand", ""),
        "ligand_type_from_excel": row.get("ligand_type_from_excel", ""),
        "ligand_site_annotation": row.get("ligand_site_annotation", ""),
        "pocket_overlap_annotation": row.get("pocket_overlap_annotation", ""),
        "overlap_amino_acid_count": total_overlap_residues,
        "hydrophobic_count": hydrophobic_count,
        "charged_count": charged_count,
        "hydrophobic_fraction": hydrophobic_fraction,
        "charged_fraction": charged_fraction
    })

desc = pd.DataFrame(descriptor_rows)

desc.to_excel(output_excel, index=False)

print("Descriptor table created:")
print(output_excel)
print("Rows:", len(desc))

# =========================
# Grafikler
# =========================

# 1. Top 15 P2Rank score barplot
top_score = desc.sort_values("score", ascending=False).head(15)

plt.figure(figsize=(10, 6))
plt.bar(top_score["PDB"], top_score["score"])
plt.xlabel("PDB")
plt.ylabel("P2Rank score")
plt.title("Top 15 pockets by P2Rank score")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "bar_top15_p2rank_score.png"), dpi=300)
plt.close()

# 2. Pocket size vs P2Rank score
plt.figure(figsize=(7, 6))
plt.scatter(desc["pocket_size_sas_points"], desc["score"])
plt.xlabel("Pocket size / SAS points")
plt.ylabel("P2Rank score")
plt.title("Pocket size vs P2Rank score")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "scatter_pocket_size_vs_score.png"), dpi=300)
plt.close()

# 3. Hydrophobic fraction vs overlap ratio
plt.figure(figsize=(7, 6))
plt.scatter(desc["hydrophobic_fraction"], desc["overlap_ratio_vs_ligand"])
plt.xlabel("Hydrophobic residue fraction")
plt.ylabel("Overlap ratio vs ligand")
plt.title("Hydrophobicity vs ligand-site overlap")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "scatter_hydrophobic_vs_overlap.png"), dpi=300)
plt.close()

# 4. Charged fraction vs overlap ratio
plt.figure(figsize=(7, 6))
plt.scatter(desc["charged_fraction"], desc["overlap_ratio_vs_ligand"])
plt.xlabel("Charged residue fraction")
plt.ylabel("Overlap ratio vs ligand")
plt.title("Charged residue fraction vs ligand-site overlap")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "scatter_charged_vs_overlap.png"), dpi=300)
plt.close()

# 5. Ligand type'a göre ortalama overlap
summary = desc.groupby("ligand_site_annotation")["overlap_ratio_vs_ligand"].mean().reset_index()

plt.figure(figsize=(8, 5))
plt.bar(summary["ligand_site_annotation"], summary["overlap_ratio_vs_ligand"])
plt.xlabel("Ligand site annotation")
plt.ylabel("Mean overlap ratio vs ligand")
plt.title("Mean ligand-site overlap by ligand type")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "bar_mean_overlap_by_ligand_type.png"), dpi=300)
plt.close()

print("Figures saved in:")
print(fig_dir)

#coletaion between pocket_size and score 
corr = desc["pocket_size_sas_points"].corr(desc["score"])

print("Pearson r =", round(corr,3))

#P2Rank yüksek score verdiği pocketlar gerçekten ligand bağlanma bölgesini daha iyi mi yakalıyor?

plt.scatter(
    desc["score"],
    desc["overlap_ratio_vs_ligand"]
)

# 6. P2Rank score vs experimental ligand overlap

plt.figure(figsize=(7,6))

plt.scatter(
    desc["score"],
    desc["overlap_ratio_vs_ligand"]
)

plt.xlabel("P2Rank score")
plt.ylabel("Overlap ratio vs ligand")
plt.title("P2Rank score vs experimental ligand overlap")

plt.tight_layout()

plt.savefig(
    os.path.join(
        fig_dir,
        "scatter_p2rank_score_vs_overlap.png"
    ),
    dpi=300
)

plt.close()


corr = desc["score"].corr(
    desc["overlap_ratio_vs_ligand"]
)

print("\nP2Rank Score vs Overlap")
print("Pearson r =", round(corr,3))



