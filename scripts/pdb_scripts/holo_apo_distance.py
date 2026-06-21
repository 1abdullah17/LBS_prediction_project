import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer

project_dir = "/Users/sevval/Desktop/kinase_project"
master_file = "/Users/sevval/Desktop/kinase_project/master_pocket_table_fixed.xlsx"
pdb_dir = "/Users/sevval/Desktop/kinase_project/data/pdb_clean"

output_excel = os.path.join(project_dir, "apo_holo_distance_comparison_FINAL.xlsx")
fig_dir = os.path.join(project_dir, "apo_holo_distance_FINAL_figures")
os.makedirs(fig_dir, exist_ok=True)

pairs = [
    {"Kinase": "KIT",   "Apo_PDB": "1T45", "Holo_PDB": "1T46"},
    {"Kinase": "ABL1",  "Apo_PDB": "6XR6", "Holo_PDB": "5MO4"},
    {"Kinase": "GSK3B", "Apo_PDB": "1H8F", "Holo_PDB": "1Q5K"},
    {"Kinase": "CDK2",  "Apo_PDB": "1HCL", "Holo_PDB": "9ESJ"},
]

parser = PDBParser(QUIET=True)

df = pd.read_excel(master_file)
df.columns = df.columns.str.strip()
df["PDB"] = df["PDB"].astype(str).str.upper().str.strip()

for col in [
    "probability", "score", "sas_points",
    "pocket_center_x", "pocket_center_y", "pocket_center_z",
    "ligand_center_x", "ligand_center_y", "ligand_center_z",
    "distance_to_ligand_A"
]:
    df[col] = pd.to_numeric(df[col], errors="coerce")


def pdb_path(pdb_id):
    return os.path.join(pdb_dir, f"{pdb_id}_clean.pdb")


def get_ca_atoms(pdb_file):
    structure = parser.get_structure("x", pdb_file)
    atoms = {}

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] == " " and "CA" in res:
                    atoms[(chain.id, res.id[1])] = res["CA"]
        break

    return atoms


def align_point_from_apo_to_holo(apo_pdb, holo_pdb, point):
    apo_file = pdb_path(apo_pdb)
    holo_file = pdb_path(holo_pdb)

    apo_atoms = get_ca_atoms(apo_file)
    holo_atoms = get_ca_atoms(holo_file)

    common = sorted(set(apo_atoms.keys()) & set(holo_atoms.keys()))

    if len(common) < 20:
        print("UYARI az ortak residue:", apo_pdb, holo_pdb, len(common))
        return point

    fixed = [holo_atoms[k] for k in common]
    moving = [apo_atoms[k] for k in common]

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    rot, tran = sup.rotran
    return np.dot(point, rot) + tran


rows = []

for pair in pairs:
    kinase = pair["Kinase"]
    apo_pdb = pair["Apo_PDB"].upper()
    holo_pdb = pair["Holo_PDB"].upper()

    apo_df = df[df["PDB"] == apo_pdb].copy()
    holo_df = df[df["PDB"] == holo_pdb].copy()

    if apo_df.empty:
        print("APO yok:", apo_pdb)
        continue
    if holo_df.empty:
        print("HOLO yok:", holo_pdb)
        continue

    # APO: probability en yüksek pocket
    apo_best = apo_df.sort_values("probability", ascending=False).iloc[0]

    # HOLO: ligand merkezine en yakın pocket
    holo_best = holo_df.sort_values("distance_to_ligand_A", ascending=True).iloc[0]

    apo_center = np.array([
        apo_best["pocket_center_x"],
        apo_best["pocket_center_y"],
        apo_best["pocket_center_z"]
    ])

    holo_center = np.array([
        holo_best["pocket_center_x"],
        holo_best["pocket_center_y"],
        holo_best["pocket_center_z"]
    ])

    ligand_center = np.array([
        holo_best["ligand_center_x"],
        holo_best["ligand_center_y"],
        holo_best["ligand_center_z"]
    ])

    # APO pocket holo yapıya hizalanıyor
    apo_center_aligned = align_point_from_apo_to_holo(
        apo_pdb,
        holo_pdb,
        apo_center
    )

    apo_distance = np.linalg.norm(apo_center_aligned - ligand_center)
    holo_distance = np.linalg.norm(holo_center - ligand_center)

    rows.append({
        "Kinase": kinase,
        "Apo_PDB": apo_pdb,
        "Holo_PDB": holo_pdb,

        "Apo_selected_pocket": apo_best["pocket_name"],
        "Holo_selected_pocket": holo_best["pocket_name"],

        "Apo_probability": apo_best["probability"],
        "Holo_probability": holo_best["probability"],

        "Apo_score": apo_best["score"],
        "Holo_score": holo_best["score"],

        "Apo_sas_points": apo_best["sas_points"],
        "Holo_sas_points": holo_best["sas_points"],

        "Apo_distance_to_holo_ligand_A": apo_distance,
        "Holo_distance_to_ligand_A": holo_distance,

        "Distance_difference_Apo_minus_Holo": apo_distance - holo_distance,

        "Holo_ligand_name": holo_best.get("ligand_name", "")
    })

out = pd.DataFrame(rows)
out.to_excel(output_excel, index=False)

print("Kaydedildi:")
print(output_excel)
print(out)

# Grafik: APO vs HOLO distance
x = np.arange(len(out))
width = 0.35

plt.figure(figsize=(8, 5))
plt.bar(x - width/2, out["Apo_distance_to_holo_ligand_A"], width, label="Apo pocket to holo ligand")
plt.bar(x + width/2, out["Holo_distance_to_ligand_A"], width, label="Holo pocket to ligand")
plt.xticks(x, out["Kinase"], rotation=45, ha="right")
plt.ylabel("Distance to ligand center (Å)")
plt.title("Apo vs Holo pocket distance to ligand center")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "apo_vs_holo_distance_to_ligand_FINAL.png"), dpi=300)
plt.close()

print("Figure kaydedildi:")
print(fig_dir)