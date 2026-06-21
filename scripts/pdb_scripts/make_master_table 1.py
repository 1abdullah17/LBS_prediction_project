import argparse
import os
import glob
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


def get_ligands_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    ligands = []

    for residue in structure.get_residues():
        hetflag = residue.id[0]

        if hetflag == " ":
            continue

        resname = residue.get_resname().strip()

        if resname in ["HOH", "WAT"]:
            continue

        coords = np.array([atom.coord for atom in residue])

        if len(coords) == 0:
            continue

        center = coords.mean(axis=0)

        ligands.append({
            "ligand_name": resname,
            "ligand_center_x": center[0],
            "ligand_center_y": center[1],
            "ligand_center_z": center[2]
        })

    return ligands


def distance_3d(x1, y1, z1, x2, y2, z2):
    return float(np.sqrt(
        (x1 - x2) ** 2 +
        (y1 - y2) ** 2 +
        (z1 - z2) ** 2
    ))


def main():
    parser = argparse.ArgumentParser(
        description="Create master P2Rank table with ligand names, ligand centers, and pocket-ligand distances."
    )

    parser.add_argument("--excel", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--pdb_dir", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    kinase_df = pd.read_excel(args.excel)
    kinase_df["Pdb_id"] = kinase_df["Pdb_id"].astype(str).str.upper().str.strip()

    all_rows = []

    prediction_files = glob.glob(os.path.join(args.results_dir, "*_predictions.csv"))

    for pred_file in prediction_files:
        try:
            pred_df = pd.read_csv(pred_file, skipinitialspace=True)
            pred_df.columns = pred_df.columns.str.strip()

            if pred_df.empty:
                continue

            pdb_id = os.path.basename(pred_file).split("_clean")[0].upper()

            pdb_file = os.path.join(args.pdb_dir, f"{pdb_id}_clean.pdb")
            ligands = get_ligands_from_pdb(pdb_file)

            meta = kinase_df[kinase_df["Pdb_id"] == pdb_id]

            if len(meta) > 0:
                meta_row = meta.iloc[0]
                kinase_name = meta_row.get("Kinase_name", "")
                uniprot_id = meta_row.get("Uniprot_id", "")
                ligand_type = meta_row.get("ligand_type", "")
                alphafold_id = meta_row.get("AlphaFold_id", "")
            else:
                kinase_name = ""
                uniprot_id = ""
                ligand_type = ""
                alphafold_id = ""

            for _, row in pred_df.iterrows():
                pocket_x = row.get("center_x", np.nan)
                pocket_y = row.get("center_y", np.nan)
                pocket_z = row.get("center_z", np.nan)

                residue_ids = row.get("residue_ids", "")
                residue_count = len(str(residue_ids).replace(",", " ").split())

                if len(ligands) == 0:
                    all_rows.append({
                        "Kinase_name": kinase_name,
                        "Uniprot_id": uniprot_id,
                        "PDB": pdb_id,
                        "ligand_type_from_excel": ligand_type,
                        "AlphaFold_id": alphafold_id,
                        "pocket_name": row.get("name", ""),
                        "rank": row.get("rank", ""),
                        "score": row.get("score", ""),
                        "probability": row.get("probability", ""),
                        "sas_points": row.get("sas_points", ""),
                        "surf_atoms": row.get("surf_atoms", ""),
                        "pocket_center_x": pocket_x,
                        "pocket_center_y": pocket_y,
                        "pocket_center_z": pocket_z,
                        "residue_count": residue_count,
                        "residue_ids": residue_ids,
                        "ligand_name": "",
                        "ligand_center_x": "",
                        "ligand_center_y": "",
                        "ligand_center_z": "",
                        "distance_to_ligand_A": ""
                    })
                else:
                    for lig in ligands:
                        distance = distance_3d(
                            pocket_x, pocket_y, pocket_z,
                            lig["ligand_center_x"],
                            lig["ligand_center_y"],
                            lig["ligand_center_z"]
                        )

                        all_rows.append({
                            "Kinase_name": kinase_name,
                            "Uniprot_id": uniprot_id,
                            "PDB": pdb_id,
                            "ligand_type_from_excel": ligand_type,
                            "AlphaFold_id": alphafold_id,
                            "pocket_name": row.get("name", ""),
                            "rank": row.get("rank", ""),
                            "score": row.get("score", ""),
                            "probability": row.get("probability", ""),
                            "sas_points": row.get("sas_points", ""),
                            "surf_atoms": row.get("surf_atoms", ""),
                            "pocket_center_x": pocket_x,
                            "pocket_center_y": pocket_y,
                            "pocket_center_z": pocket_z,
                            "residue_count": residue_count,
                            "residue_ids": residue_ids,
                            "ligand_name": lig["ligand_name"],
                            "ligand_center_x": lig["ligand_center_x"],
                            "ligand_center_y": lig["ligand_center_y"],
                            "ligand_center_z": lig["ligand_center_z"],
                            "distance_to_ligand_A": distance
                        })

        except Exception as e:
            print(f"Skipped {pred_file}: {e}")

    master_df = pd.DataFrame(all_rows)

    master_df = master_df.sort_values(
        by=["PDB", "rank", "distance_to_ligand_A"],
        ascending=[True, True, True]
    )

    master_df.to_excel(args.output, index=False)

    print("Master table with ligands created:")
    print(args.output)
    print("Total rows:", len(master_df))


if __name__ == "__main__":
    main()
    


    
    
    
    
    
    
    