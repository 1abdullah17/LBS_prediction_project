import os
import shutil

project_dir = "/Users/sevval/Desktop/kinase_project"

folders_to_clean = [
    "data/pdb_raw",
    "data/pdb_clean",
    "results/crystal"
]

for folder in folders_to_clean:
    path = os.path.join(project_dir, folder)
    if os.path.exists(path):
        shutil.rmtree(path)
        print("Deleted:", path)

print("Eski raw, clean ve P2Rank sonuçları silindi.")





#task1 pdb structure, cleaned 
import os
import pandas as pd
import requests
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT_DIR = "/Users/sevval/Desktop/kinase_project"
os.chdir(PROJECT_DIR)

excel_file = os.path.join(PROJECT_DIR, "kinases_list.xlsx")

raw_dir = os.path.join(PROJECT_DIR, "data", "pdb_raw")
clean_dir = os.path.join(PROJECT_DIR, "data", "pdb_clean")
dataset_file = os.path.join(PROJECT_DIR, "data", "datasets", "crystal.ds")

os.makedirs(raw_dir, exist_ok=True)
os.makedirs(clean_dir, exist_ok=True)
os.makedirs(os.path.dirname(dataset_file), exist_ok=True)

df = pd.read_excel(excel_file)

class CleanProtein(Select):
    def accept_residue(self, residue):
        # Sadece suyu çıkar, ligandları ve proteinleri tut
        if residue.get_resname() == "HOH":
            return False
        return True

parser = PDBParser(QUIET=True)
clean_paths = []

for i, row in df.iterrows():
    pdb_id = str(row["Pdb_id"]).upper().strip()

    print(f"\nProcessing {pdb_id}")

    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    raw_path = os.path.join(raw_dir, f"{pdb_id}.pdb")
    clean_path = os.path.join(clean_dir, f"{pdb_id}_clean.pdb")

    if not os.path.exists(raw_path):
        try:
            r = requests.get(pdb_url)
            r.raise_for_status()
            with open(raw_path, "wb") as f:
                f.write(r.content)
            print(f"Downloaded {pdb_id}")
        except Exception as e:
            print(f"Download failed: {pdb_id} | {e}")
            continue

    try:
        structure = parser.get_structure(pdb_id, raw_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(clean_path, CleanProtein())
        print(f"Cleaned {pdb_id}")
        clean_paths.append(clean_path)
    except Exception as e:
        print(f"Cleaning failed: {pdb_id} | {e}")

with open(dataset_file, "w") as f:
    for path in clean_paths:
        f.write(path + "\n")

print("\nDONE ✅ dataset hazır:", dataset_file)
print("Cleaned PDB count:", len(clean_paths))

#task 2 p2rank
import os
import subprocess

PROJECT_DIR = "/Users/sevval/Desktop/kinase_project"
P2RANK_DIR = "/Users/sevval/Desktop/p2rank_2.5.1"

clean_dir = os.path.join(PROJECT_DIR, "data", "pdb_clean")
dataset_file = os.path.join(PROJECT_DIR, "data", "datasets", "crystal.ds")
results_dir = os.path.join(PROJECT_DIR, "results", "crystal")

os.makedirs(os.path.dirname(dataset_file), exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# 1. Tüm temiz kristal yapılarını .ds dosyasına yaz
clean_pdbs = sorted([
    os.path.join(clean_dir, f)
    for f in os.listdir(clean_dir)
    if f.endswith("_clean.pdb")
])

with open(dataset_file, "w") as f:
    for pdb_path in clean_pdbs:
        f.write(pdb_path + "\n")

print(f"Dataset file created: {dataset_file}")
print(f"Number of crystal structures: {len(clean_pdbs)}")

# 2. P2Rank default model ile çalıştır
prank_path = os.path.join(P2RANK_DIR, "prank")

cmd = [
    prank_path,
    "predict",
    dataset_file,
    "-o",
    results_dir
]

print("Running P2Rank...")
subprocess.run(cmd, check=True)

print("P2Rank finished successfully.")
print(f"Results saved in: {results_dir}")


#task 3
import pandas as pd
import glob
import os

folder = "/Users/sevval/Desktop/kinase_project/results/crystal"

rows = []

for file in glob.glob(os.path.join(folder, "*_predictions.csv")):

    try:
        df = pd.read_csv(file, skipinitialspace=True)

        df.columns = [c.strip() for c in df.columns]

        top = df.sort_values("rank").iloc[0]

        rows.append({
            "PDB": os.path.basename(file).split("_")[0],
            "rank": top["rank"],
            "score": top["score"],
            "probability": top["probability"],
            "center_x": top["center_x"],
            "center_y": top["center_y"],
            "center_z": top["center_z"]
        })

    except Exception as e:
        print(file, e)

summary = pd.DataFrame(rows)

summary.to_excel(
    "/Users/sevval/Desktop/kinase_project/top_rank_pockets_summary.xlsx",
    index=False
)

print(summary.head())
print("Top pocket summary oluşturuldu.")



#prediction results
import os
import pandas as pd

results_dir = "/Users/sevval/Desktop/kinase_project/results/crystal"

total_files = 0
with_pocket = []
empty_files = []

for f in os.listdir(results_dir):
    if f.endswith("_predictions.csv"):
        total_files += 1
        path = os.path.join(results_dir, f)
        df = pd.read_csv(path, skipinitialspace=True)
        df.columns = df.columns.str.strip()

        pdb = f.replace("_clean.pdb_predictions.csv", "")

        if len(df) > 0:
            with_pocket.append((pdb, len(df)))
        else:
            empty_files.append(pdb)

print("Prediction file count:", total_files)
print("Proteins with at least 1 pocket:", len(with_pocket))
print("Proteins with no pocket rows:", len(empty_files))

print("\nNo-pocket proteins:")
print(empty_files)

#p2rank bulamadığı
import os
import pandas as pd

excel = "/Users/sevval/Desktop/kinase_project/kinases_list.xlsx"
results_dir = "/Users/sevval/Desktop/kinase_project/results/crystal"

df = pd.read_excel(excel)
excel_ids = set(df["Pdb_id"].astype(str).str.upper().str.strip())

pred_ids = {
    f.replace("_clean.pdb_predictions.csv", "")
    for f in os.listdir(results_dir)
    if f.endswith("_predictions.csv")
}

missing = sorted(excel_ids - pred_ids)

print("Excel PDB count:", len(excel_ids))
print("P2Rank result count:", len(pred_ids))
print("Missing count:", len(missing))
print("Missing PDBs:")
print(missing)

print("\nMissing kinase information:")
print(df[df["Pdb_id"].astype(str).str.upper().str.strip().isin(missing)])











