import requests
import pandas as pd
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select


# =========================================================
# DOWNLOAD ALPHAFOLD STRUCTURE
# =========================================================
def download_AF_model(alpha_fold_id):

    base_url = "https://alphafold.ebi.ac.uk/files/"
    filename = alpha_fold_id + "-model_v6.pdb"

    url = base_url + filename

    try:
        response = requests.get(url, timeout=30)

        if response.status_code == 200:
            return response.text

        else:
            return f"ERROR: {alpha_fold_id} not downloaded"

    except Exception:
        return f"REQUEST FAILED: {alpha_fold_id}"


# =========================================================
# GROUP CONSECUTIVE RESIDUES
# =========================================================
def group_regions(residue_numbers):

    if len(residue_numbers) == 0:
        return []

    residue_numbers = sorted(residue_numbers)

    regions = []

    start = residue_numbers[0]
    end = residue_numbers[0]

    for r in residue_numbers[1:]:

        if r == end + 1:
            end = r

        else:
            regions.append((start, end))

            start = r
            end = r

    regions.append((start, end))

    return regions


# =========================================================
# CLASSIFY REGION
# =========================================================
def classify_region(start, end, protein_length):

    n_limit = protein_length * 0.10
    c_limit = protein_length * 0.90

    if end <= n_limit:
        return "N-terminal"

    elif start >= c_limit:
        return "C-terminal"

    else:
        return "Internal region"


# =========================================================
# PDB FILTER
# =========================================================
class PLDDTFilter(Select):

    def __init__(self, threshold=50):

        self.threshold = threshold

        self.total_residues = 0
        self.kept_residues = 0

        self.plddt_values = []

        self.low_conf_residues = []

        self.all_residue_numbers = []

    def accept_residue(self, residue):

        self.total_residues += 1

        residue_number = residue.get_id()[1]

        self.all_residue_numbers.append(
            residue_number
        )

        b_factors = [
            atom.get_bfactor()
            for atom in residue
            if atom.get_bfactor() is not None
            and 0 <= atom.get_bfactor() <= 100
        ]

        if len(b_factors) == 0:
            return 0

        avg_plddt = np.mean(b_factors)

        self.plddt_values.append(avg_plddt)

        if avg_plddt < self.threshold:

            self.low_conf_residues.append(
                residue_number
            )

            return 0

        self.kept_residues += 1

        return 1


# =========================================================
# CLEAN STRUCTURE
# =========================================================
def clean_structure(
    input_pdb,
    output_pdb,
    threshold=50
):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        input_pdb
    )

    io = PDBIO()

    io.set_structure(structure)

    selector = PLDDTFilter(threshold)

    io.save(
        output_pdb,
        selector
    )

    if len(selector.plddt_values) == 0:
        return None

    plddt_array = np.array(
        selector.plddt_values
    )

    total = selector.total_residues

    kept = selector.kept_residues

    removed = total - kept

    avg_plddt = float(
        np.mean(plddt_array)
    )

    min_plddt = float(
        np.min(plddt_array)
    )

    max_plddt = float(
        np.max(plddt_array)
    )

    low_conf_percent = float(
        np.sum(plddt_array < 50)
        / len(plddt_array)
        * 100
    )

    medium_conf_percent = float(
        np.sum(
            (plddt_array >= 50)
            &
            (plddt_array < 70)
        )
        / len(plddt_array)
        * 100
    )

    high_conf_percent = float(
        np.sum(plddt_array >= 70)
        / len(plddt_array)
        * 100
    )

    protein_length = max(
        selector.all_residue_numbers
    )

    grouped_regions = group_regions(
        selector.low_conf_residues
    )

    region_types = []

    for start, end in grouped_regions:

        region_name = classify_region(
            start,
            end,
            protein_length
        )

        region_types.append(
            f"{region_name} ({start}-{end})"
        )

    region_text = "; ".join(
        region_types
    )

    return {

        "total": total,
        "kept": kept,
        "removed": removed,

        "avg_plddt": avg_plddt,
        "min_plddt": min_plddt,
        "max_plddt": max_plddt,

        "low_conf_percent":
            low_conf_percent,

        "medium_conf_percent":
            medium_conf_percent,

        "high_conf_percent":
            high_conf_percent,

        "low_conf_region_type":
            region_text
    }


# =========================================================
# MAIN
# =========================================================
def main():

    kinase_table = pd.read_excel(
        "data/kinase_list.xlsx"
    )

    kinase_ids = (
        kinase_table["AlphaFold_id"]
        .dropna()
        .unique()
    )

    raw_dir = "data/alphafold_raw"
    clean_dir = "data/alphafold_clean"

    os.makedirs(
        raw_dir,
        exist_ok=True
    )

    os.makedirs(
        clean_dir,
        exist_ok=True
    )

    results = []

    failed_kinases = []

    for kinase_id in kinase_ids:

        kinase_id = str(
            kinase_id
        ).strip()

        raw_path = os.path.join(
            raw_dir,
            kinase_id + ".pdb"
        )

        clean_path = os.path.join(
            clean_dir,
            kinase_id + "_clean.pdb"
        )

        if not os.path.exists(raw_path):

            model = download_AF_model(
                kinase_id
            )

            if (
                model.startswith("ERROR")
                or
                model.startswith("REQUEST")
            ):

                print(model)

                failed_kinases.append(
                    kinase_id
                )

                continue

            with open(
                raw_path,
                "w"
            ) as file:

                file.write(model)

        try:

            stats = clean_structure(
                raw_path,
                clean_path,
                threshold=50
            )

            if stats is None:

                failed_kinases.append(
                    kinase_id
                )

                continue

            percent_removed = (
                stats["removed"]
                /
                stats["total"]
                * 100
            )

            results.append({

                "kinase":
                    kinase_id,

                "total_residues":
                    stats["total"],

                "kept_residues":
                    stats["kept"],

                "removed_residues":
                    stats["removed"],

                "percent_removed":
                    round(
                        percent_removed,
                        2
                    ),

                "avg_plddt":
                    round(
                        stats["avg_plddt"],
                        2
                    ),

                "min_plddt":
                    round(
                        stats["min_plddt"],
                        2
                    ),

                "max_plddt":
                    round(
                        stats["max_plddt"],
                        2
                    ),

                "low_conf_percent":
                    round(
                        stats["low_conf_percent"],
                        2
                    ),

                "medium_conf_percent":
                    round(
                        stats["medium_conf_percent"],
                        2
                    ),

                "high_conf_percent":
                    round(
                        stats["high_conf_percent"],
                        2
                    ),

                "low_conf_region_type":
                    stats[
                        "low_conf_region_type"
                    ]
            })

        except Exception as e:

            print(
                f"Cleaning failed for "
                f"{kinase_id}: {e}"
            )

            failed_kinases.append(
                kinase_id
            )

    results_df = pd.DataFrame(
        results
    )

    results_df.to_excel(
        "results/cleaning_summary.xlsx",
        index=False
    )

    if len(failed_kinases) > 0:

        pd.DataFrame(
            {
                "failed_kinase":
                failed_kinases
            }
        ).to_csv(
            "data/failed_downloads.csv",
            index=False
        )

    print(
        "\nCleaning summary saved:"
    )

    print(
        "results/cleaning_summary.xlsx"
    )


# =========================================================
# RUN
# =========================================================
if __name__ == "__main__":

    main()