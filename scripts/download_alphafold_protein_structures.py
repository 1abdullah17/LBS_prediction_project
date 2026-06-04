import requests
import pandas as pd
import os
from Bio.PDB import PDBParser, PDBIO, Select


# =========================================================
# DOWNLOAD ALPHAFOLD STRUCTURE
# =========================================================
def download_AF_model(alpha_fold_id):
    """
    Download AlphaFold protein structure using AlphaFold ID.
    """

    base_url = "https://alphafold.ebi.ac.uk/files/"
    filename = alpha_fold_id + "-model_v6.pdb"

    url = base_url + filename

    try:
        response = requests.get(url)

        if response.status_code == 200:
            return response.text

        else:
            return f"ERROR: {alpha_fold_id} not downloaded"

    except Exception:
        return f"REQUEST FAILED: {alpha_fold_id}"


# =========================================================
# GROUP CONSECUTIVE RESIDUES INTO REGIONS
# =========================================================
def group_regions(residue_numbers):
    """
    Example:
        [1,2,3,10,11,20]

    becomes:
        [(1,3), (10,11), (20,20)]
    """

    if len(residue_numbers) == 0:
        return []

    regions = []

    start = residue_numbers[0]
    end = residue_numbers[0]

    for r in residue_numbers[1:]:

        # consecutive residue
        if r == end + 1:
            end = r

        # new region
        else:
            regions.append((start, end))

            start = r
            end = r

    # add last region
    regions.append((start, end))

    return regions


# =========================================================
# CLASSIFY REGION TYPE
# =========================================================
def classify_region(start, end, total_residues):
    """
    Classify low-confidence regions as:

        - N-terminal
        - C-terminal
        - Internal region
    """

    n_terminal_limit = total_residues * 0.10
    c_terminal_limit = total_residues * 0.90

    # beginning of protein
    if end <= n_terminal_limit:
        return "N-terminal"

    # end of protein
    elif start >= c_terminal_limit:
        return "C-terminal"

    # middle of protein
    else:
        return "Internal region"


# =========================================================
# FILTER LOW-CONFIDENCE RESIDUES
# =========================================================
class PLDDTFilter(Select):
    """
    Remove residues with low pLDDT confidence.

    AlphaFold stores pLDDT values inside B-factor column.
    """

    def __init__(self, threshold=50):

        self.threshold = threshold

        # counters
        self.total_residues = 0
        self.kept_residues = 0

        # store all pLDDT values
        self.plddt_values = []

        # store low-confidence residue numbers
        self.low_conf_residues = []

    def accept_residue(self, residue):

        self.total_residues += 1

        # residue number
        residue_number = residue.get_id()[1]

        # collect pLDDT values from atoms
        b_factors = [atom.get_bfactor() for atom in residue]

        # skip empty residues
        if len(b_factors) == 0:
            return 0

        # average pLDDT for residue
        avg_plddt = sum(b_factors) / len(b_factors)

        # store pLDDT
        self.plddt_values.append(avg_plddt)

        # low-confidence residue
        if avg_plddt < self.threshold:

            self.low_conf_residues.append(residue_number)

            return 0

        # keep residue
        self.kept_residues += 1

        return 1


# =========================================================
# CLEAN STRUCTURE + ANALYZE CONFIDENCE
# =========================================================
def clean_structure(input_pdb, output_pdb, threshold=50):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", input_pdb)

    io = PDBIO()

    io.set_structure(structure)

    # apply filtering
    selector = PLDDTFilter(threshold)

    io.save(output_pdb, selector)

    # -----------------------------------------------------
    # BASIC STATISTICS
    # -----------------------------------------------------
    total = selector.total_residues

    kept = selector.kept_residues

    removed = total - kept

    avg_plddt = sum(selector.plddt_values) / len(selector.plddt_values)

    min_plddt = min(selector.plddt_values)

    max_plddt = max(selector.plddt_values)

    # -----------------------------------------------------
    # CONFIDENCE PERCENTAGES
    # -----------------------------------------------------
    low_conf_percent = (
        len(selector.low_conf_residues) / total * 100
    )

    medium_conf_percent = (
        sum(50 <= x < 70 for x in selector.plddt_values)
        / total
        * 100
    )

    high_conf_percent = (
        sum(x >= 70 for x in selector.plddt_values)
        / total
        * 100
    )

    # -----------------------------------------------------
    # CLASSIFY LOW-CONFIDENCE REGION TYPES
    # -----------------------------------------------------
    grouped_regions = group_regions(
        selector.low_conf_residues
    )

    region_types = []

    for start, end in grouped_regions:

        region_name = classify_region(
            start,
            end,
            total
        )

        region_types.append(
            f"{region_name} ({start}-{end})"
        )

    # Example:
    # N-terminal (1-25); Internal region (330-360)
    region_type_text = "; ".join(region_types)

    return {

        "total": total,
        "kept": kept,
        "removed": removed,

        "avg_plddt": avg_plddt,
        "min_plddt": min_plddt,
        "max_plddt": max_plddt,

        "low_conf_percent": low_conf_percent,
        "medium_conf_percent": medium_conf_percent,
        "high_conf_percent": high_conf_percent,

        # ONLY descriptive column kept
        "low_conf_region_type": region_type_text
    }

 