# Machine Learning-Based Pocket Prediction and Characterization on Protein Kinase Conformational Ensembles

## Overview

Protein kinases are among the most important therapeutic targets in modern drug discovery due to their central roles in cellular signaling, growth, metabolism, and differentiation. Identification of ligand-binding pockets is a fundamental step in structure-based drug design because these regions determine ligand recognition, binding affinity, and inhibitor selectivity.

This project evaluates the performance of **P2Rank**, a machine learning-based ligand-binding pocket prediction tool, on a dataset of human protein kinases. Pocket predictions generated from **AlphaFold models** and **experimentally determined crystal structures** are compared using structural alignment, ligand information, and pocket annotation analyses.

The workflow integrates structure retrieval, pocket prediction, ligand coordinate extraction, pocket annotation, and correspondence analysis to assess whether binding-site predictions derived from AlphaFold models reproduce experimentally observed binding regions.

This work was developed as part of the **BSB 680: Current Topics in Systems Biology** course project.

---

## Project Information

**Course:** BSB 680 – Current Topics in Systems Biology

**Supervisor:** Onur Serçinoğlu

**Project Members:**

* Abdullah Uğurlu
* Şevval Dik

---

## Objectives

The primary objectives of this study are:

* Predict ligand-binding pockets in human protein kinases using P2Rank.
* Identify experimentally validated ligand-binding regions from crystal structures.
* Annotate predicted pockets according to ligand proximity.
* Compare binding-site predictions obtained from AlphaFold and crystal structures.
* Evaluate the agreement between predicted and experimentally observed binding sites.
* Assess the reliability of AlphaFold models for binding-site characterization.

---

## Dataset

The dataset consists of **50 human protein kinases** selected from the KinCore database.

For each kinase:

* A ligand-bound experimental crystal structure was obtained from the Protein Data Bank (PDB).
* A corresponding AlphaFold structure was retrieved from the AlphaFold Protein Structure Database.

Both well-characterized and relatively understudied ("dark") kinases were included.

---

## Software Requirements

### Core Software

* Python 3.11
* P2Rank 2.5
* PyMOL

### Python Packages

* pandas
* numpy
* biopython
* matplotlib

Install required packages:

```bash
pip install pandas numpy biopython matplotlib
```

---

## Workflow

The complete analysis workflow is shown below.

![Pipeline Overview](pipeline.png)

### Pipeline Description

The workflow consists of two parallel analysis branches.

The **AlphaFold branch** begins with downloading and cleaning AlphaFold kinase models. Binding pockets are predicted using P2Rank and summarized in a master pocket table. Experimental ligand coordinates obtained from crystal structures are then aligned into the AlphaFold coordinate system. The predicted pocket closest to the transformed ligand centroid is assigned as the ligand-binding pocket, and the resulting annotations are integrated into a final AlphaFold pocket dataset.

The **crystal structure branch** processes experimentally determined kinase structures. Crystal structures are downloaded, analyzed using P2Rank, and transformed into the AlphaFold reference frame to enable direct comparison with AlphaFold-derived pocket predictions. This produces an annotated crystal structure pocket dataset.

The final stage of the workflow identifies corresponding pockets between AlphaFold and crystal structures. Pocket centers are compared using Euclidean distance, and the nearest pocket pairs are assigned as corresponding pockets. Annotation labels from both datasets are then compared to evaluate agreement between AlphaFold-derived and experimentally derived binding-site predictions.

---

## Output Files

### Alphafold_master_annotated

Annotated pocket dataset generated from AlphaFold structures.

Contains:

* Pocket rank
* Pocket score
* Prediction probability
* Pocket coordinates
* Pocket descriptors
* Ligand-binding annotations

---

### crystal_master_annotated

Annotated pocket dataset generated from experimentally determined crystal structures.

Contains:

* Pocket coordinates
* Pocket descriptors
* Pocket annotations

---

### pocket_correspondence

Mapping of corresponding pockets between AlphaFold and crystal structures.

Contains:

* Kinase identifiers
* Pocket pairs
* Pocket-center distances

---

### Annotation_agreement

Final comparison dataset used to evaluate consistency between AlphaFold and crystal structure annotations.

Contains:

* Matched pocket pairs
* Annotation labels
* Agreement status
* Pocket prediction metrics

---

## Methodological Summary

### Structure Preparation

AlphaFold models were filtered using residue-level pLDDT scores. Residues with average pLDDT values below 50 were removed prior to analysis.

Crystal structures were cleaned by retaining relevant protein chains and removing non-essential molecules.

### Pocket Prediction

Binding pockets were predicted using:

* P2Rank default model for crystal structures
* P2Rank AlphaFold-specific model for AlphaFold structures

For each predicted pocket, information including pocket rank, score, probability, center coordinates, and residue composition was extracted.

### Ligand Mapping

Experimental ligand coordinates were obtained from crystal structures and transformed into the AlphaFold coordinate frame through structural alignment.

Ligand centroids were then used to identify the predicted pocket most closely associated with the experimentally observed binding site.

### Pocket Correspondence Analysis

Pocket centers from AlphaFold and crystal structures were compared after structural alignment.

For each AlphaFold pocket, the nearest crystal pocket was identified based on Euclidean distance.

The resulting correspondence table was used to compare annotation labels and evaluate binding-site recovery.

---

## Key Findings

* AlphaFold kinase models generally exhibited high structural confidence and were suitable for binding-site analysis.
* P2Rank successfully identified biologically relevant binding regions in both AlphaFold and crystal structures.
* Pocket size showed a positive relationship with P2Rank score.
* Pocket score alone was not strongly associated with agreement with experimentally observed ligand-binding sites.
* Ligand-bound (holo) structures produced pocket predictions that were substantially closer to experimentally observed ligand positions than apo structures.
* Nearly half of the pockets identified in AlphaFold and crystal structures could be matched through spatial correspondence analysis.
* Most shared pockets exhibited consistent biological annotations between AlphaFold and crystal structure datasets.

---

## Repository Structure

```text
project/
│
├── scripts/
│   ├── download_af_model.py
│   ├── parse_p2rank.py
│   ├── ligand_coordinates.py
│   ├── annotate_atp_pockets.py
│   ├── annotation_final.py
│   ├── download_pdb.py
│   ├── pdb_p2rank.py
│   ├── transforming_crystal.py
│   ├── shared_pockets.py
│   └── annotation_final2.py
│
├── data/
│
├── results/
│
├── figures/
│   └── pipeline.png
│
└── README.md
```

---

## References

1. Roskoski, R. Jr. (2024). Properties of FDA-approved small molecule protein kinase inhibitors: a 2024 update. *Pharmacological Research*, 200, 107059.

2. Anderson, B., Rosston, P., Ong, H. W., Hossain, M. A., Davis-Gilbert, Z. W., & Drewry, D. H. (2023). How many kinases are druggable? A review of our current understanding. *Biochemical Journal*, 480(16), 1331–1363.

3. Krivák, R., & Hoksza, D. (2018). P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. *Journal of Cheminformatics*, 10(1), 39.

4. Modi, V., & Dunbrack, R. L. Jr. (2022). KinCore: a web resource for structural classification of protein kinases and their inhibitors. *Nucleic Acids Research*, 50(D1), D654–D664.

5. Polák, L., Škoda, P., Riedlová, K., Krivák, R., Novotný, M., & Hoksza, D. (2025). PrankWeb 4: a modular web server for protein–ligand binding site prediction and downstream analysis. *Nucleic Acids Research*, 53(W1), W466–W471.
