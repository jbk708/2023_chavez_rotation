# 2023_chavez_rotation
Rotation project folder with Lukas Chavez - Fall 2023

## Overview

This project asks whether structural variants (SVs) in medulloblastoma change the expression of nearby genes.

1. **SV calling:** Inter-chromosomal SVs were called from Hi-C data with HiSV (`analysis/expression_analysis/data/inter_sv_data/`). HiSV calls were compared against breakfinder, and SVs were viewed in HiGlass (`analysis/jupyter_notebooks/`).
2. **RNA-seq annotation:** Normalized RSEM gene counts for 39 medulloblastoma samples were annotated with GENCODE gene coordinates (`analysis/RNAseq/`, `rnaseq_anno.ipynb`).
3. **Expression analysis:** For each sample, `analysis/expression_analysis/expression_analysis/expression.py` takes the genes near that sample's SV breakpoints. It compares their expression with the mean across the RNA-seq cohort and writes genes with outlier (elevated) expression to `data/gene_data/{SAMPLE}_genes.csv`.

## Samples

These 13 samples have both Hi-C SV calls and RNA-seq and were run through the expression pipeline. Subgroup information comes from `analysis/RNAseq/sample.info.archer.txt`.

| # | Sample | Consensus subgroup | Methylation subgroup | Subtype | Sex | Output in `data/gene_data/` |
|---|---|---|---|---|---|---|
| 1 | MB102 | SHH | SHH (adult) | SHHa | F | ✅ |
| 2 | MB106 | GR3 | II | G3a | M | — |
| 3 | MB164 | GR3 | III | G3b | M | — |
| 4 | MB174 | GR4 | VI | G4 | M | — (separate `data/MB174.csv`) |
| 5 | MB199 | GR4 | V | G4 | M | — |
| 6 | MB227 | GR4 | VI | G4 | M | — |
| 7 | MB234 | SHH | SHH (adult) | SHHa | F | ✅ |
| 8 | MB244 | SHH | SHH (infant) | SHHa | M | — |
| 9 | MB248 | GR3 | II | G3a | F | — |
| 10 | MB264 | GR4 | VII | G4 | F | ✅ |
| 11 | MB268 | SHH | SHH (infant) | SHHa | F | ✅ |
| 12 | MB274 | SHH | SHH (infant) | SHHa | M | — |
| 13 | MB277 | GR3 | III | G3b | F | ✅ |

Notes:
- MB275 has only an intra-chromosomal SV file, so it was not included. RCMB56 has no RNA-seq and was used only in the HiGlass and HiSV-vs-breakfinder notebooks.
- The background population is the 39 samples in `medullo_rnaseq_annotated.csv`.
- **Known issue:** `expression.py` passes a hardcoded `sample_name="MB174"` to `filter_outliers_by_sample`. As a result, per-sample outlier results may reflect MB174's expression instead of each sample's own. This has not been fixed yet.
