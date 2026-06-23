# Software Versions

Updated: 2026-06-10
This file documents the software environment used in the analysis outputs. The versions are recorded run-time versions for the current outputs.

## Core environment

- R: `R version 4.4.2 (2024-10-31)`
- Python: `Python 3.13.5`

## R packages used in `ipf_integrative_analysis.R`

| Package | Version |
|---|---:|
| `clusterProfiler` | `4.14.4` |
| `edgeR` | `4.4.1` |
| `fgsea` | `1.32.0` |
| `ggplot2` | `3.5.1` |
| `ggrepel` | `0.9.6` |
| `limma` | `3.62.1` |
| `msigdbr` | `26.1.0` |
| `org.Hs.eg.db` | `3.20.0` |
| `pheatmap` | `1.0.12` |

## Additional R packages required by helper scripts

| Script | Required packages |
|---|---|
| `make_sensitivity_figure.R` | `ggplot2`, `grid` |
| `make_overlap_upset.R` | `UpSetR` |
| `make_network_figure.R` | `igraph`, `tidygraph`, `ggraph`, `ggplot2`, `scales` |

## Python packages used in `ipf_classifier_analysis.py` and `study_imbalance_analysis.py`

| Package | Version |
|---|---:|
| `matplotlib` | `3.10.7` |
| `numpy` | `2.3.4` |
| `pandas` | `2.3.3` |
| `scikit-learn` | `1.7.2` |