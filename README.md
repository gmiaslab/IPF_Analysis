# IPF Integrated Transcriptomics Analysis

## Overview

This repository contains the analysis code, harmonized inputs, and generated outputs associated with the manuscript Integrated RNA sequencing reanalysis reveals reproducible matrix-immune signatures in idiopathic
pulmonary fibrosis, by Nandimandalam et al. 

The study reanalyzed publicly available bulk RNA-seq datasets using a common processing and statistical framework to identify disease-associated transcriptional programs that are reproducible across independent cohorts while accounting for study effects, repeated sampling, age, and sex.

All manuscript figures, differential-expression results, pathway analyses, sensitivity analyses, exploratory network analyses, and classifier outputs can be regenerated from the harmonized count matrix and metadata included in the repository. *Note*: The downstream analysis is self-contained. Upstream raw-read processing steps (SRA retrieval, FASTQ preprocessing, STAR mapping, strandedness assessment, and per-study matrix assembly) are documented with example HPCC scripts in `hpcc_example_scripts/`.


## Scientific Workflow

The workflow is summarized below:

```text
Public SRA datasets
        ↓
FASTQ generation
        ↓
Adapter/read-position trimming
        ↓
STAR alignment and gene counting
        ↓
Study-level count matrices
        ↓
Harmonized count matrix + metadata
        ↓
Integrative differential-expression analysis
        ↓
Sensitivity analyses
        ↓
Classifier and imbalance analyses
        ↓
Manuscript figures, tables, and outputs
```


## Repository Contents

Documentation:

- `software_versions_used.md`: list of software and libraries used at runtime.
- `Study_Selection_2026_06_05.docx`: study-selection workflow diagram documenting GEO record screening and study inclusion.

Main data inputs:

- `combined_metadata.csv`: harmonized sample-level metadata used by the analysis models.
- `combined_study_matrix.csv`: harmonized gene-by-sample count matrix used by the analysis models.
- `ensembl_to_gene.csv`: gene annotation table used to attach gene symbols and Entrez identifiers.

Analysis scripts:

- `ipf_integrative_analysis.R`: main limma-voom differential-expression, GSEA, QC, visualization, and exploratory network workflow.
- `ipf_sensitivity_analysis.R`: sensitivity analyses, including leave-one-study-out and no-age/excluded-study checks.
- `make_sensitivity_figure.R`: summary figure helper for sensitivity-analysis outputs.
- `make_network_figure.R`: exploratory network figure helper using the network tables from the main analysis.
- `make_overlap_upset.R`: sensitivity-overlap UpSet figure helper.
- `ipf_classifier_analysis.py`: fixed-signature elastic-net leave-one-study-out IPF classifier analysis.
- `study_imbalance_analysis.py`: study-composition and imbalance diagnostics focused on GSE213001.

Raw-read processing HPCC/SLURM example scripts:

- `hpcc_example_scripts/eg_SRA_run.sb`: example SLURM script for SRA read retrieval.
- `hpcc_example_scripts/eg_trimMap.sb`: example SLURM submission wrapper for trimming and STAR mapping jobs.
- `hpcc_example_scripts/eg_TrimSTAR.sh`: example trimming and STAR mapping script.
- `hpcc_example_scripts/eg_strandedness_count_selection.R`: STAR strandedness/count-column selection helper from the upstream assembly workflow.

These example scripts document the scientific sequence of study raw data acquisition and processing steps used prior to the downstream analysis: raw data from SRA retrieval, FASTQ generation, adapter/read-position trimming, STAR mapping with gene-count output, and strandedness determination. Users may need to adapt module names, reference paths, accession-list names, read-layout settings, trimming parameters, resource requests, and scheduler/email settings for use in their own computing environment.

Generated analysis outputs:

- `ipf_analysis_results/`: main differential-expression tables, GSEA tables, QC summaries, PCA PDFs, volcano/heatmap/pathway PDFs, and exploratory network tables and figures.
- `ipf_sensitivity_results/`: sensitivity-analysis DE tables, metadata subsets, overlap summaries, and summary PDF figure.
- `ipf_classifier_results/`: fixed-signature classifier panel definitions, predictions, coefficients, tuning records, performance tables, text summary, and PDF figures.
- `study_imbalance_results/`: study-composition tables, imbalance metrics, text summary, and PDF love plot.

All outputs referenced by the manuscript are contained in these folders.

## Downstream Reproduction

The downstream scripts use relative paths with respect to the script location. To run these, the R and Python dependencies listed in `software_versions_used.md` must be first satisfied. The workflow can the be executed from the repository root.

The scripts were ran in the following order:

1. `Rscript ipf_integrative_analysis.R`
2. `Rscript ipf_sensitivity_analysis.R`
3. `Rscript make_sensitivity_figure.R`
4. `Rscript make_network_figure.R`
5. `Rscript make_overlap_upset.R`
6. `python ipf_classifier_analysis.py`
7. `python study_imbalance_analysis.py`

## References
If you use this code and/or results please cite the associated manuscript:
- Sneha Nandimandalam, Jin He and George I. Mias, *Integrated RNA sequencing reanalysis reveals reproducible matrix-immune signatures in idiopathic pulmonary fibrosis.* Journal TBD AND DOI.

## Contact Information
* Project Principal Investigator: Dr George I. Mias (gmias@msu.edu)
* [G Mias Lab](https://georgemias.org)