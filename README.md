# Daphnia ATAC-seq Pipeline

This repository contains the complete, automated workflow for analyzing ATAC-seq data in *Daphnia*, extending from raw reads to functional Master Regulator Gene Regulatory Networks (GRNs).

## Pipeline Structure

The pipeline is organized sequentially into five main steps:

### `01_Environment_Setup`
Scripts for building isolated, reproducible Conda environments used across the analysis, including tools for motif discovery and R-based functional enrichment packages (like `clusterProfiler`).

### `02_Preprocessing_and_PeakCalling`
Core ATAC-seq pre-processing steps:
- Read trimming, alignment, and shifting (accounting for Tn5 insertion).
- Peak calling using MACS2.
- Generating Consensus Peaks and evaluating TSS enrichment profiles.
- Defining functional regulatory tiers using `bedtools` (Promoter, Proximal, and Distal regions).

### `03_Motif_Discovery`
Identification of Transcription Factor binding sites within the open chromatin peaks utilizing HOMER and FIMO, mapped across the spatially tiered regulatory regions.

### `04_GRN_Construction`
R pipelines and execution wrappers to build co-expression networks and identify Master Regulators:
- **Phase 1**: Correlating TF motif presence with target gene expression expression.
- **Phase 2**: Pruning the network and identifying topological hubs.
- **Phase 3**: Identifying global Master Regulators and formatting multi-tier outputs for Cytoscape visualization.

### `05_Functional_Enrichment`
Hypergeometric GO enrichment mapping of the target genes bound by top Master Regulators, validated across multiple annotation databases (EggNOG, InterPro, Fantasia).

### `06_Misc_Scratch_Scripts`
Additional quality control, formatting, and intermediate test scripts used during pipeline structuring and testing.

## Execution
Scripts are numbered in their required execution order. Environmental setup (`01_Environment_Setup`) must be successfully completed before executing jobs in the SLURM HPC environment.
