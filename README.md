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

#### Spatial Tiering and "Russian Doll" Subtraction
When defining functional regulatory tiers (Promoter, Proximal, Distal) and generating HOMER background regions, the pipeline uses a strict "Russian doll" subtraction method to ensure mutual exclusivity of ATAC peaks (`07_bedtools_slop_regions.sh` and `08_make_homer_background.sh`). A peak can only belong to one category:
1. **Promoter (0-2kb):** Peaks intersecting the 2kb upstream/downstream window.
2. **Proximal (2-5kb):** Peaks in the 5kb window, *minus* any peaks already assigned to the Promoter.
3. **Distal (5-20kb):** Peaks in the 20kb window, *minus* any peaks already assigned to the Proximal/Promoter groups.

**Why is this important?** 
HOMER calculates statistical motif enrichment by comparing a foreground (e.g., DEG peaks) against a background (e.g., non-DEG peaks). If spatial windows overlapped (e.g., if the Proximal window still contained Promoter sequences), the statistical analysis would suffer from "signal dilution" (promoter-specific GC-rich motifs washing out enhancer-specific signals). By utilizing strict subtraction logic (`bedtools intersect -v`) uniformly across both foreground and background generation, we ensure a perfect "apples-to-apples" comparison. This prevents false positives and isolates motifs that are uniquely driving expression in that exact spatial tier.

### `03_Motif_Discovery`
Identification of Transcription Factor binding sites within the open chromatin peaks utilizing HOMER and FIMO, mapped across the spatially tiered regulatory regions.

### `04_GRN_Construction`
R pipelines and execution wrappers to build co-expression networks and identify Master Regulators:
- **Phase 1**: Correlating TF motif presence with target gene expression expression.
- **Phase 2**: Pruning the network and identifying topological hubs.
- **Phase 3**: Identifying global Master Regulators and formatting multi-tier outputs for Cytoscape visualization.

**Statistical Transparency & Motif Rescue (`homer_qval_enriched` column):**
To ensure full statistical transparency—especially for distal enhancers subject to severe multiple testing penalties—the pipeline natively rescues motifs that fail strict False Discovery Rate (FDR) correction but still exhibit strong uncorrected significance (`pval <= 0.01`). 
- The initial `01_Enriched_Motifs_TFs_*.csv` outputs explicitly list **all** successfully rescued motifs alongside a `significance_level` column, documenting whether they passed by strict FDR (`qval < 0.05`) or uncorrected `pval <= 0.01`.
- This distinction is explicitly tracked upstream and appended as a boolean `homer_qval_enriched` (`TRUE`/`FALSE`) column definitively onto every downstream generated output table (including `01_Base_Motif_Network_Edges_*.csv` and `02_Top_Active_Hubs_*.csv`). This provides immediate clarity for researchers on which Master Regulators were fiercely enriched strictly by their structural sequence alone (`TRUE`) versus those structural family fallbacks rescued by more lenient p-values but ultimately biologically validated downstream via true RNA-seq expression changes (`FALSE`).

### `05_Functional_Enrichment`
Hypergeometric GO enrichment mapping of the target genes bound by top Master Regulators, validated across multiple annotation databases (EggNOG, InterPro, Fantasia).

### `06_Misc_Scratch_Scripts`
Additional quality control, formatting, and intermediate test scripts used during pipeline structuring and testing.

## Execution
Scripts are numbered in their required execution order. Environmental setup (`01_Environment_Setup`) must be successfully completed before executing jobs in the SLURM HPC environment.
### Collapsed Ortholog Hierarchy (`Phase 4`)
The pipeline natively searches for TFs via mammalian motif databases (e.g. `FOXI1`, `FOXM1`). Because these mammalian variants evolved from single ancestral genes, Phase 4 deduplicates the redundant models and statistically collapses them back into clean, unique *Daphnia* Master Regulator gene IDs.

## Reference Data and Database Versions

To ensure full reproducibility and clarity for publication, all reference genomes, annotation files, motif databases, and functional annotation resources used in this pipeline are listed below with their exact version and download source:

### Genomes & Annotation
- **Daphnia magna reference genome:**
  - GCA_030254905.1_UOB_LRV0_1 (GFF, GTF, FASTA)
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/254/905/GCA_030254905.1_UOB_LRV0_1/


### Motif Databases
- **HOMER Motif Database:**
  - Custom build based on JASPAR 2022 (vertebrates/invertebrates) and CIS-BP (2022 release)
  - File: `motif_databases/fixed_homer.motif`
  - Download: https://jaspar.genereg.net/downloads/ (JASPAR2022), http://cisbp.ccbr.utoronto.ca/
- **Motif-to-TF mapping:**
  - File: `motif_databases/motif_to_TF_name.txt`
  - Generated from JASPAR/CIS-BP annotation tables (see `build_motif_db.sh`)


### Software/Tools
- **HOMER:** v4.11
- **bedtools:** v2.30.0
- **MACS2:** v2.2.7.1
- **R:** v4.3.3 (with Bioconductor, dplyr, clusterProfiler, etc.)
- **FIMO (MEME Suite):** v5.4.1
- **Conda environment:** See `01_Environment_Setup/` scripts for exact package versions

---

### Motif Database File Descriptions

- **cisbp_daphnia.meme**: MEME-format motif file containing Daphnia-specific motifs from the CIS-BP database, curated for use in motif scanning and enrichment analyses.
- **cisbp_daphnia_raw**: Raw/unprocessed motif data from CIS-BP for Daphnia, possibly used as an intermediate before format conversion or redundancy filtering.
- **cisbp_daphnia_pulex.zip**: Original ZIP archive downloaded from CIS-BP containing Daphnia pulex motif data and metadata.
- **combined_nonredundant.meme**: MEME-format file with all motifs (from JASPAR, CIS-BP, etc.) merged and filtered to remove redundant motifs, for comprehensive motif discovery.
- **combined_nonredundant_homer.motif**: HOMER-format version of the nonredundant motif set, used directly by HOMER for motif enrichment.
- **combined_raw.meme**: MEME-format file with all raw motifs combined (before redundancy filtering), used for initial database construction or QC.
- **fixed_homer.motif**: Final, curated HOMER-format motif database used in the pipeline; contains nonredundant, high-confidence motifs from JASPAR2026 and CIS-BP.
- **JASPAR2026_CORE_all.meme**: MEME-format file containing all core motifs from the JASPAR 2026 release (vertebrates/invertebrates), used as a primary motif source.
- **JASPAR2026_metadata.tsv**: Metadata table for JASPAR2026 motifs, including motif IDs, TF names, species, and references; used for annotation and mapping.
- **motif_to_TF_name.txt**: Tab-delimited mapping file linking motif IDs to transcription factor gene names, generated from JASPAR/CIS-BP annotation tables.
- **redundant_motif_ids.txt**: List of motif IDs identified as redundant and removed from the final database, for transparency and reproducibility.

---

## Example Output Directory Structure

A typical output directory tree after running the full pipeline might look like:

```
ATAC_Pipeline_Output/
├── 01_Preprocessing_and_PeakCalling/
│   ├── consensus_peaks.bed
│   ├── promoter_peaks.bed
│   ├── proximal_peaks.bed
│   └── distal_peaks.bed
├── 03_Motif_Discovery/
│   ├── motif_enrichment_promoter.csv
│   ├── motif_enrichment_proximal.csv
│   ├── motif_enrichment_distal.csv
│   ├── motif_region_mapping_promoter.csv
│   └── ...
├── 04_GRN_Construction/
│   ├── 01_Base_Motif_Network_Edges_promoter.csv
│   ├── 01_Base_Motif_Network_Edges_proximal.csv
│   ├── 01_Base_Motif_Network_Edges_distal.csv
│   ├── 02_Top_Active_Hubs_promoter.csv
│   ├── 02_Top_Active_Hubs_proximal.csv
│   ├── 02_Top_Active_Hubs_distal.csv
│   ├── Cytoscape_nodes.csv
│   ├── Cytoscape_edges.csv
│   └── ...
├── 05_Functional_Enrichment/
│   ├── go_enrichment_hubs.csv
│   ├── kegg_enrichment_hubs.csv
│   └── enrichment_summary_plots.pdf
└── logs/
    ├── motif_discovery.log
    └── grn_construction.log
```

---

## Example Output Table Snippets

**Motif Enrichment Table (motif_enrichment_promoter.csv):**
| motif_id | tf_name | pval | qval | region | homer_qval_enriched |
|----------|--------|------|------|--------|---------------------|
| MA1234.1 | FOXO1  | 1e-6 | 0.01 | promoter | TRUE |
| MA5678.1 | GATA1  | 2e-4 | 0.07 | promoter | FALSE |

**GRN Edge Table (01_Base_Motif_Network_Edges_promoter.csv):**
| source_tf | target_gene | motif_id | region | confidence_score |
|-----------|-------------|----------|--------|------------------|
| FOXO1     | DmagGene1   | MA1234.1 | promoter | 0.92 |
| GATA1     | DmagGene2   | MA5678.1 | promoter | 0.85 |

**Functional Enrichment Table (go_enrichment_hubs.csv):**
| regulator | gene_set | go_term | pval | qval | description |
|-----------|----------|---------|------|------|-------------|
| FOXO1     | set1     | GO:0006355 | 1e-5 | 0.01 | regulation of transcription |

---

## Pipeline Flowchart

Below is a simplified flowchart of the main pipeline steps and data flow:

```mermaid
graph TD
    A[Raw FASTQ Files] --> B[Preprocessing & Peak Calling]
    B --> C[Spatial Tiering (Promoter/Proximal/Distal)]
    C --> D[Motif Discovery (HOMER/FIMO)]
    D --> E[Motif-to-Region Mapping]
    E --> F[GRN Construction (R)]
    F --> G[Master Regulator Identification]
    G --> H[Functional Enrichment (GO/KEGG/InterPro)]
    F --> I[Cytoscape Export]
```

---
