# Daphnia ATAC-seq Pipeline

This repository contains the complete, automated workflow for analyzing ATAC-seq data in *Daphnia*, extending from raw reads to functional Master Regulator Gene Regulatory Networks (GRNs).

## Pipeline Structure

The pipeline is organized sequentially into the following main steps. Each step is performed manually from the command line and designed specifically for Daphnia magna.

---

## Input Requirements

- **Raw FASTQ files** (paired-end ATAC-seq libraries from *Daphnia magna*)
- **Genome and annotation files:**
  - Genome assembly: GCA_030254905.1_UOB_LRV0_1 (GFF, GTF, FASTA)
    - Download: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/254/905/GCA_030254905.1_UOB_LRV0_1/
- **Motif databases** (custom, Daphnia-focused):
  - Derived from JASPAR2026 vertebrate/invertebrate (https://jaspar.genereg.net/downloads/) and CIS-BP 2022 (http://cisbp.ccbr.utoronto.ca/)
  - Processed using in-house scripts for nonredundancy and Daphnia curation.
  - Key files required:
    - motif_databases/fixed_homer.motif
    - motif_databases/combined_nonredundant.meme
    - motif_databases/motif_to_TF_name.txt
    - motif_databases/JASPAR2026_metadata.tsv
- **Protein BLAST database:**
  - Daphnia magna protein FASTA (from reference above)
- **RNA-seq DEG output**: (e.g., `DEGs.csv`) for network correlation

---

### `01_Environment_Setup`
Scripts for building isolated, reproducible Conda environments used across the analysis, including tools for motif discovery and R-based functional enrichment packages (like `clusterProfiler`).
Run these scripts before any other step to ensure all dependencies are available.

---

### `02_Preprocessing_and_PeakCalling`
Core ATAC-seq pre-processing steps for Daphnia magna:
- Read trimming, alignment (against the supplied reference), Tn5 shift correction
- Peak calling using MACS2
- Generation of Consensus Peaks
- Evaluation of TSS enrichment profiles and assignment of regulatory tiers via `bedtools`

Output: Chromatin peak BED files for subsequent motif and regulatory analysis.

---

#### Spatial Tiering and "Russian Doll" Subtraction
When defining functional regulatory tiers (Promoter, Proximal, Distal) and generating HOMER background regions, the pipeline uses a strict "Russian doll" subtraction method to ensure mutual exclusivity:
1. **Promoter (0-2kb):** Peaks intersecting the 2kb upstream/downstream window.
2. **Proximal (2-5kb):** Peaks in the 5kb window, *minus* any peaks already assigned to the Promoter.
3. **Distal (5-50kb):** Peaks in the 50kb window, *minus* any peaks already assigned to the Proximal/Promoter groups.

**Why is this important?** 
HOMER calculates statistical motif enrichment by comparing a foreground (e.g., DEG peaks) against a background (e.g., non-DEG peaks). If spatial windows overlapped (e.g., if the Proximal window still contained Promoter sequences), the statistical analysis would suffer from "signal dilution" (promoter-specific GC-rich motifs washing out enhancer-specific signals). By utilizing strict subtraction logic (`bedtools intersect -v`) uniformly across both foreground and background generation, we ensure a perfect "apples-to-apples" comparison. This prevents false positives and isolates motifs that are uniquely driving expression in that exact spatial tier.

---

### `03_Motif_Discovery`
Identification of Transcription Factor binding sites within the open chromatin peaks utilizing HOMER and FIMO (using custom Daphnia motif sets), mapped across the spatially tiered regulatory regions.
- Input: Tiered peak BED files, custom motif databases.
- Output: Table of site/region-to-motif mappings, motif enrichment summaries per tier.

---

### `04_GRN_Construction`
R scripts and wrappers to build co-expression networks, map motifs to actual Daphnia gene IDs, identify Master Regulators, and collapse orthologs/hubs. Typical order:

- **00_grn_phase0_ortholog_mapping.R** —
  - Maps motif names to Daphnia gene orthologs using Uniprot ID linkage and BLASTP.
  - Produces `tf_ortholog_mapping.csv`. **Run this first!**
- **05_grn_phase2_promoter.R**, **06_grn_phase2_proximal.R**, **07_grn_phase2_distal.R** —
  - Builds the biological gene regulatory network for each spatial tier by relating motif hits to DEG-corrected gene IDs.
  - Requires `tf_ortholog_mapping.csv` and motif mapping output.
- **08_grn_phase4_collapse.R** —
  - Collapses redundant mammalian motif models onto unique Daphnia gene IDs for master regulator identification.
- **rename_outputs.sh**, **rename_pipeline_scripts.sh** — (optional) Helpers to standardize file naming throughout phases.

#### Main Inputs/Outputs in this phase:
- Inputs: Motif region mapping results, `tf_ortholog_mapping.csv`, DEGs.csv
- Outputs: Per-tier GRN network .csv files, master regulator tables, Cytoscape node/edge tables, collapsed hub files

**Statistical Transparency & Motif Rescue (`homer_qval_enriched` column):**
To ensure full statistical transparency (especially for distal enhancers subject to strong multiple testing), the pipeline records if motifs fail standard FDR correction and tracks them for downstream reference (`homer_qval_enriched` in all output tables).

---

### `05_Functional_Enrichment`
Hypergeometric (GO/KEGG/interpro) enrichment tests for target genes bound by top Master Regulators.
- Input: Hub/GRN target gene lists from previous step
- Output: Enrichment results tables, summary plots per regulator/gene set

---

### `06_Misc_Scratch_Scripts`
Additional scripts for QC, formatting, and ad-hoc investigation.
- Example: `find_missing_mr.R` examines master regulators not mapped to DEGs.
- `run_mapping.sh`: Example SLURM job and conda activation for reproducibility

---

## Execution
Scripts are numbered and named for manual, ordered execution. You must:
1. Complete environment setup (`01_Environment_Setup`).
2. Run all preprocessing and peak calling scripts for your raw ATAC-seq data (see `02_` directory scripts; adapt paths/input as needed).
3. Proceed to motif discovery; ensure required motif DBs and genome are in place (see `03_`).
4. Move to GRN construction scripts in order (see above) — always run ortholog mapping first after motif discovery.
5. After network construction and hub collapsing, proceed to 05_Functional_Enrichment for GO analysis.
6. Use supporting scripts as needed for QC/update/renaming.

Each script may require some adaptation for path and file locations; check the script headers for expected files.

---

### Example Workflow (Manual):

```bash
# 1. Setup
cd 01_Environment_Setup
bash setup_conda_envs.sh # or your custom script here

# 2. Preprocessing & Peak Calling
cd ../02_Preprocessing_and_PeakCalling
bash run_preprocessing.sh # or your provided main script

# 3. Motif Discovery
cd ../03_Motif_Discovery
Rscript run_motif_discovery.R

# 4. GRN Construction
cd ../04_GRN_Construction
Rscript 00_grn_phase0_ortholog_mapping.R
Rscript 05_grn_phase2_promoter.R
Rscript 06_grn_phase2_proximal.R
Rscript 07_grn_phase2_distal.R
Rscript 08_grn_phase4_collapse.R

# 5. Functional Enrichment
cd ../05_Functional_Enrichment
Rscript run_enrichment.R
```

---

## Typical Output Structure
- Peak BED and motif enrichment files per regulatory tier
- Motif region/sequence mapping tables
- GRN network .csv tables (per regulatory region)
- Master regulator results / Cytoscape node and edge tables
- Functional enrichment tables and summary plots

---

## Daphnia-Specific Reference Data and Database Versions

- **Genome/annotation:** GCA_030254905.1_UOB_LRV0_1 — always download latest from NCBI link above
- **Motif db:** Built from JASPAR2026, CIS-BP 2022, subject to unique curation with Daphnia motif extraction, redundancy filtering, standardized motif naming
- All database versions, download URLs, and file hashes (where possible) are tracked in the pipeline scripts and this README.

For full reproducibility, see the `01_Environment_Setup` and `build_motif_db.sh` scripts in the repo.

---

For questions or contributions, please open an issue or contact the repository owner.
