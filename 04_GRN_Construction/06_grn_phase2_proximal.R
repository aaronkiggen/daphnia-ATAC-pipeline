#!/usr/bin/env Rscript
# Phase 2: Integrating RNA-seq expression and Orthologue mappings to determine Directionality

library(dplyr)
library(tidyr)
# library(igraph)  # Missing in conda env
# library(ggraph)
# library(tidygraph)
# library(ggplot2)
# library(ggrepel)

BASE_DIR  <- "/scratch/leuven/354/vsc35429/ATAC/TF_analysis"
setwd(BASE_DIR)
OUTPUT_DIR <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"

cat("==== Starting Phase 2: Biological GRN Integration ====\n")

# 1. READ INPUTS
cat("Loading Phase 1 Structural Network and RNA-seq data...\n")
grn_edges_raw <- read.csv(file.path(OUTPUT_DIR, "01_Base_Motif_Network_Edges_proximal.csv"), stringsAsFactors = FALSE)
degs <- read.csv("DEGs.csv", stringsAsFactors = FALSE)

# 2. SET UP ORTHOLOGUE MAPPING
# For true expression inference, TF names (like GATA3 or cg) from the databases MUST 
# be mapped to Daphnia gene IDs (like gene-OUZ56_00377).
ORTHO_FILE <- "tf_ortholog_mapping.csv" # Replace with your actual orthologue mapping path if different

if (file.exists(ORTHO_FILE)) {
  cat("\nOrthologue mapping file found! Using it to map TF Names -> Daphnia IDs.\n")
  orthologue_map <- read.csv(ORTHO_FILE, stringsAsFactors = FALSE)
  # Ensure columns are named properly: tf_name, daphnia_gene_id
  colnames(orthologue_map) <- c("tf_name", "daphnia_gene_id")
} else {
  cat("\nWARNING: No orthologue mapping file ('tf_ortholog_mapping.csv') found!\n")
  cat("We will construct a fallback map assuming the TF name exactly matches the Daphnia gene ID.\n")
  cat("This will likely only capture native Daphnia TFs from the CIS-BP database.\n")
  orthologue_map <- data.frame(
    tf_name = unique(grn_edges_raw$tf_name),
    daphnia_gene_id = unique(grn_edges_raw$tf_name),
    stringsAsFactors = FALSE
  )
}

# 3. BUILD BIOLOGICAL GRN
cat("\nLoading Expression Background...\n")
bg_genes_raw <- readLines("/user/leuven/354/vsc35429/DATA/PhD/ATAC/C85_background_genes.txt")
bg_genes <- trimws(gsub("gene-", "", bg_genes_raw))

cat("\nIntegrating TF Expression & Assigning Directionality...\n")

# Prepare a clean DEGs list
degs_clean <- degs %>% mutate(clean_gene_id = gsub("gene-", "", gene_id)) 

biological_grn <- grn_edges_raw %>%
  # Clean up GRN gene queries (gene-OUZ56_XXXXX to OUZ56_XXXXX)
  mutate(
    gene_id = gsub("gene-", "", gene_id),
    # Default database_source if somehow missing
    database_source = if("database_source" %in% names(.)) database_source else "Unknown"
  ) %>%
  
  # Remove the broken Phase-1 target variables
  select(-any_of(c("target_lfc", "target_padj"))) %>%
  
  # Attach the true RNA-seq log2FoldChange for the TARGET GENE
  left_join(degs_clean %>% select(clean_gene_id, target_lfc = log2FoldChange, target_padj = padj), 
            by = c("gene_id" = "clean_gene_id")) %>%
  
  # Attach the physical Daphnia Gene ID to the transcription factor
  left_join(orthologue_map, by = "tf_name") %>%
  
  # Pull the RNA-seq log2FoldChange for the TRANSCRIPTION FACTOR ITSELF
  left_join(degs_clean %>% select(clean_gene_id, tf_lfc = log2FoldChange, tf_padj = padj), 
            by = c("daphnia_gene_id" = "clean_gene_id")) %>%
  
  # Determine if the TF is actually physically expressed and calculable
  mutate(
    tf_mapped       = !is.na(daphnia_gene_id),
    tf_has_exp_data = !is.na(tf_lfc),
    tf_is_expressed = daphnia_gene_id %in% bg_genes,
    tf_is_deg       = tf_has_exp_data & !is.na(tf_padj) & tf_padj < 0.05,
    
    # Calculate the Biological Edge Score using the magnitude of TF expression
    # If unmapped/no expression, fallback to structural motif score alone to avoid dropping CIS-BP
    biological_edge_score = if_else(
        tf_has_exp_data, 
        motif_score * abs(tf_lfc), 
        motif_score
    ),
    
    # Assign activation vs repression based on matched expression directions
    direction = case_when(
      !tf_has_exp_data ~ "unmapped",
      tf_lfc > 0 & target_lfc > 0 ~ "activation",
      tf_lfc > 0 & target_lfc < 0 ~ "repression",
      tf_lfc < 0 & target_lfc > 0 ~ "repression",
      tf_lfc < 0 & target_lfc < 0 ~ "activation",
      TRUE ~ "unknown"
    )
  )

# 4. FILTER FOR ONLY THE HIGH-CONFIDENCE BIOLOGICAL NETWORK
# To make clean biological claims, we should only trust:
# 1. TF touched the immediate promoter region (Tier1)
# 2. Relaxed filter: If it IS mapped to Daphnia (JASPAR), it MUST be expressed in the bg_genes array. 
#    If it is unmapped (CIS-BP), we keep it based on structural motif strength alone so we don't drop them.
clean_grn <- biological_grn %>%
  filter(tier == "Tier2_Proximal") %>%
  filter((tf_mapped == TRUE & tf_is_expressed == TRUE) | tf_mapped == FALSE) %>%
  arrange(desc(biological_edge_score))

# 4b. MASTER REGULATORS NETWORK (TFs that are DEGs themselves)
master_regulators_grn <- clean_grn %>%
  filter(tf_is_deg == TRUE)

cat("\n--- High-Confidence Biological Network Summary ---\n")
cat("Total mapped and expressed edges:      ", sum(biological_grn$tf_mapped), "\n")
cat("Proximal active edges (Tier 2):  ", nrow(clean_grn), "\n")
cat("Unique active Regulators found:        ", length(unique(clean_grn$tf_name)), "\n")
cat("MASTER REGULATORS (TF is a DEG):       ", length(unique(master_regulators_grn$tf_name)), " (representing", length(unique(master_regulators_grn$daphnia_gene_id)), "Daphnia DEG TFs)\n")

write.csv(biological_grn, file.path(OUTPUT_DIR, "02_Filtered_Bio_Network_Edges_All_proximal.csv"), row.names = FALSE)
write.csv(clean_grn, file.path(OUTPUT_DIR, "02_Filtered_Bio_Network_Edges_HighConf_proximal.csv"), row.names = FALSE)
write.csv(master_regulators_grn, file.path(OUTPUT_DIR, "02_Network_Centrality_Scores_proximal.csv"), row.names = FALSE)

# 5. SUMMARIZE ACTIVE HUBS
active_hubs <- clean_grn %>%
  group_by(tf_name, database_source, daphnia_gene_id) %>%
  summarise(
    n_tier1_targets     = n_distinct(gene_id),
    n_target_up         = n_distinct(gene_id[!is.na(target_lfc) & target_lfc > 0]),
    n_target_down       = n_distinct(gene_id[!is.na(target_lfc) & target_lfc < 0]),
    n_activation_edges  = sum(direction == "activation", na.rm=TRUE),
    n_repression_edges  = sum(direction == "repression", na.rm=TRUE),
    mean_bio_edge_score = mean(biological_edge_score, na.rm=TRUE),      homer_qval_enriched = first(homer_qval_enriched),    tf_lfc              = first(tf_lfc),
    tf_is_expressed     = first(tf_is_expressed),
    tf_is_deg           = first(tf_is_deg),
    .groups = "drop"
  ) %>%
  arrange(desc(n_tier1_targets))

write.csv(active_hubs, file.path(OUTPUT_DIR, "02_Top_Active_Hubs_proximal.csv"), row.names = FALSE)

# 6. PLOT THE CLEAN NETWORK
if(nrow(clean_grn) > 0) {
  cat("\nGraphing libraries missing from 'atac_motif_env'. Output CSVs are generated safely!\n")
  cat("Please plot the network map locally in standard RStudio using the CSV files.\n")
}
cat("\n==== Phase 2 Finished Successfully! ====\n")

