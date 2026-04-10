#!/usr/bin/env Rscript
# Phase 3: Paper-level Visualization & Cytoscape Export (Multi-Tier)

library(dplyr)
library(tidyr)
library(ggplot2)

BASE_DIR  <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output"
setwd(BASE_DIR)

cat("==== Starting Phase 3: Paper Visualizations & Cytoscape Export ====\n")

tiers <- c("promoter", "proximal", "distal")

for (tier in tiers) {
  cat(sprintf("\n--- Processing %s tier ---\n", toupper(tier)))
  
  # 1. LOAD PHASE 2 OUTPUTS
  # Check if the distal file exists without the suffix or with it
  if(tier == "distal") {
     hubs_file <- "top_TF_Biological_Hubs.csv"
     edges_file <- "GRN_Biological_Edges_Tier1_Only.csv"
  } else {
     hubs_file <- paste0("top_TF_Biological_Hubs_", tier, ".csv")
     edges_file <- paste0("GRN_Biological_Edges_Tier1_Only_", tier, ".csv")
  }
  
  if (!file.exists(hubs_file) | !file.exists(edges_file)) {
    cat(sprintf("Missing Phase 2 outputs for %s. Skipping...\n", tier))
    next
  }
  
  hubs  <- read.csv(hubs_file, stringsAsFactors = FALSE)
  edges <- read.csv(edges_file, stringsAsFactors = FALSE)
  
  # ------------------------------------------------------------------------------
  # VISUAL 1: The Master Regulator Bubble Plot (ggplot2)
  # ------------------------------------------------------------------------------
  cat("Generating Master Regulator Bubble Plot...\n")
  
  mr_hubs <- hubs %>% filter(tf_is_deg == TRUE) %>% head(15) 
  if(nrow(mr_hubs) == 0) {
    mr_hubs <- hubs %>% filter(tf_is_expressed == TRUE) %>% head(15)
  }
  
  if(nrow(mr_hubs) > 0) {
    p_bubble <- ggplot(mr_hubs, aes(x = n_tier1_targets, y = reorder(tf_name, n_tier1_targets))) +
      geom_point(aes(size = mean_bio_edge_score, fill = tf_lfc), shape = 21, color="black", alpha = 0.8) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "TF Log2FC") +
      scale_size_continuous(range = c(3, 10), name = "Mean Edge Score") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        axis.text.y = element_text(face = "bold")
      ) +
      labs(
        title = sprintf("Top Master Regulators (%s)", tools::toTitleCase(tier)),
        x = "Number of Target Genes",
        y = "Transcription Factor"
      )
    
    ggsave(paste0("Fig1_Master_Regulators_Bubble_", tier, ".pdf"), p_bubble, width = 9, height = 7)
  }
  
  # ------------------------------------------------------------------------------
  # VISUAL 2: Cytoscape Export (Nodes & Edges) with GO info
  # ------------------------------------------------------------------------------
  cat("Generating Cytoscape Export Files...\n")
  
  top_tfs_for_network <- mr_hubs$tf_name[1:min(5, nrow(mr_hubs))]
  
  cyto_edges <- edges %>%
    filter(tf_name %in% top_tfs_for_network) %>%
    select(source = tf_name, target = gene_id, interaction = direction, weight = biological_edge_score) %>%
    mutate(norm_weight = weight / max(weight, na.rm=TRUE))
  
  write.csv(cyto_edges, paste0("Cytoscape_Edges_Top5_MRs_", tier, ".csv"), row.names = FALSE)
  
  sources <- edges %>%
    filter(tf_name %in% top_tfs_for_network) %>%
    select(node_id = tf_name, lfc = tf_lfc, padj = tf_padj) %>%
    mutate(node_type = "Transcription Factor", is_deg = TRUE) %>%
    distinct(node_id, .keep_all = TRUE)
  
  targets <- edges %>%
    filter(tf_name %in% top_tfs_for_network) %>%
    select(node_id = gene_id, lfc = target_lfc, padj = target_padj) %>%
    mutate(node_type = "Target Gene", is_deg = (!is.na(padj) & padj < 0.05)) %>%
    distinct(node_id, .keep_all = TRUE)
  
  cyto_nodes <- bind_rows(sources, targets) %>% group_by(node_id) %>% slice(1) %>% ungroup()
  
  # Load GO results if available to attach to nodes
  go_file <- "GO_enrichment_results/Master_Regulator_GO_Enrichment_All_Tools.csv"
  if(file.exists(go_file)) {
      go_res <- read.csv(go_file)
      # Collapse top 3 GO terms for each gene that is targeted by these TFs
      # (clusterProfiler output lists geneID as string 'gene1/gene2/gene3')
      # We extract mapping to add a 'Top_GO_Terms' column to the Nodes table.
      
      # Extract mapping from GO back to genes
      go_mapping <- go_res %>% 
        filter(TF %in% top_tfs_for_network) %>%
        select(TF, Description, geneID) %>%
        separate_rows(geneID, sep = "/") %>%
        group_by(geneID) %>%
        summarise(Top_GO_Terms = paste(unique(Description)[1:min(3, length(unique(Description)))], collapse = " | "))
        
      cyto_nodes <- cyto_nodes %>%
         left_join(go_mapping, by = c("node_id" = "geneID")) %>%
         mutate(Top_GO_Terms = ifelse(is.na(Top_GO_Terms), "Unknown", Top_GO_Terms))
         
      cat(" -> Injected GO Enrichment terms into Nodes table!\n")
  }
  
  write.csv(cyto_nodes, paste0("Cytoscape_Nodes_Top5_MRs_", tier, ".csv"), row.names = FALSE)
}

cat("\n==== Phase 3 Complete! ====\n")
