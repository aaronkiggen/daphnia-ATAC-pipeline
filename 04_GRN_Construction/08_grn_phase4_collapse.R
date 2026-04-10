# Phase 4: Collapse redundant mammalian TF motifs into unique Daphnia Master Regulators
# This script deduplicates the GRN by mapping multiple motif hits to their single Daphnia ortholog.

library(dplyr)

tiers <- c("promoter", "proximal", "distal")
output_dir <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output"

for (tier in tiers) {
  cat("\n--- Collapsing", tier, "tier ---\n")
  
  edge_file <- file.path(output_dir, paste0("GRN_Biological_Edges_All_", tier, ".csv"))
  if (!file.exists(edge_file)) {
    cat("Missing", edge_file, "- skipping.\n")
    next
  }
  
  edges <- read.csv(edge_file, stringsAsFactors = FALSE)
  
  # 1. Create a true collapsed TF source name: Map strictly to Daphnia gene ID if available
  edges <- edges %>%
    mutate(collapsed_source = ifelse(!is.na(daphnia_gene_id) & daphnia_gene_id != "NA", 
                                     daphnia_gene_id, 
                                     tf_name))
  
  # 2. Collapse the redundant edges by taking the maximum biological edge score
  collapsed_edges <- edges %>%
    group_by(collapsed_source, target = gene_id) %>%
    summarise(
      weight = max(biological_edge_score, na.rm = TRUE),
      interaction = first(direction), # Keep activation/repression
      # Concatenate the original mammalian names just for reference metadata
      original_motifs = paste(unique(tf_name), collapse = "|"),
      .groups = "drop"
    ) %>%
    # Normalize weights for Cytoscape (scale 0-1)
    mutate(norm_weight = weight / max(weight))
  
  # 3. Identify the true Top Collapsed Hubs (Most Targets)
  hub_counts <- collapsed_edges %>%
    group_by(collapsed_source) %>%
    summarise(target_count = n()) %>%
    arrange(desc(target_count))
  
  top_hubs <- head(hub_counts$collapsed_source, 5)
  cat("Top 5 Collapsed Hubs:\n")
  print(hub_counts[1:5, ])
  
  # 4. Filter Edges for only Top 5 Hubs
  top_edges <- collapsed_edges %>% filter(collapsed_source %in% top_hubs)
  
  # 5. Create Cytoscape Node Table
  sources <- data.frame(node_id = unique(top_edges$collapsed_source), node_type = "Collapsed Transcription Factor")
  targets <- data.frame(node_id = unique(top_edges$target), node_type = "Target Gene")
  nodes <- bind_rows(sources, targets) %>% distinct(node_id, .keep_all = TRUE)
  
  # Inject target DEG log2FoldChange and padj
  target_metadata <- edges %>% 
    select(gene_id, target_lfc, target_padj) %>% 
    distinct()
  nodes <- nodes %>%
    left_join(target_metadata, by = c("node_id" = "gene_id")) %>%
    mutate(lfc = target_lfc, padj = target_padj) %>%
    select(node_id, node_type, lfc, padj)
  
  # Export final files
  write.csv(top_edges, file.path(output_dir, paste0("Cytoscape_Collapsed_Edges_Top5_", tier, ".csv")), row.names = FALSE)
  write.csv(nodes, file.path(output_dir, paste0("Cytoscape_Collapsed_Nodes_Top5_", tier, ".csv")), row.names = FALSE)
  cat("Exported Collapsed Cytoscape Networks!\n")
}
cat("\n==== Phase 4 Complete! ====\n")
