# Phase 4: Collapse redundant mammalian TF motifs into unique Daphnia Master Regulators
# This script deduplicates the GRN by mapping multiple motif hits to their single Daphnia ortholog.

library(dplyr)

tiers <- c("promoter", "proximal", "distal")
output_dir <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"

for (tier in tiers) {
  cat("\n--- Collapsing", tier, "tier ---\n")
  
  edge_file <- file.path(output_dir, paste0("02_Filtered_Bio_Network_Edges_All_", tier, ".csv"))
  if (!file.exists(edge_file)) {
    cat("Missing", edge_file, "- skipping.\n")
    next
  }
  
  edges <- read.csv(edge_file, stringsAsFactors = FALSE)
  
  # 1. Create a true collapsed TF source name
  edges <- edges %>%
    mutate(collapsed_source = ifelse(!is.na(daphnia_gene_id) & daphnia_gene_id != "NA", 
                                     daphnia_gene_id, 
                                     tf_name))
  
  # 2. Collapse the redundant edges
  collapsed_edges <- edges %>%
    group_by(collapsed_source, target = gene_id) %>%
    arrange(desc(direction != "unmapped"), desc(biological_edge_score)) %>%
    summarise(
      weight = max(biological_edge_score, na.rm = TRUE),
      interaction = first(direction), # Now grabs activation/repression over unmapped
      original_motifs = paste(unique(tf_name), collapse = "|"),
      .groups = "drop"
    ) %>%
    mutate(norm_weight = weight / max(weight))
  
  # 3. Identify the true Top Collapsed Hubs (MUST BE DEGs!)
  deg_tfs <- edges %>% filter(tf_is_deg == TRUE) %>% pull(collapsed_source) %>% unique()
  
  hub_counts <- collapsed_edges %>%
    filter(collapsed_source %in% deg_tfs) %>%
    group_by(collapsed_source) %>%
    summarise(target_count = n()) %>%
    arrange(desc(target_count))
  
  top_hubs <- hub_counts$collapsed_source
  cat("All Collapsed DEG Hubs:\n")
  print(hub_counts)
  
  # 4. Filter Edges for only Top 5 Hubs
  top_edges <- collapsed_edges %>% filter(collapsed_source %in% top_hubs)
  
  # NEW: Tag edges between TFs
  top_edges <- top_edges %>%
    mutate(edge_type = ifelse(target %in% deg_tfs, "TF-TF", "TF-Target"))
  
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
  
  # Export final files WITHOUT quotes to fix Cytoscape strict matching
  write.csv(top_edges, file.path(output_dir, paste0("04_Cyto_Collapsed_Edges_All_DEGTFs_", tier, ".csv")), row.names = FALSE, quote = FALSE)
  write.csv(nodes, file.path(output_dir, paste0("04_Cyto_Collapsed_Nodes_All_DEGTFs_", tier, ".csv")), row.names = FALSE, quote = FALSE)
  cat("Exported Collapsed Cytoscape Networks!\n")
}
cat("\n==== Phase 4 Complete! ====\n")
