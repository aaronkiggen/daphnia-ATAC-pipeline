# Phase 3 fixing nodes
library(dplyr)
hubs_file <- "top_TF_Biological_Hubs.csv"
edges_file <- "GRN_Biological_Edges_Tier1_Only.csv"
hubs  <- read.csv(hubs_file, stringsAsFactors = FALSE)
edges <- read.csv(edges_file, stringsAsFactors = FALSE)

mr_hubs <- hubs %>% filter(tf_is_deg == TRUE) %>% head(15) 
top_tfs_for_network <- mr_hubs$tf_name[1:min(5, nrow(mr_hubs))]

# Subset edges
cyto_edges <- edges %>%
  filter(tf_name %in% top_tfs_for_network) %>%
  select(
    source = tf_name, 
    target = gene_id,   # Using the target gene ID! It is `gene_id`
    interaction = direction, 
    weight = biological_edge_score
  ) %>%
  mutate(norm_weight = weight / max(weight, na.rm=TRUE))

write.csv(cyto_edges, "Cytoscape_Edges_Top5_MRs.csv", row.names = FALSE)
cat(" -> Saved Cytoscape_Edges_Top5_MRs.csv (", nrow(cyto_edges), " edges)\n")

# Build the Nodes table
sources <- edges %>%
  filter(tf_name %in% top_tfs_for_network) %>%
  select(node_id = tf_name, lfc = tf_lfc, padj = tf_padj) %>%
  mutate(node_type = "Transcription Factor", is_deg = TRUE) %>%
  distinct(node_id, .keep_all = TRUE)

targets <- edges %>%
  filter(tf_name %in% top_tfs_for_network) %>%
  select(node_id = gene_id, lfc = target_lfc, padj = target_padj) %>%
  mutate(node_type = "Target Gene", is_deg = (padj < 0.05)) %>%
  distinct(node_id, .keep_all = TRUE)

cyto_nodes <- bind_rows(sources, targets) %>%
  group_by(node_id) %>%
  slice(1) %>%
  ungroup()

write.csv(cyto_nodes, "Cytoscape_Nodes_Top5_MRs.csv", row.names = FALSE)
cat(" -> Saved Cytoscape_Nodes_Top5_MRs.csv (", nrow(cyto_nodes), " nodes)\n")
