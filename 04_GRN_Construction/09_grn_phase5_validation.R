#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)

tiers <- c("promoter", "proximal", "distal")
output_dir <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output"
base_dir <- "/scratch/leuven/354/vsc35429/ATAC/TF_analysis"

motif_names <- read.table("/scratch/leuven/354/vsc35429/ATAC/motif_databases/motif_to_TF_name.txt", sep="", header=FALSE, col.names=c("matrix_id", "tf_name"), stringsAsFactors=FALSE)
jaspar_meta <- read.table("/scratch/leuven/354/vsc35429/ATAC/motif_databases/JASPAR2026_metadata.tsv", sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

tf_map <- inner_join(motif_names, jaspar_meta %>% select(matrix_id, uniprot_ids), by="matrix_id") %>%
  filter(uniprot_ids != "" & uniprot_ids != "-") %>%
  separate_rows(uniprot_ids, sep=",|::") %>%
  mutate(uniprot_accession = trimws(uniprot_ids)) %>% filter(uniprot_accession != "")

blast_out <- file.path(base_dir, "tf_blastp_daphnia.txt")
blast_res <- read.table(blast_out, sep="", col.names=c("qseqid", "sseqid", "pident", "evalue", "bitscore"), stringsAsFactors=FALSE)
blast_res$uniprot_accession <- sapply(strsplit(blast_res$qseqid, "\\|"), function(x) x[2])

# Link SSEQID (KAK IDs) to Daphnia OUZ56 IDs using the headers
headers <- readLines(file.path(base_dir, "daphnia_headers.txt"))
header_df <- data.frame(full_header = headers, stringsAsFactors=FALSE)
header_df$sseqid <- str_extract(header_df$full_header, "^>([^ ]+)", group=1)
header_df$daphnia_id <- str_extract(header_df$full_header, "(OUZ56_[0-9]+)")
header_df <- header_df %>% filter(!is.na(daphnia_id))

blast_res <- left_join(blast_res, header_df, by="sseqid") %>% filter(!is.na(daphnia_id))

for (tier in tiers) {
  cat("\n--- Validating Top 20 TFs for", tier, "tier ---\n")
  edge_file <- file.path(output_dir, paste0("GRN_Biological_Edges_All_", tier, ".csv"))
  if (!file.exists(edge_file)) next
  
  edges <- read.csv(edge_file, stringsAsFactors = FALSE)
  edges <- edges %>%
    mutate(collapsed_source = ifelse(!is.na(daphnia_gene_id) & daphnia_gene_id != "NA", daphnia_gene_id, tf_name))
  
  # Identify top 20 TFs (DEG TFs first, then by target count)
  tf_stats <- edges %>%
    group_by(collapsed_source) %>%
    summarise(
      target_count = n_distinct(gene_id),
      tf_is_deg = first(tf_is_deg),
      avg_motif_score = median(motif_score, na.rm=TRUE),
      max_motif_score = ifelse(all(is.na(motif_score)), NA, max(motif_score, na.rm=TRUE)),
      .groups = "drop"
    ) %>%
    arrange(desc(tf_is_deg), desc(target_count)) %>%
    slice_head(n = 20)
    
  # Link to BLAST metrics
  # For each TF, it came from some tf_name(s). What are the average BLAST metrics?
  # We extract the edges for these top 20 sources to see which Jaspar tfs map to them
  hub_edges <- edges %>% filter(collapsed_source %in% tf_stats$collapsed_source)
  
  tf_blast_info <- hub_edges %>% 
    select(collapsed_source, tf_name) %>% 
    distinct() %>%
    left_join(tf_map, by="tf_name") %>%
    left_join(blast_res, by=c("uniprot_accession", "collapsed_source" = "daphnia_id")) %>%
    filter(!is.na(bitscore)) %>%
    group_by(collapsed_source) %>%
    summarise(
      ortholog_blast_bitscore = ifelse(all(is.na(bitscore)), NA, max(bitscore, na.rm=TRUE)),
      ortholog_blast_pident = ifelse(all(is.na(pident)), NA, max(pident, na.rm=TRUE)),
      ortholog_blast_evalue = ifelse(all(is.na(evalue)), NA, min(evalue, na.rm=TRUE)),
      original_tf_names = paste(unique(tf_name), collapse="|"),
      .groups="drop"
    )
    
  # DEG target info
  deg_targets <- hub_edges %>%
    filter(!is.na(target_padj) & target_padj < 0.05) %>%
    group_by(collapsed_source) %>%
    summarise(
      deg_target_count = n_distinct(gene_id),
      median_motif_score_to_degs = median(motif_score, na.rm=TRUE),
      highest_motif_score_to_degs = ifelse(all(is.na(motif_score)), NA, max(motif_score, na.rm=TRUE)),
      .groups="drop"
    )
    
  final_validation <- tf_stats %>%
    left_join(tf_blast_info, by="collapsed_source") %>%
    left_join(deg_targets, by="collapsed_source") %>%
    arrange(desc(tf_is_deg), desc(target_count))
    
  write.csv(final_validation, file.path(output_dir, paste0("Validation_Top20_TFs_", tier, ".csv")), row.names=FALSE, quote=FALSE)
  cat("Wrote Validation_Top20_TFs_", tier, ".csv\n", sep="")
  print(head(final_validation, 5))
}
