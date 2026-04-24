#!/usr/bin/env Rscript
# ----------------------------------------------------------
# Phase 0: Map discovered TF motifs to Daphnia gene orthologs
# This step is a prerequisite for all GRN construction scripts using TF to gene mapping.
# Script adapted from original: 06_Misc_Scratch_Scripts/make_tf_ortholog_mapping.R
#
# This step uses JASPAR metadata to map motif IDs to Uniprot IDs, runs BLASTP, and merges results
# to produce tf_ortholog_mapping.csv needed in later pipeline steps (Phase 2+)
# ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

BASE_DIR <- "/scratch/leuven/354/vsc35429/ATAC/TF_analysis"
setwd(BASE_DIR)

# 1. READ MOTIFS 
motif_names <- read.table("../motif_databases/motif_to_TF_name.txt", sep="", header=FALSE, col.names=c("matrix_id", "tf_name"), stringsAsFactors=FALSE)
cat("Found", nrow(motif_names), "motifs in ATAC FIMO analysis.\n")

# 2. READ JASPAR METADATA to get UNIPROT IDs
cat("Reading JASPAR metadata...\n")
jaspar_meta <- read.table("../motif_databases/JASPAR2026_metadata.tsv", sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

# Combine mapping
tf_map <- inner_join(motif_names, jaspar_meta %>% select(matrix_id, uniprot_ids), by="matrix_id")

# Filter out empty uniprot IDs
tf_map <- tf_map %>% filter(uniprot_ids != "" & uniprot_ids != "-")
cat("Found", nrow(tf_map), "motifs with Uniprot IDs to query.\n")

# Split by comma and :: 
tf_map <- tf_map %>% 
  separate_rows(uniprot_ids, sep=",|::") %>%
  mutate(uniprot_ids = trimws(uniprot_ids)) %>%
  filter(uniprot_ids != "")

uniprots <- unique(tf_map$uniprot_ids)
cat("Expanded to", length(uniprots), "unique Uniprot proteins...\n")

fasta_out <- "tf_proteins_new.fasta"
cat("Using existing", fasta_out, "which has already been safely downloaded.\n")

# 3. RUN BLAST
protein_fa <- "/user/leuven/354/vsc35429/DATA/PhD/QTL_BLAST/ncbi_dataset/data/GCA_030254905.1/protein.faa"

blast_out <- "tf_blastp_daphnia.txt"
cat("Running BLASTP...\n")
# Use max_target_seqs 1 for top Daphnia ortholog
cmd <- paste0("blastp -query ", fasta_out, " -db ", protein_fa, " -outfmt '6 qseqid sseqid pident evalue bitscore' -evalue 1e-5 -max_target_seqs 1 -num_threads 4 > ", blast_out)
system(cmd)

cat("BLASTP finished.\n")

# 4. PARSE BLAST AND MERGE
blast_res <- read.table(blast_out, sep="", col.names=c("qseqid", "sseqid", "pident", "evalue", "bitscore"), stringsAsFactors=FALSE)

# qseqid in uniprot fasta is format 'sp|P29383|AGL3_ARATH'
blast_res$uniprot_accession <- sapply(strsplit(blast_res$qseqid, "\\|"), function(x) x[2])

cat("Parsing Daphnia FASTA headers to map Accession to OUZ56 IDs...\n")
system(paste0("grep '^>' ", protein_fa, " > daphnia_headers.txt"))
headers <- readLines("daphnia_headers.txt")
header_df <- data.frame(full_header = headers, stringsAsFactors=FALSE)
header_df$sseqid <- str_extract(header_df$full_header, "^>([^ ]+)", group=1)
header_df$daphnia_id <- str_extract(header_df$full_header, "(OUZ56_[0-9]+)")
header_df <- header_df %>% filter(!is.na(daphnia_id))

# Merge blast with header info
blast_res <- left_join(blast_res, header_df, by="sseqid")
blast_res <- blast_res %>% filter(!is.na(daphnia_id))

# Take the best bitscore per uniprot accession
blast_res <- blast_res %>% group_by(uniprot_accession) %>% slice_max(order_by=bitscore, n=1, with_ties=FALSE) %>% ungroup()

# Merge back to tf_map to get TF Name to Daphnia ID
tf_ortholog_mapping <- tf_map %>% 
  left_join(blast_res %>% select(uniprot_accession, daphnia_gene_id = daphnia_id), by=c("uniprot_ids" = "uniprot_accession")) %>% 
  filter(!is.na(daphnia_gene_id)) %>%
  select(tf_name, daphnia_gene_id) %>%
  distinct()

cat("Successfully mapped", nrow(tf_ortholog_mapping), "TF names to Daphnia IDs.\n")
write.csv(tf_ortholog_mapping, "tf_ortholog_mapping.csv", row.names=FALSE, quote=FALSE)
cat("Wrote tf_ortholog_mapping.csv\n")
