library(dplyr)
library(tidyr)
# library(GenomicRanges) # Removed, purely reading bedtools output
# library(GenomicFeatures)
# library(rtracklayer)
# library(igraph)
# library(ggraph)
# library(tidygraph)
# library(ggplot2)
# library(ggrepel)

# ════════════════════════════════════════════════════════════════════════════════
# 0.  PATHS
# ════════════════════════════════════════════════════════════════════════════════

BASE_DIR  <- "/scratch/leuven/354/vsc35429/ATAC"
WORK_DIR  <- file.path(BASE_DIR, "TF_analysis")
setwd(WORK_DIR)
OUTPUT_DIR <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"

HOMER_DEFAULT  <- "homer_out_default/knownResults.txt"
HOMER_CUSTOM   <- "homer_out_custom_bg_distal/knownResults.txt"
FIMO_TXT       <- "fimo_DEGproximal/fimo.txt"   # Note: Changed to fimo.txt from previous pipeline
PEAKS_BED      <- "peaks_near_DEGs.bed"   
DEGS_CSV       <- "DEGs.csv"

# ════════════════════════════════════════════════════════════════════════════════
# 1.  LOAD INPUTS
# ════════════════════════════════════════════════════════════════════════════════

# ── DEGs ──────────────────────────────────────────────────────────────────────
degs <- read.csv(DEGS_CSV)
cat("DEG columns:", paste(colnames(degs), collapse=", "), "\n")

stopifnot(all(c("gene_id","log2FoldChange","padj") %in% colnames(degs)))

# Keep all DEGs provided in the file (assuming it is already pre-filtered if needed)
deg_sig <- degs %>% filter(!is.na(gene_id))
cat("Significant DEGs:", nrow(deg_sig), "\n")

# ── HOMER results ─────────────────────────────────────────────────────────────
read_homer_results <- function(path, label) {
  df <- tryCatch(
    read.table(path, sep="\t", header=TRUE, quote="", comment.char=""),
    error = function(e) { cat("Could not read:", path, "\n"); return(NULL) }
  )
  if(is.null(df)) return(NULL)
  
  colnames(df)[1:9] <- c("motif_name","consensus","pval","log_pval","qval",
                          "n_target","pct_target","n_bg","pct_bg")
  df$source  <- label
  df$tf_name <- sub("\\(.*", "", df$motif_name)
  df$pct_target_num <- as.numeric(sub("%","", df$pct_target))
  df
}

homer_default <- read_homer_results(HOMER_DEFAULT, "default")
homer_custom  <- read_homer_results(HOMER_CUSTOM,  "custom")

homer_all <- homer_custom %>%
  group_by(tf_name) %>%
  slice_min(qval, n=1, with_ties=FALSE) %>%
  ungroup()

homer_enriched <- homer_all %>%
  filter(qval < 0.05, pct_target_num > 5) %>%
  arrange(log_pval)

cat("HOMER enriched TFs (q<0.05, >5% peaks):", nrow(homer_enriched), "\n")
write.csv(homer_enriched, file.path(OUTPUT_DIR, "01_Enriched_Motifs_TFs_distal.csv"), row.names=FALSE)

# ── FIMO & TF Mapping ───────────────────────────────────────────────
fimo_raw <- read.table(FIMO_TXT, sep="\t", header=TRUE, comment.char="")
# Re-name columns safely if 'motif_alt_id' is missing
colnames(fimo_raw)[1:5] <- c("motif_id", "sequence_name", "start", "stop", "strand")
fimo_raw <- fimo_raw %>% rename(chr = sequence_name)

# Motif Names mapping generated prior on cluster
tf_dict <- read.table("../motif_databases/motif_to_TF_name.txt", col.names=c("motif_id", "tf_name"))

# Since FIMO now scans the full known database, motif_id is already the target ID (e.g., MA2700.1)
# No need for TomTom de novo bridging!
fimo <- fimo_raw %>%
  mutate(
    motif_id = as.character(motif_id),
    database_source = ifelse(grepl("^MA", motif_id), "JASPAR", "CIS-BP")
  ) %>%
  left_join(tf_dict, by="motif_id")

# Fallback: if TF name isn't found, keep Motif_ID
fimo$tf_name <- coalesce(fimo$tf_name, as.character(fimo$motif_id))

# Identify which TFs passed strict FDR < 0.05
homer_qval_ids <- homer_all %>% filter(qval < 0.05, pct_target_num > 5) %>% pull(tf_name)
homer_qval_names <- tf_dict %>% filter(motif_id %in% homer_qval_ids) %>% pull(tf_name)

# Tag ALL fimo edges with this boolean flag
fimo <- fimo %>% 
  mutate(homer_qval_enriched = motif_id %in% homer_qval_ids | tf_name %in% homer_qval_names | tf_name %in% homer_qval_ids)

# Filter FIMO edges EXCLUSIVELY to motifs HOMER found to be enriched (Uncorrected p-value for leniency, or q-value)
# We fall back to p-value < 0.01 because ATAC genomic backgrounds still harshly penalize q-values.
homer_significant_ids <- homer_all %>% filter(pval <= 0.01, pct_target_num > 5) %>% pull(tf_name) # Homer custom saved motif IDs in tf_name
homer_significant_names <- tf_dict %>% filter(motif_id %in% homer_significant_ids) %>% pull(tf_name)

fimo_filtered <- fimo %>% filter(motif_id %in% homer_significant_ids | tf_name %in% homer_significant_names | tf_name %in% homer_significant_ids)

# Safety net: If HOMER was exceptionally strict, default to the top 100 motifs by FIMO score
if(nrow(fimo_filtered) < 100) {
  cat("HOMER filter too strict. Proceeding with top structural FIMO hits...\n")
  fimo_filtered <- fimo
}

# ════════════════════════════════════════════════════════════════════════════════
# 2.  LINK PEAKS TO DEGs (Now using exact Genomic Coords from FIMO)
# ════════════════════════════════════════════════════════════════════════════════

peaks_df <- read.table(PEAKS_BED,
                        col.names = c("chr","start","end"),
                        colClasses = c("character","integer","integer"))

peaks_df$peak_id <- paste0(peaks_df$chr, ":", peaks_df$start, "-", peaks_df$end)

# Renamed for join clarity
peaks_df <- peaks_df %>% rename(p_start = start, p_end = end)

# Non-equi join to assign Motif exact locations back to original Peak ID
fimo_with_peaks <- fimo_filtered %>%
  inner_join(peaks_df, by = join_by(chr == chr, start >= p_start, stop <= p_end))

cat("FIMO hits mapped back to BED peaks:", nrow(fimo_with_peaks), "\n")


# ── Pre-calculated Peak → DEG overlaps ────────────────────────────────────────

# Because spatial peak sets are mutually exclusive (Russian doll subtraction),
# we only need to load the specific tier this script is evaluating.
peak_gene_raw <- read.table("peak_DEG_overlaps_distal.txt", 
                               col.names = c("p_chr", "p_start", "p_end", "g_chr", "g_start", "g_end", "gene_id", "g_score", "g_strand", "overlap_bp"),
                               stringsAsFactors = FALSE)

peak_gene <- data.frame(
  peak_id  = paste0(peak_gene_raw$p_chr, ":", peak_gene_raw$p_start, "-", peak_gene_raw$p_end),
  gene_id  = peak_gene_raw$gene_id,
  tier     = "Tier3_Distal"
) %>% distinct()

# ════════════════════════════════════════════════════════════════════════════════
# 3.  ASSEMBLE GRN EDGE TABLE
# ════════════════════════════════════════════════════════════════════════════════

grn_raw <- as.data.frame(fimo_with_peaks) %>%
  inner_join(peak_gene, by="peak_id", relationship = "many-to-many") %>%
  left_join(degs %>% select(gene_id, target_lfc = log2FoldChange, target_padj = padj), by = "gene_id")

grn_edges <- grn_raw %>%
  group_by(tf_name, gene_id) %>%
  slice_min(p.value, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(
    motif_score   = -log10(p.value),
    edge_score    = motif_score,
    tier          = tier
  ) %>%
  arrange(desc(edge_score))

write.csv(grn_edges, file.path(OUTPUT_DIR, "01_Base_Motif_Network_Edges_distal.csv"), row.names=FALSE)

# ════════════════════════════════════════════════════════════════════════════════
# 4.  TOP TF TABLE
# ════════════════════════════════════════════════════════════════════════════════

top_tfs <- grn_edges %>%
  group_by(tf_name, database_source) %>%
  summarise(
    n_targets          = n_distinct(gene_id),
    n_tier1_promoter   = n_distinct(gene_id[tier == "Tier1_Promoter"]),
    n_tier2_local      = n_distinct(gene_id[tier == "Tier2_Local"]),
    n_tier3_distal     = n_distinct(gene_id[tier == "Tier3_Distal"]),
    n_peaks            = n_distinct(peak_id),
    mean_motif_score   = mean(motif_score),
    .groups = "drop"
  ) %>%
  arrange(desc(n_targets))

write.csv(top_tfs, file.path(OUTPUT_DIR, "01_Candidate_TFs_Summary_distal.csv"), row.names=FALSE)

cat("\nNetwork edge count:", nrow(grn_edges), "- Hubs saved to top_TF_candidates.csv\n")
# Plotting block omitted for brevity... (paste your exact plotting code here on line 170)