library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ggplot2)

options(stringsAsFactors = FALSE)

# 1. Paths
BASE_DIR  <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"
GO_FILE   <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/combined_go_annotations.csv"
OUT_DIR   <- file.path(BASE_DIR, "GO_enrichment_results")

dir.create(OUT_DIR, showWarnings = FALSE)

# 2. Load Top TFs and Edges
hubs <- read.csv(file.path(BASE_DIR, "top_TF_Biological_Hubs_promoter.csv"))
edges <- read.csv(file.path(BASE_DIR, "GRN_Biological_Edges_Tier1_Only_promoter.csv"))

# Let's select Top 5 Master Regulators (TFs that are DEGs)
top_mrs <- hubs %>% filter(tf_is_deg == TRUE) %>% head(5) %>% pull(tf_name)

# If there are less than 5 MRs, fallback to Top 5 overall Hubs
if(length(top_mrs) < 5) {
   cat("Not enough DE TFs found. Using the top 5 hubs globally...\n")
   top_mrs <- hubs %>% head(5) %>% pull(tf_name)
}

cat("Top TFs selected for enrichment:", paste(top_mrs, collapse=", "), "\n")

# 3. Load and format GO annotations for each tool
cat("Loading GO annotations...\n")
go_data <- read.csv(GO_FILE, stringsAsFactors = FALSE)

# We have 3 tools: eggnog_GO, interpro_GO, fantasia_GO
# Let's melt them into term2gene format
build_term2gene <- function(df, tool_col) {
  df %>%
    select(gene, go_terms = !!sym(tool_col)) %>%
    filter(!is.na(go_terms) & go_terms != "") %>%
    separate_rows(go_terms, sep = ";") %>%
    select(term = go_terms, gene = gene) %>%
    distinct()
}

tools <- c("eggnog_GO", "interpro_GO", "fantasia_GO")

# 4. Perform Enrichment for each TF x Tool
all_results <- list()

for(tf in top_mrs) {
  cat("\nProcessing TF:", tf, "\n")
  
  # Get target genes for this TF
  target_genes <- edges %>% 
    filter(tf_name == tf) %>% 
    pull(gene_id) %>% 
    unique()
    
  if(length(target_genes) < 5) {
     cat("Skipping", tf, "- less than 5 targets.\n")
     next
  }
  
  for(tool in tools) {
     cat("  -> Running tool:", tool, "\n")
     
     term2gene <- build_term2gene(go_data, tool)
     # Universe is all genes annotated by this tool
     universe <- unique(term2gene$gene)
     
     if (nrow(term2gene) == 0) next
     
     # clusterProfiler enricher
     res <- enricher(
       gene = target_genes,
       pvalueCutoff = 0.05,
       qvalueCutoff = 0.05,
       universe = universe,
       TERM2GENE = term2gene
     )
     
     if(!is.null(res) && nrow(res) > 0) {
        res_df <- as.data.frame(res)
        res_df$TF <- tf
        res_df$Tool <- tool
        all_results[[paste(tf, tool, sep="_")]] <- res_df
        
        # Plot
        p <- dotplot(res, showCategory=15) + 
             ggtitle(paste("GO Enrichment:", tf, "(Targets) -", gsub("_GO", "", tool)))
             
        plot_file <- file.path(OUT_DIR, paste0(tf, "_", tool, "_dotplot.pdf"))
        ggsave(plot_file, p, width=8, height=6)
        cat("     Saved plot to", plot_file, "\n")
     } else {
        cat("     No significant enrichment found.\n")
     }
  }
}

# 5. Save all results to a single CSV
if(length(all_results) > 0) {
  final_df <- bind_rows(all_results)
  write.csv(final_df, file.path(OUT_DIR, "Master_Regulator_GO_Enrichment_All_Tools.csv"), row.names=FALSE)
  cat("\nAll significant results saved to Master_Regulator_GO_Enrichment_All_Tools.csv\n")
} else {
  cat("\nNo significant enrichment found for any tool/TF combination.\n")
}
