library(dplyr)
# Load homer enriched IDs
HOMER_DEFAULT <- "homer_out_default/knownResults.txt"
HOMER_CUSTOM  <- "homer_out_custom/knownResults.txt"

read_homer <- function(path) {
  df <- tryCatch(read.table(path, sep="\t", header=TRUE, quote="", comment.char=""), error=function(e) NULL)
  if(is.null(df)) return(NULL)
  colnames(df)[1:9] <- c("motif_name","consensus","pval","log_pval","qval","n_target","pct_target","n_bg","pct_bg")
  df$tf_name <- sub("\\(.*", "", df$motif_name)
  df$pct_target_num <- as.numeric(sub("%","", df$pct_target))
  df
}
homer_all <- bind_rows(read_homer(HOMER_DEFAULT), read_homer(HOMER_CUSTOM))
homer_enriched <- homer_all %>% filter(qval < 0.05, pct_target_num > 5)

# Pull IDs
# we need mapping from tf_name to Motif IDs for the default HOMER
# For CUSTOM, tf_name IS the Motif ID (MA...)
custom_ids <- homer_enriched %>% filter(grepl("^MA", tf_name)) %>% pull(tf_name)

# Default IDs: have to lookup JASPAR
tf_dict <- read.table("../motif_databases/motif_to_TF_name.txt", col.names=c("motif_id","tf_name"))
default_names <- homer_enriched %>% filter(!grepl("^MA", tf_name)) %>% pull(tf_name)
default_ids <- tf_dict %>% filter(toupper(tf_name) %in% toupper(default_names)) %>% pull(motif_id)

all_ids <- unique(c(custom_ids, default_ids))

writeLines(all_ids, "homer_enriched_motif_ids.txt")
