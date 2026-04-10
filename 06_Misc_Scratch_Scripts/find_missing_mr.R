library(dplyr)
degs <- read.csv("DEGs.csv") %>% filter(!is.na(log2FoldChange), padj < 0.05)
map <- read.csv("tf_ortholog_mapping.csv", stringsAsFactors=F)

# degs in tf database
mr <- map %>% filter(daphnia_gene_id %in% degs$gene_id)
print(paste("Total MRs in ortholog map:", length(unique(mr$tf_name))))
print(unique(mr$tf_name))
