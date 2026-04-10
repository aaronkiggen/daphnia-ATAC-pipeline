library(dplyr)
edges <- read.csv("GRN_Biological_Edges_Tier1_Only.csv", stringsAsFactors=F)
print(head(edges %>% filter(tf_name == "FoxI-b")))
