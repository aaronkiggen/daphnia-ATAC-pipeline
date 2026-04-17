import sys

with open("/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/07_grn_phase3_visuals.R", "r") as f:
    content = f.read()

old_block = """  # Check if the distal file exists without the suffix or with it
  if(tier == "distal") {
     hubs_file <- "top_TF_Biological_Hubs.csv"
     edges_file <- "GRN_Biological_Edges_Tier1_Only.csv"
  } else {
     hubs_file <- paste0("top_TF_Biological_Hubs_", tier, ".csv")
     edges_file <- paste0("GRN_Biological_Edges_Tier1_Only_", tier, ".csv")
  }"""

new_block = """  hubs_file <- paste0("02_Top_Active_Hubs_", tier, ".csv")
  edges_file <- paste0("02_Filtered_Bio_Network_Edges_HighConf_", tier, ".csv")"""

content = content.replace(old_block, new_block)

# Also rename the outputs in Phase 3
content = content.replace('paste0("Fig1_Master_Regulators_Bubble_", tier, ".pdf")', 'paste0("03_Plot_Hub_Enrichment_", tier, ".pdf")')
content = content.replace('paste0("Cytoscape_Edges_Top5_MRs_", tier, ".csv")', 'paste0("03_Cyto_RawNetwork_Edges_Top5_", tier, ".csv")')
content = content.replace('paste0("Cytoscape_Nodes_Top5_MRs_", tier, ".csv")', 'paste0("03_Cyto_RawNetwork_Nodes_Top5_", tier, ".csv")')

with open("/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/07_grn_phase3_visuals.R", "w") as f:
    f.write(content)
print("Phase 3 Distal Loader Updated!")

