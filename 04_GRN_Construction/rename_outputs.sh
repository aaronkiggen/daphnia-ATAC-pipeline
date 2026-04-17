#!/bin/bash
DIR="/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"
cd $DIR

# Helper to rename with a clean suffix
rnm() {
  local OLD=$1
  local NEW=$2
  
  if [ -f "${OLD}.csv" ]; then mv "${OLD}.csv" "${NEW}_distal.csv"; fi
  if [ -f "${OLD}_promoter.csv" ]; then mv "${OLD}_promoter.csv" "${NEW}_promoter.csv"; fi
  if [ -f "${OLD}_proximal.csv" ]; then mv "${OLD}_proximal.csv" "${NEW}_proximal.csv"; fi
  if [ -f "${OLD}.pdf" ]; then mv "${OLD}.pdf" "${NEW}_distal.pdf"; fi
  if [ -f "${OLD}_promoter.pdf" ]; then mv "${OLD}_promoter.pdf" "${NEW}_promoter.pdf"; fi
  if [ -f "${OLD}_proximal.pdf" ]; then mv "${OLD}_proximal.pdf" "${NEW}_proximal.pdf"; fi
}

# Phase 1
rnm "homer_enriched_TFs" "01_Enriched_Motifs_TFs"
rnm "GRN_edges" "01_Base_Motif_Network_Edges"
rnm "top_TF_candidates" "01_Candidate_TFs_Summary"

# Phase 2
rnm "GRN_Biological_Edges_All" "02_Filtered_Bio_Network_Edges_All"
rnm "GRN_Biological_Edges_Tier1_Only" "02_Filtered_Bio_Network_Edges_HighConf"
rnm "GRN_Master_Regulators" "02_Network_Centrality_Scores"
rnm "top_TF_Biological_Hubs" "02_Top_Active_Hubs"

# Phase 3
rnm "Fig1_Master_Regulators_Bubble" "03_Plot_Hub_Enrichment"
rnm "Cytoscape_Nodes_Top5_MRs" "03_Cyto_RawNetwork_Nodes_Top5"
rnm "Cytoscape_Edges_Top5_MRs" "03_Cyto_RawNetwork_Edges_Top5"

# Phase 4
rnm "Cytoscape_Collapsed_Nodes_All_DEG" "04_Cyto_Collapsed_Nodes_All_DEGTFs"
rnm "Cytoscape_Collapsed_Edges_All_DEG" "04_Cyto_Collapsed_Edges_All_DEGTFs"

# Phase 5
rnm "Validation_Top20_TFs" "05_Ortholog_Validation_Top20_TFs"

echo "Files renamed successfully!"
