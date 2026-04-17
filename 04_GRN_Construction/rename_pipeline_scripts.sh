#!/bin/bash

cd /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction

# Phase 1 & 2 Literal Replacements (Distal)
sed -i 's|"homer_enriched_TFs.csv"|"01_Enriched_Motifs_TFs_distal.csv"|g' 01_grn_phase1_distal.R
sed -i 's|"GRN_edges.csv"|"01_Base_Motif_Network_Edges_distal.csv"|g' 01_grn_phase1_distal.R
sed -i 's|"top_TF_candidates.csv"|"01_Candidate_TFs_Summary_distal.csv"|g' 01_grn_phase1_distal.R

sed -i 's|"GRN_edges.csv"|"01_Base_Motif_Network_Edges_distal.csv"|g' 04_grn_phase2_distal.R
sed -i 's|"GRN_Biological_Edges_All.csv"|"02_Filtered_Bio_Network_Edges_All_distal.csv"|g' 04_grn_phase2_distal.R
sed -i 's|"GRN_Biological_Edges_Tier1_Only.csv"|"02_Filtered_Bio_Network_Edges_HighConf_distal.csv"|g' 04_grn_phase2_distal.R
sed -i 's|"GRN_Master_Regulators.csv"|"02_Network_Centrality_Scores_distal.csv"|g' 04_grn_phase2_distal.R
sed -i 's|"top_TF_Biological_Hubs.csv"|"02_Top_Active_Hubs_distal.csv"|g' 04_grn_phase2_distal.R

# Phase 1 & 2 Literal Replacements (Promoter)
sed -i 's|"homer_enriched_TFs_promoter.csv"|"01_Enriched_Motifs_TFs_promoter.csv"|g' 02_grn_phase1_promoter.R
sed -i 's|"GRN_edges_promoter.csv"|"01_Base_Motif_Network_Edges_promoter.csv"|g' 02_grn_phase1_promoter.R
sed -i 's|"top_TF_candidates_promoter.csv"|"01_Candidate_TFs_Summary_promoter.csv"|g' 02_grn_phase1_promoter.R

sed -i 's|"GRN_edges_promoter.csv"|"01_Base_Motif_Network_Edges_promoter.csv"|g' 05_grn_phase2_promoter.R
sed -i 's|"GRN_Biological_Edges_All_promoter.csv"|"02_Filtered_Bio_Network_Edges_All_promoter.csv"|g' 05_grn_phase2_promoter.R
sed -i 's|"GRN_Biological_Edges_Tier1_Only_promoter.csv"|"02_Filtered_Bio_Network_Edges_HighConf_promoter.csv"|g' 05_grn_phase2_promoter.R
sed -i 's|"GRN_Master_Regulators_promoter.csv"|"02_Network_Centrality_Scores_promoter.csv"|g' 05_grn_phase2_promoter.R
sed -i 's|"top_TF_Biological_Hubs_promoter.csv"|"02_Top_Active_Hubs_promoter.csv"|g' 05_grn_phase2_promoter.R

# Phase 1 & 2 Literal Replacements (Proximal)
sed -i 's|"homer_enriched_TFs_proximal.csv"|"01_Enriched_Motifs_TFs_proximal.csv"|g' 03_grn_phase1_proximal.R
sed -i 's|"GRN_edges_proximal.csv"|"01_Base_Motif_Network_Edges_proximal.csv"|g' 03_grn_phase1_proximal.R
sed -i 's|"top_TF_candidates_proximal.csv"|"01_Candidate_TFs_Summary_proximal.csv"|g' 03_grn_phase1_proximal.R

sed -i 's|"GRN_edges_proximal.csv"|"01_Base_Motif_Network_Edges_proximal.csv"|g' 06_grn_phase2_proximal.R
sed -i 's|"GRN_Biological_Edges_All_proximal.csv"|"02_Filtered_Bio_Network_Edges_All_proximal.csv"|g' 06_grn_phase2_proximal.R
sed -i 's|"GRN_Biological_Edges_Tier1_Only_proximal.csv"|"02_Filtered_Bio_Network_Edges_HighConf_proximal.csv"|g' 06_grn_phase2_proximal.R
sed -i 's|"GRN_Master_Regulators_proximal.csv"|"02_Network_Centrality_Scores_proximal.csv"|g' 06_grn_phase2_proximal.R
sed -i 's|"top_TF_Biological_Hubs_proximal.csv"|"02_Top_Active_Hubs_proximal.csv"|g' 06_grn_phase2_proximal.R


# Phase 3 Hardcoded tier replacing
# Remove the old if-else distal logic to use unified naming
cat << 'PHASE3' > patch_phase3.R
  # 1. LOAD PHASE 2 OUTPUTS
  hubs_file <- paste0("02_Top_Active_Hubs_", tier, ".csv")
  edges_file <- paste0("02_Filtered_Bio_Network_Edges_HighConf_", tier, ".csv")
PHASE3
sed -i '/if(tier == "distal") {/,/  }/c\'"$(cat patch_phase3.R)" 07_grn_phase3_visuals.R
rm patch_phase3.R

sed -i 's|"Fig1_Master_Regulators_Bubble_"|"03_Plot_Hub_Enrichment_"|g' 07_grn_phase3_visuals.R
sed -i 's|"Cytoscape_Edges_Top5_MRs_"|"03_Cyto_RawNetwork_Edges_Top5_"|g' 07_grn_phase3_visuals.R
sed -i 's|"Cytoscape_Nodes_Top5_MRs_"|"03_Cyto_RawNetwork_Nodes_Top5_"|g' 07_grn_phase3_visuals.R

# Phase 4 Hardcoded string replacing
sed -i 's|"GRN_Biological_Edges_All_"|"02_Filtered_Bio_Network_Edges_All_"|g' 08_grn_phase4_collapse.R
sed -i 's|"Cytoscape_Collapsed_Edges_All_DEG_"|"04_Cyto_Collapsed_Edges_All_DEGTFs_"|g' 08_grn_phase4_collapse.R
sed -i 's|"Cytoscape_Collapsed_Nodes_All_DEG_"|"04_Cyto_Collapsed_Nodes_All_DEGTFs_"|g' 08_grn_phase4_collapse.R

# Phase 5 Hardcoded string replacing
sed -i 's|"GRN_Biological_Edges_All_"|"02_Filtered_Bio_Network_Edges_All_"|g' 09_grn_phase5_validation.R
sed -i 's|Validation_Top20_TFs_|05_Ortholog_Validation_Top20_TFs_|g' 09_grn_phase5_validation.R

echo "Done modifying scripts!"
