#!/bin/bash
for file in 01_grn_phase1_distal.R 02_grn_phase1_promoter.R 03_grn_phase1_proximal.R; do
  sed -i 's|"homer_out_custom_promoter/knownResults.txt"|"homer_out_custom_bg_promoter/knownResults.txt"|g' /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/$file
  sed -i 's|"homer_out_custom_proximal/knownResults.txt"|"homer_out_custom_bg_proximal/knownResults.txt"|g' /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/$file
  sed -i 's|"homer_out_custom/knownResults.txt"|"homer_out_custom_bg_distal/knownResults.txt"|g' /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/$file
  sed -i 's|homer_all <- bind_rows(homer_default, homer_custom)|homer_all <- homer_custom|g' /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/$file
done
