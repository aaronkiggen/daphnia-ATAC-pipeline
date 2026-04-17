#!/bin/bash
DIR="/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction"
for file in 01_grn_phase1_distal.R 02_grn_phase1_promoter.R 03_grn_phase1_proximal.R 04_grn_phase2_distal.R 05_grn_phase2_promoter.R 06_grn_phase2_proximal.R; do
  # Add OUTPUT_DIR variable
  sed -i '/setwd(WORK_DIR)/a OUTPUT_DIR <- "/user/leuven/354/vsc35429/DATA/PhD/ATAC/output_custom_background"' $DIR/$file
  # Change write.csv(..., "file.csv") to write.csv(..., file.path(OUTPUT_DIR, "file.csv"))
  sed -i -E 's|write\.csv\(([^,]+),\s*"([^"]+\.csv)"|write.csv(\1, file.path(OUTPUT_DIR, "\2")|g' $DIR/$file
done
