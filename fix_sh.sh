#!/bin/bash
DIR="/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction"
sed -i 's|/user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_promoter.R|'$DIR'/02_grn_phase1_promoter.R|g' $DIR/run_grn_phase1_promoter.sh
sed -i 's|/user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase2_promoter.R|'$DIR'/05_grn_phase2_promoter.R|g' $DIR/run_grn_phase1_promoter.sh

sed -i 's|/user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_proximal.R|'$DIR'/03_grn_phase1_proximal.R|g' $DIR/run_grn_phase1_proximal.sh
sed -i 's|/user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase2_proximal.R|'$DIR'/06_grn_phase2_proximal.R|g' $DIR/run_grn_phase1_proximal.sh

sed -i 's|/user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase3_Visuals.R|'$DIR'/07_grn_phase3_visuals.R|g' $DIR/run_grn_phase3.sh
