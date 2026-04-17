#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=PHASE3
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction
Rscript 07_grn_phase3_visuals.R
