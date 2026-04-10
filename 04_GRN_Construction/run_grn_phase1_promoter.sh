#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=R_PROM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /scratch/leuven/354/vsc35429/ATAC/TF_analysis

Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_promoter.R
Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase2_promoter.R
