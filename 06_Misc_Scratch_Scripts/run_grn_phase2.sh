#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=ATAC_GRN_P2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /scratch/leuven/354/vsc35429/ATAC/TF_analysis
Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase2.R
