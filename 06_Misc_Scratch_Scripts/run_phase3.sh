#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=GRN_P3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:15:00
#SBATCH -A lp_svbelleghem
#SBATCH --output=slurm_phase3.out

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /scratch/leuven/354/vsc35429/ATAC/TF_analysis
Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase3_Visuals.R
