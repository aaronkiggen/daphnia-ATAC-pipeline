#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=REMAP_TF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH -A lp_svbelleghem
#SBATCH --output=slurm_mapping.out

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /scratch/leuven/354/vsc35429/ATAC/TF_analysis
Rscript make_tf_ortholog_mapping.R
