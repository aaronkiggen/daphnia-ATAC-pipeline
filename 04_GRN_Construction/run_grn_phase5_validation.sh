#!/bin/bash
#SBATCH --job-name=grn_val
#SBATCH --output=/scratch/leuven/354/vsc35429/ATAC/logs/grn_validation_%j.out
#SBATCH --error=/scratch/leuven/354/vsc35429/ATAC/logs/grn_validation_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate go_env_final

Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/09_grn_phase5_validation.R
