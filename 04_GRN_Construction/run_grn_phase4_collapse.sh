#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=GRN_Phase4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate go_env_final

cd /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction
Rscript 08_grn_phase4_collapse.R
