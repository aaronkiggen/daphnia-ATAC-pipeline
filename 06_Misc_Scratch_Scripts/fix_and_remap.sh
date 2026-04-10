#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=REMAP_FIX
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --output=slurm_remap_fix.out

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

cd /scratch/leuven/354/vsc35429/ATAC/TF_analysis

echo "1. Fetching missing 342 fastas (incl. FOXA1)..."
python3 fetch_missing.py

echo "2. Rebuilding tf_ortholog_mapping.csv with fixed split complexes..."
Rscript make_tf_ortholog_mapping.R

echo "3. Re-running Phase 2 Script to update GRN outputs..."
Rscript /user/leuven/354/vsc35429/DATA/PhD/ATAC/atac_TF_GRN_Phase2.R

echo "ALL DONE!"
