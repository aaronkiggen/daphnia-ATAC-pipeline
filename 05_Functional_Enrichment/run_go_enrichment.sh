#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=GO_ENR
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=12G
#SBATCH --time=4:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"

# Delete broken env
rm -rf /data/leuven/354/vsc35429/miniconda3/envs/go_env_final

# Build clean env for plotting and GO
conda create -y -n go_env_final -c conda-forge -c bioconda r-base=4.3.0 bioconductor-clusterprofiler r-dplyr r-tidyr r-ggplot2 pango
conda activate go_env_final

cd /user/leuven/354/vsc35429/DATA/PhD/ATAC/
Rscript run_go_enrichment_hubs.R
