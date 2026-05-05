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

# We assume go_env_final is already built. If not, uncomment below:
# conda create -y -n go_env_final -c conda-forge -c bioconda r-base=4.3.0 bioconductor-clusterprofiler r-dplyr r-tidyr r-ggplot2 pango
conda activate go_env_final

cd /user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/05_Functional_Enrichment
Rscript 01_go_enrichment_hubs.R
