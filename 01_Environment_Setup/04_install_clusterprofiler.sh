#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=INSTALL
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

conda install -y -c bioconda bioconductor-clusterprofiler
# if it fails, fallback to R BiocManager
Rscript -e 'if (!require("clusterProfiler", quietly = TRUE)) { if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org"); BiocManager::install("clusterProfiler", ask=FALSE, update=FALSE) }'
