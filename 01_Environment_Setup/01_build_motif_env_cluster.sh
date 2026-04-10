#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=MOTIF_ENV
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
conda activate base

# Force conda to store its bulky caching on scratch/data instead of tiny /user home
export CONDA_PKGS_DIRS=/data/leuven/354/vsc35429/miniconda3/pkgs
mkdir -p $CONDA_PKGS_DIRS

echo "Building isolated motif environment purely on /data mount..."
conda create -p /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env -y \
  -c bioconda -c conda-forge \
  meme \
  homer \
  bioconductor-universalmotif \
  r-dplyr \
  r-readr \
  r-tidyr \
  r-stringr \
  pandas \
  wget \
  unzip \
  bedtools

echo "Environment cleanly built at: /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env"
