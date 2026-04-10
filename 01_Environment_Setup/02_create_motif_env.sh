#!/bin/bash
source ~/.bashrc
eval "$(conda shell.bash hook)"

echo "Creating highly specific environment for Motif discovery..."
mamba create -p /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env -y -c bioconda -c conda-forge \
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

echo "Environment created correctly in miniconda path!"
