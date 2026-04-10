#!/bin/bash
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate atac_seq_fixed

echo "Creating the install log..."
mamba install -v -y -c bioconda -c conda-forge meme homer bioconductor-universalmotif r-dplyr r-readr r-tidyr r-stringr > mamba_install.log 2>&1
echo "Mamba install finished with code $?" >> mamba_install.log
