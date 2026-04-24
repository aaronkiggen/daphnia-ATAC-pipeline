#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=MOTIF_REG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

BASE_DIR="/scratch/leuven/354/vsc35429/ATAC"
WORK_DIR="${BASE_DIR}/TF_analysis"
cd $WORK_DIR

MOTIF_DB_HOMER="${BASE_DIR}/motif_databases/fixed_homer.motif"

# Run Promoters
if [ -f "regions_DEG_promoters.fa" ] && [ -f "regions_NON_DEG_promoters.fa" ]; then
    mkdir -p homer_out_regions_promoter
    echo "Running HOMER vs Region Agnostic Background (Promoter)"
    findMotifs.pl regions_DEG_promoters.fa fasta homer_out_regions_promoter/ \
      -fastaBg regions_NON_DEG_promoters.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_regions_promoter/homer.log
fi

# Run Proximal
if [ -f "regions_DEG_proximal.fa" ] && [ -f "regions_NON_DEG_proximal.fa" ]; then
    mkdir -p homer_out_regions_proximal
    echo "Running HOMER vs Region Agnostic Background (Proximal)"
    findMotifs.pl regions_DEG_proximal.fa fasta homer_out_regions_proximal/ \
      -fastaBg regions_NON_DEG_proximal.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_regions_proximal/homer.log
fi

# Run Distal
if [ -f "regions_DEG_distal.fa" ] && [ -f "regions_NON_DEG_distal.fa" ]; then
    mkdir -p homer_out_regions_distal
    echo "Running HOMER vs Region Agnostic Background (Distal)"
    findMotifs.pl regions_DEG_distal.fa fasta homer_out_regions_distal/ \
      -fastaBg regions_NON_DEG_distal.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_regions_distal/homer.log
fi

echo "Done running HOMER with Region Agnostic Target/Background mapping!"
