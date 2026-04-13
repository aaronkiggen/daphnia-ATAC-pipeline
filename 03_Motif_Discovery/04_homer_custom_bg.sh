#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=MOTIF_BG
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
if [ -f "peaks_near_DEG_promoters.fa" ] && [ -f "peaks_near_NON_DEG_promoters.fa" ]; then
    mkdir -p homer_out_custom_bg_promoter
    echo "Running HOMER vs Custom Background (Promoter)"
    findMotifs.pl peaks_near_DEG_promoters.fa fasta homer_out_custom_bg_promoter/ \
      -fastaBg peaks_near_NON_DEG_promoters.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_custom_bg_promoter/homer.log
fi

# Run Proximal
if [ -f "peaks_near_DEG_proximal.fa" ] && [ -f "peaks_near_NON_DEG_proximal.fa" ]; then
    mkdir -p homer_out_custom_bg_proximal
    echo "Running HOMER vs Custom Background (Proximal)"
    findMotifs.pl peaks_near_DEG_proximal.fa fasta homer_out_custom_bg_proximal/ \
      -fastaBg peaks_near_NON_DEG_proximal.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_custom_bg_proximal/homer.log
fi

# Run Distal
if [ -f "peaks_near_DEG_distal.fa" ] && [ -f "peaks_near_NON_DEG_distal.fa" ]; then
    mkdir -p homer_out_custom_bg_distal
    echo "Running HOMER vs Custom Background (Distal)"
    findMotifs.pl peaks_near_DEG_distal.fa fasta homer_out_custom_bg_distal/ \
      -fastaBg peaks_near_NON_DEG_distal.fa \
      -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 \
      2> homer_out_custom_bg_distal/homer.log
fi

echo "Done running HOMER with Strict Target/Background mapping!"
