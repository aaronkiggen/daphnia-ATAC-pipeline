#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=MOTIF_PROM
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

INPUT_FA="peaks_near_DEG_promoters.fa"
if [ ! -f "$INPUT_FA" ]; then echo "Missing $INPUT_FA"; exit 1; fi

MOTIF_DB_MEME="${BASE_DIR}/motif_databases/combined_nonredundant.meme"
MOTIF_DB_HOMER="${BASE_DIR}/motif_databases/fixed_homer.motif"

mkdir -p homer_out_default_promoter homer_out_custom_promoter memechip_out_promoter fimo_promoter

findMotifs.pl $INPUT_FA fasta homer_out_default_promoter/ -p 8 -len 8,10,12 2> homer_out_default_promoter/homer.log
findMotifs.pl $INPUT_FA fasta homer_out_custom_promoter/ -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 2> homer_out_custom_promoter/homer.log

meme-chip -oc memechip_out_promoter -db $MOTIF_DB_MEME $INPUT_FA

fimo --oc fimo_promoter --parse-genomic-coord $MOTIF_DB_MEME $INPUT_FA
echo "Done Promoter"
