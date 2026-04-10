#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=ATAC_MOTIF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# Activate environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

BASE_DIR="/scratch/leuven/354/vsc35429/ATAC"
WORK_DIR="${BASE_DIR}/TF_analysis"
cd $WORK_DIR

echo "==== Motif Discovery Pipeline ===="

if [ ! -f "peaks_near_DEGs.fa" ] || [ ! -f "all_consensus_peaks.fa" ]; then
    echo "ERROR: Missing FASTA files. Run atac_TF.sh first."
    exit 1
fi

# Fix bedtools `::` naming artifact in FASTA headers
sed -i 's/^>::/>/' peaks_near_DEGs.fa all_consensus_peaks.fa

MOTIF_DB_MEME="${BASE_DIR}/motif_databases/combined_nonredundant.meme"
MOTIF_DB_HOMER="${BASE_DIR}/motif_databases/fixed_homer.motif"

echo "Using MEME Database: $MOTIF_DB_MEME"
echo "Using HOMER Database: $MOTIF_DB_HOMER"

# ------------------------------------------------------------------------------
# 1. HOMER (findMotifs.pl) - Hypothesis Generator
# ------------------------------------------------------------------------------
echo "Running HOMER..."

# A. Run with HOMER's default database (e.g., standard insect/vertebrate motifs)
# Using automatic genomic background instead of peak background
echo " -> HOMER with Default Native Database"
mkdir -p homer_out_default
findMotifs.pl peaks_near_DEGs.fa fasta homer_out_default/ -p 8 -len 8,10,12 2> homer_out_default/homer.log

# B. Run with our Custom JASPAR/CIS-BP Daphnia merged database
# Using automatic genomic background instead of peak background
echo " -> HOMER with Custom Daphnia/JASPAR Database"
mkdir -p homer_out_custom
findMotifs.pl peaks_near_DEGs.fa fasta homer_out_custom/ -p 8 -mknown $MOTIF_DB_HOMER -len 8,10,12 2> homer_out_custom/homer.log

# ------------------------------------------------------------------------------
# 2. MEME-ChIP - Validation & De Novo Discovery
# ------------------------------------------------------------------------------
echo "Running MEME-ChIP..."
# MEME-ChIP performs de novo discovery AND scans against the known db
meme-chip -oc memechip_out -db $MOTIF_DB_MEME peaks_near_DEGs.fa

# ------------------------------------------------------------------------------
# 3. FIMO - Scan Prioritized Peaks
# ------------------------------------------------------------------------------
echo "Running FIMO..."

# Instead of limiting FIMO to just the de novo repeats from MEME-ChIP,
# we scan the peaks straight against the COMPLETE Daphnia/known Motif database
# (the R script will filter these thousands of hits down to just the ones
# HOMER identified as statistically enriched).

SHORTLIST_MOTIFS=$MOTIF_DB_MEME

fimo --oc fimo_DEGproximal --parse-genomic-coord $SHORTLIST_MOTIFS peaks_near_DEGs.fa

echo "Motif discovery script complete."
