#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=BG_GEN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate atac_seq_fixed

BASE_DIR="/scratch/leuven/354/vsc35429/ATAC"
WORK_DIR="${BASE_DIR}/TF_analysis"
cd $WORK_DIR

CONSENSUS="${BASE_DIR}/consensus_analysis/zorimar_consensus.bed"
FASTA="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fasta"
GFF="/user/leuven/354/vsc35429/DATA/PhD/ATAC/GCA_030254905.1_UOB_LRV0_1_genomic.gff"
CHROM_SIZES="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fa.sizes"
DEG_LIST="DEGs_genes_only.txt"

if [ ! -f "$DEG_LIST" ]; then
    echo "Missing $DEG_LIST"; exit 1;
fi

echo "Step 1: Extracting Non-DEG windows..."
# Get all genes, EXCLUDE the DEGs
awk '$3 == "gene"' $GFF | grep -wF -v -f $DEG_LIST > nonDEGs.gff

awk -v OFS='\t' '{
    match($9, /ID=([^;]+)/, arr); 
    gene_name = arr[1]; 
    print $1, $4-1, $5, gene_name, ".", $7
}' nonDEGs.gff > nonDEGs.bed

# Generate background windows
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -l 2000 -r 2000 -s > nonDEG_promoters.bed
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -l 5000 -r 2000 -s > nonDEG_proximal.bed
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -b 50000 > nonDEG_distal.bed

echo "Step 2: Intersecting and thoroughly filtering out DEG overlaps..."
# Find peaks near Non-DEGs, but use -v against DEG windows to ensure they are strictly background!
bedtools intersect -a $CONSENSUS -b nonDEG_promoters.bed -u | bedtools intersect -a - -b DEG_promoters.bed -v > peaks_near_NON_DEG_promoters.bed
bedtools intersect -a $CONSENSUS -b nonDEG_proximal.bed -u | bedtools intersect -a - -b DEG_proximal.bed -v > peaks_near_NON_DEG_proximal.bed
bedtools intersect -a $CONSENSUS -b nonDEG_distal.bed -u | bedtools intersect -a - -b DEG_distal.bed -v > peaks_near_NON_DEG_distal.bed

echo "Step 3: Extracting background sequences..."
bedtools getfasta -fi $FASTA -bed peaks_near_NON_DEG_promoters.bed -name > peaks_near_NON_DEG_promoters.fa
bedtools getfasta -fi $FASTA -bed peaks_near_NON_DEG_proximal.bed -name > peaks_near_NON_DEG_proximal.fa
bedtools getfasta -fi $FASTA -bed peaks_near_NON_DEG_distal.bed -name > peaks_near_NON_DEG_distal.fa

# Fix format for meme/homer
sed -i 's/^>::/>/' peaks_near_NON_DEG_promoters.fa peaks_near_NON_DEG_proximal.fa peaks_near_NON_DEG_distal.fa

echo "Completed!"
