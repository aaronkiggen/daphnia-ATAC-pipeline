#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=ATAC_TF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate atac_seq_fixed

BASE_DIR="/scratch/leuven/354/vsc35429/ATAC"
WORK_DIR="${BASE_DIR}/TF_analysis"
mkdir -p $WORK_DIR
cd $WORK_DIR

CONSENSUS="${BASE_DIR}/consensus_analysis/zorimar_consensus.bed"
FASTA="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fasta"
GFF="/user/leuven/354/vsc35429/DATA/PhD/ATAC/GCA_030254905.1_UOB_LRV0_1_genomic.gff"
CHROM_SIZES="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fa.sizes"
DEG_LIST="DEGs_genes_only.txt"

if [ ! -f "$DEG_LIST" ]; then
    echo "gene1" > $DEG_LIST
fi

echo "Step 1: Building gene windows..."
awk '$3 == "gene"' $GFF | grep -wF -f $DEG_LIST > DEGs.gff

awk -v OFS='\t' '{
    match($9, /ID=([^;]+)/, arr); 
    gene_name = arr[1]; 
    print $1, $4-1, $5, gene_name, ".", $7
}' DEGs.gff > DEGs.bed

# Generate windows
bedtools slop -i DEGs.bed -g $CHROM_SIZES -l 2000 -r 2000 -s > DEG_promoters.bed
bedtools slop -i DEGs.bed -g $CHROM_SIZES -l 5000 -r 2000 -s > DEG_proximal.bed
bedtools slop -i DEGs.bed -g $CHROM_SIZES -b 50000 > DEG_distal.bed

echo "Step 2: Intersecting..."
bedtools intersect -a $CONSENSUS -b DEG_promoters.bed -u > peaks_near_DEG_promoters.bed

bedtools intersect -a $CONSENSUS -b DEG_proximal.bed -u > all_proximal_peaks.bed
# Filter out peaks that are already assigned to promoters
bedtools intersect -v -a all_proximal_peaks.bed -b peaks_near_DEG_promoters.bed > peaks_near_DEG_proximal.bed

bedtools intersect -a $CONSENSUS -b DEG_distal.bed -u > all_distal_peaks.bed
# Filter out peaks that are already assigned to proximal (which includes promoters)
bedtools intersect -v -a all_distal_peaks.bed -b all_proximal_peaks.bed > peaks_near_DEG_distal.bed

# Get coordinate overlaps (from the final filtered sets)
bedtools intersect -a peaks_near_DEG_promoters.bed -b DEG_promoters.bed -wo > peak_DEG_overlaps_promoter.txt
bedtools intersect -a peaks_near_DEG_proximal.bed -b DEG_proximal.bed -wo > peak_DEG_overlaps_proximal.txt
bedtools intersect -a peaks_near_DEG_distal.bed -b DEG_distal.bed -wo > peak_DEG_overlaps_distal.txt

echo "Step 3: Extracting sequences..."
bedtools getfasta -fi $FASTA -bed peaks_near_DEG_promoters.bed -name > peaks_near_DEG_promoters.fa
bedtools getfasta -fi $FASTA -bed peaks_near_DEG_proximal.bed -name > peaks_near_DEG_proximal.fa
bedtools getfasta -fi $FASTA -bed peaks_near_DEG_distal.bed -name > peaks_near_DEG_distal.fa
bedtools getfasta -fi $FASTA -bed $CONSENSUS -name > all_consensus_peaks.fa

# Fix format for meme 
sed -i 's/^>::/>/' peaks_near_DEG_promoters.fa peaks_near_DEG_proximal.fa peaks_near_DEG_distal.fa all_consensus_peaks.fa

echo "Completed!"
