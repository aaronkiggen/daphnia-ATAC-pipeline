#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=REG_GEN
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

FASTA="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fasta"
GFF="/user/leuven/354/vsc35429/DATA/PhD/ATAC/GCA_030254905.1_UOB_LRV0_1_genomic.gff"
CHROM_SIZES="/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fa.sizes"
DEG_LIST="DEGs_genes_only.txt"

if [ ! -f "$DEG_LIST" ]; then
    echo "Missing $DEG_LIST"; exit 1;
fi

echo "Step 1: Extracting basic gene footprints..."
# 1. Base gene regions
awk '$3 == "gene"' $GFF | grep -wF -f $DEG_LIST > DEGs.gff
awk -v OFS='\t' '{ match($9, /ID=([^;]+)/, arr); print $1, $4-1, $5, arr[1], ".", $7 }' DEGs.gff > DEGs.bed

awk '$3 == "gene"' $GFF | grep -wF -v -f $DEG_LIST > nonDEGs.gff
awk -v OFS='\t' '{ match($9, /ID=([^;]+)/, arr); print $1, $4-1, $5, arr[1], ".", $7 }' nonDEGs.gff > nonDEGs.bed

echo "Step 2: Defining regional windows via slop..."
# 2. Slop to create base tier windows (Foreground/DEGs)
bedtools slop -i DEGs.bed -g $CHROM_SIZES -l 2000 -r 2000 -s > DEG_promoters.bed
bedtools slop -i DEGs.bed -g $CHROM_SIZES -l 5000 -r 2000 -s > DEG_proximal_full.bed
bedtools slop -i DEGs.bed -g $CHROM_SIZES -b 20000 > DEG_distal_full.bed

# 3. Slop to create base tier windows (Background/non-DEGs)
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -l 2000 -r 2000 -s > nonDEG_promoters.bed
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -l 5000 -r 5000 -s > nonDEG_proximal_full.bed
bedtools slop -i nonDEGs.bed -g $CHROM_SIZES -b 20000 > nonDEG_distal_full.bed

echo "Step 3: Strictly policing Russian Doll tier overlap..."
# 4. Russian Doll Subtraction on the Regions Themselves
# Because we are no longer selecting punctuted ATAC peaks, we use 'subtract' to literally 
# chisel the BED coordinates into perfectly mutually exclusive zones, splitting contiguous 
# regions on the boundaries.

# For Foregrounds (DEGs) - strictly extract the exclusive outer regions
cp DEG_promoters.bed REGION_DEG_promoters.bed
bedtools subtract -a DEG_proximal_full.bed -b REGION_DEG_promoters.bed > REGION_DEG_proximal.bed
bedtools subtract -a DEG_distal_full.bed -b DEG_proximal_full.bed > REGION_DEG_distal.bed

# For Backgrounds (non-DEGs) - strictly extract the exclusive regions AND remove overlap with DEGs entirely to keep them clean "backgrounds"
# First, Russian doll the tiers internally
cp nonDEG_promoters.bed temp_nonDEG_promoters.bed
bedtools subtract -a nonDEG_proximal_full.bed -b temp_nonDEG_promoters.bed > temp_nonDEG_proximal.bed
bedtools subtract -a nonDEG_distal_full.bed -b nonDEG_proximal_full.bed > temp_nonDEG_distal.bed

# Second, ensure a non-DEG background region does NOT cross functionally into a DEG's territory
bedtools subtract -a temp_nonDEG_promoters.bed -b DEG_promoters.bed > REGION_nonDEG_promoters.bed
bedtools subtract -a temp_nonDEG_proximal.bed -b DEG_proximal_full.bed > REGION_nonDEG_proximal.bed
bedtools subtract -a temp_nonDEG_distal.bed -b DEG_distal_full.bed > REGION_nonDEG_distal.bed

echo "Step 4: Compiling full genomic FASTA sequences for motif scanning..."
# 5. Extract FASTA sequences for HOMER/FIMO
bedtools getfasta -fi $FASTA -bed REGION_DEG_promoters.bed -name > regions_DEG_promoters.fa
bedtools getfasta -fi $FASTA -bed REGION_DEG_proximal.bed -name > regions_DEG_proximal.fa
bedtools getfasta -fi $FASTA -bed REGION_DEG_distal.bed -name > regions_DEG_distal.fa

bedtools getfasta -fi $FASTA -bed REGION_nonDEG_promoters.bed -name > regions_NON_DEG_promoters.fa
bedtools getfasta -fi $FASTA -bed REGION_nonDEG_proximal.bed -name > regions_NON_DEG_proximal.fa
bedtools getfasta -fi $FASTA -bed REGION_nonDEG_distal.bed -name > regions_NON_DEG_distal.fa

# Clean FASTA sequence headers for downstream tools
sed -i 's/^>::/>/' regions_DEG_*.fa
sed -i 's/^>::/>/' regions_NON_DEG_*.fa

echo "Successfully generated entirely Peak-Agnostic region sets!" 
