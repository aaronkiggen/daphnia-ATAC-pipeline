#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name ATAC_LDM_downstream
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=0-10

# Load DeepTools safely from Conda
source /data/leuven/354/vsc35429/miniconda3/etc/profile.d/conda.sh
module --force purge  # clear any system modules before conda
conda activate atac_seq_fixed
# 1. Dynamically find all concatenated BAM files for LDM and ZORIMAR
# This creates an array of all absolute paths to the *filtered.sorted.nd.bam files.
BAM_LIST=( $(find /scratch/leuven/354/vsc35429/ATAC -name "*.filtered.sorted.nd.bam" | grep "concatenated" | sort) )

# 2. Select the input BAM based on the array task ID
INPUT=${BAM_LIST[$SLURM_ARRAY_TASK_ID]}
WORK_DIR=$(dirname $INPUT)
SAMPLE=$(basename $INPUT .filtered.sorted.nd.bam)

cd $WORK_DIR
echo "Processing $SAMPLE in $WORK_DIR"

# ── 1. Index the dedup BAM (if not done) ─────────────────────────────────────
if [ ! -f "${INPUT}.bai" ]; then
    samtools index -@ 8 $INPUT
fi

# ── 3. Tn5 shift (+4 bp on +strand, -5 bp on -strand) ───────────────────────
alignmentSieve \
    --ATACshift \
    --bam ${INPUT} \
    --outFile ${SAMPLE}.shifted.bam \
    -p 8

samtools sort -@ 8 ${SAMPLE}.shifted.bam -o ${SAMPLE}.shifted.sorted.bam
samtools index -@ 8 ${SAMPLE}.shifted.sorted.bam
rm ${SAMPLE}.shifted.bam   # save space

# ── 4. Re-generate BigWig from shifted reads ────────────────────────
FASTA=/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fasta
SIZES=${FASTA%.fasta}.fa.sizes

bamCoverage \
    --bam ${SAMPLE}.shifted.sorted.bam \
    --outFileName ${SAMPLE}.shifted.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --extendReads \
    -p 8

# ── 5. MACS2 peak calling per sample ─────────────────────────────────────────
mkdir -p peaks
conda deactivate  # MACS2 is in the system modules, not conda
module --force purge
module load cluster/wice/batch
module load MACS2/2.2.9.1-foss-2023a

# Daphnia magna genome size ~250 Mb; use -g 2.5e8
# --nomodel because we already shifted; -f BAMPE for paired-end
macs2 callpeak \
    -t ${SAMPLE}.shifted.sorted.bam \
    -f BAMPE \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    --keep-dup all \
    -g 2.5e8 \
    -n $SAMPLE \
    --outdir peaks/ \
    2> peaks/${SAMPLE}_macs2.log

# ── 6. TSS enrichment profile (deepTools) ────────────────────────────────────
module --force purge
conda activate atac_seq_fixed
# Requires a BED file of TSS positions for Daphnia magna
# Build from GFF/GTF: awk '$3=="gene"' annotation.gff | ... or use your annotation
TSS_BED=/user/leuven/354/vsc35429/DATA/PhD/ATAC/Dmagna_TSS.bed

if [ -f "$TSS_BED" ]; then
    mkdir -p tss_qc
    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName ${SAMPLE}.shifted.bw \
        --regionsFileName $TSS_BED \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --skipZeros \
        -o tss_qc/${SAMPLE}_tss_matrix.gz \
        -p 8

    plotProfile \
        -m tss_qc/${SAMPLE}_tss_matrix.gz \
        --plotTitle "${SAMPLE} TSS enrichment" \
        --outFileName tss_qc/${SAMPLE}_TSS_profile.pdf
else
    echo "WARNING: TSS BED not found, skipping TSS enrichment for ${SAMPLE}"
fi
