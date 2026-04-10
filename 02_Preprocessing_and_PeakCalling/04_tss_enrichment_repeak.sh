#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name LDM_repeak
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=02:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=0-5

# LDM clones only have 6 samples (indices 0 to 5)

source /data/leuven/354/vsc35429/miniconda3/etc/profile.d/conda.sh
module --force purge

# 1. Find the already shifted BAM files specifically for LDM clones
BAM_LIST=( $(find /scratch/leuven/354/vsc35429/ATAC/LDM_clones/concatenated -name "*.shifted.sorted.bam" | sort) )

# 2. Select the input BAM based on the array task ID
INPUT=${BAM_LIST[$SLURM_ARRAY_TASK_ID]}
WORK_DIR=$(dirname $INPUT)
SAMPLE=$(basename $INPUT .shifted.sorted.bam)

cd $WORK_DIR
echo "Processing $SAMPLE in $WORK_DIR for relaxed peak calling"

# ── 1. Relaxed MACS2 peak calling ────────────────────────────────────────────
mkdir -p peaks_relaxed
module load cluster/wice/batch
module load MACS2/2.2.9.1-foss-2023a

# Relaxed parameters for shallow sequenced clones: 
# -p 0.01 instead of default q-value 
# -g 2.0e8 for tighter effective genome estimation
macs2 callpeak \
    -t ${SAMPLE}.shifted.sorted.bam \
    -f BAMPE \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    --keep-dup all \
    -g 2.0e8 \
    -p 0.01 \
    -n $SAMPLE \
    --outdir peaks_relaxed/ \
    2> peaks_relaxed/${SAMPLE}_macs2_relaxed.log

# ── 2. TSS enrichment profile (deepTools) ────────────────────────────────────
module --force purge
conda activate atac_seq_fixed

TSS_BED=/user/leuven/354/vsc35429/DATA/PhD/ATAC/Dmagna_TSS.bed

if [ -f "$TSS_BED" ]; then
    mkdir -p tss_qc_relaxed
    
    # We reuse the BigWig (.bw) files generated from the previous run
    BW_FILE="${SAMPLE}.shifted.bw"
    
    if [ -f "$BW_FILE" ]; then
        computeMatrix reference-point \
            --referencePoint TSS \
            --scoreFileName $BW_FILE \
            --regionsFileName $TSS_BED \
            --beforeRegionStartLength 2000 \
            --afterRegionStartLength 2000 \
            --skipZeros \
            -o tss_qc_relaxed/${SAMPLE}_tss_matrix.gz \
            -p 8

        plotProfile \
            -m tss_qc_relaxed/${SAMPLE}_tss_matrix.gz \
            --plotTitle "${SAMPLE} (Relaxed) TSS enrichment" \
            --outFileName tss_qc_relaxed/${SAMPLE}_TSS_profile.pdf
    else
        echo "WARNING: BigWig file $BW_FILE not found in $WORK_DIR"
    fi
else
    echo "WARNING: TSS BED not found, skipping TSS enrichment."
fi
