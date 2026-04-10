#!/bin/bash -l 
#SBATCH --cluster=wice 
#SBATCH --job-name ATAC_LDM 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --time=24:10:00 
#SBATCH -A lp_svbelleghem

# Activate conda
source /vsc-hard-mounts/leuven-data/354/vsc35429/miniconda3/etc/profile.d/conda.sh
conda activate atac_seq

# Set TMPDIR to scratch to avoid filling DATA filesystem with bowtie2 temp files
export TMPDIR=/scratch/leuven/354/vsc35429/ATAC/tmp
mkdir -p $TMPDIR

ID=$((SLURM_ARRAY_TASK_ID - 1))

samples=(GC177147_ATGTAGCG-CGCTACAT_S1 GC177147_CCAAGCAA-TTGCTTGG_S2 GC177147_CGTGATCA-TGATCACG_S3 GC177147_GACTTAGG-CCTAAGTC_S4 GC177147_GATAGGCT-AGCCTATC_S5 GC177147_TTGACAGG-CCTGTCAA_S6)

WORK_DIR=/scratch/leuven/354/vsc35429/ATAC/LDM_clones/concatenated
cd $WORK_DIR

# -------------------------
# Trimming
# -------------------------
trimmomatic PE \
${samples[$ID]}_R1_001.fastq.gz ${samples[$ID]}_R2_001.fastq.gz \
${samples[$ID]}_R1_trimmed.fastq.gz ${samples[$ID]}_R1_unpaired.fastq.gz \
${samples[$ID]}_R2_trimmed.fastq.gz ${samples[$ID]}_R2_unpaired.fastq.gz \
ILLUMINACLIP:/user/leuven/354/vsc35429/DATA/PhD/ATAC/TruSeq3-PE.fa:2:30:10:2:true \
SLIDINGWINDOW:4:20 \
LEADING:3 TRAILING:3 \
MINLEN:30

# -------------------------
# Reference preparation
# -------------------------
REF=/user/leuven/354/vsc35429/DATA/PhD/ATAC/LRV01.fna.gz
FASTA=/user/leuven/354/vsc35429/DATA/PhD/ATAC/DmagnaLRV01.fasta

# Unzip only if not already present
if [ ! -f "$FASTA" ]; then
    gunzip -c $REF > $FASTA
fi

# Index fasta (only if needed)
if [ ! -f "${FASTA}.fai" ]; then
    samtools faidx $FASTA
fi

# Chrom sizes
SIZES=${FASTA%.fasta}.fa.sizes
if [ ! -f "$SIZES" ]; then
    cut -f1,2 ${FASTA}.fai > $SIZES
fi

# Bowtie2 index must be pre-built and available
BOWTIE_INDEX=$WORK_DIR/DmagnaLRV01

if ! [ -f "${BOWTIE_INDEX}.1.bt2" ] || ! [ -f "${BOWTIE_INDEX}.rev.2.bt2" ]; then
    echo "ERROR: Bowtie2 index not found at $BOWTIE_INDEX"
    echo "Please build the index manually before running this script"
    exit 1
fi

# -------------------------
# Mapping
# -------------------------
bowtie2 -t -k 2 -p 8 --local -x $BOWTIE_INDEX \
-1 ${samples[$ID]}_R1_trimmed.fastq.gz \
-2 ${samples[$ID]}_R2_trimmed.fastq.gz | \
samtools view -bS - | \
samtools addreplacerg -r "@RG\tID:${samples[$ID]}\tSM:${samples[$ID]}" -o ${samples[$ID]}.bam -

# -------------------------
# Filtering + sorting
# -------------------------
samtools view -f 0x02 -q 20 -b ${samples[$ID]}.bam > ${samples[$ID]}.filtered.bam

samtools sort ${samples[$ID]}.filtered.bam -o ${samples[$ID]}.filtered.sorted.bam

# -------------------------
# Duplicate removal
# -------------------------
picard MarkDuplicates \
I=${samples[$ID]}.filtered.sorted.bam \
O=${samples[$ID]}.filtered.sorted.nd.bam \
REMOVE_DUPLICATES=true \
M=${samples[$ID]}_dup_metrics.txt \
ASSUME_SORTED=true

# -------------------------
# Coverage
# -------------------------
bedtools genomecov \
-ibam ${samples[$ID]}.filtered.sorted.nd.bam -bg \
> ${samples[$ID]}.filtered.sorted.nd.bdg

LC_COLLATE=C sort -k1,1 -k2,2n \
${samples[$ID]}.filtered.sorted.nd.bdg \
> ${samples[$ID]}.filtered.sorted.nd.collate.bdg

# -------------------------
# BigWig
# -------------------------
bedGraphToBigWig \
${samples[$ID]}.filtered.sorted.nd.collate.bdg \
$SIZES \
${samples[$ID]}.filtered.sorted.nd.collate.bw