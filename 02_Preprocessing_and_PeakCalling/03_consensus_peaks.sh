#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name ATAC_consensus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00
#SBATCH -A lp_svbelleghem

source /vsc-hard-mounts/leuven-data/354/vsc35429/miniconda3/etc/profile.d/conda.sh
conda activate atac_seq_fixed

# Directories
ATAC_DIR=/scratch/leuven/354/vsc35429/ATAC
OUT_DIR=$ATAC_DIR/consensus_analysis
mkdir -p $OUT_DIR
cd $OUT_DIR

echo "==================================================="
echo "1. ZORIMAR Consensus (>= 3 of 5)"
echo "==================================================="
ZORIMAR_PEAKS=(
    $ATAC_DIR/ZORIMAR_03112025/concatenated/peaks/GC172193_peaks.narrowPeak
    $ATAC_DIR/ZORIMAR_03112025/concatenated/peaks/GC172194_peaks.narrowPeak
    $ATAC_DIR/ZORIMAR_03112025/concatenated/peaks/GC172195_peaks.narrowPeak
    $ATAC_DIR/ZORIMAR_29092025/concatenated/peaks/GC170733_peaks.narrowPeak
    $ATAC_DIR/ZORIMAR_29092025/concatenated/peaks/GC170734_peaks.narrowPeak
)

cat "${ZORIMAR_PEAKS[@]}" | awk 'OFS="\t" {print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > zorimar_raw_union.bed
bedtools intersect -a zorimar_raw_union.bed -b "${ZORIMAR_PEAKS[@]}" -c | awk '$4 >= 3' | cut -f1-3 > zorimar_consensus.bed
echo "ZORIMAR consensus peaks: $(wc -l < zorimar_consensus.bed)"

echo "==================================================="
echo "2. LDM Consensus (>= 2 of 3 per condition)"
echo "==================================================="
# Based on Metadata: Control = S1, S2, S4 ; Fish = S3, S5, S6
LDM_DIR=$ATAC_DIR/LDM_clones/concatenated/peaks_relaxed

NOFISH_PEAKS=(
    $LDM_DIR/GC177147_ATGTAGCG-CGCTACAT_S1_peaks.narrowPeak
    $LDM_DIR/GC177147_CCAAGCAA-TTGCTTGG_S2_peaks.narrowPeak
    $LDM_DIR/GC177147_GACTTAGG-CCTAAGTC_S4_peaks.narrowPeak
)

FISH_PEAKS=(
    $LDM_DIR/GC177147_CGTGATCA-TGATCACG_S3_peaks.narrowPeak
    $LDM_DIR/GC177147_GATAGGCT-AGCCTATC_S5_peaks.narrowPeak
    $LDM_DIR/GC177147_TTGACAGG-CCTGTCAA_S6_peaks.narrowPeak
)

cat "${NOFISH_PEAKS[@]}" "${FISH_PEAKS[@]}" | awk 'OFS="\t" {print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > ldm_raw_union.bed

bedtools intersect -a ldm_raw_union.bed -b "${NOFISH_PEAKS[@]}" -c | awk '$4 >= 2' | cut -f1-3 > ldm_control_consensus.bed
bedtools intersect -a ldm_raw_union.bed -b "${FISH_PEAKS[@]}" -c | awk '$4 >= 2' | cut -f1-3 > ldm_fish_consensus.bed

echo "LDM Control (no-fish) consensus peaks: $(wc -l < ldm_control_consensus.bed)"
echo "LDM Fish consensus peaks: $(wc -l < ldm_fish_consensus.bed)"

cat ldm_control_consensus.bed ldm_fish_consensus.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ldm_overall_consensus.bed
echo "LDM Overall consensus peaks: $(wc -l < ldm_overall_consensus.bed)"

echo "==================================================="
echo "3. LDM Condition-Specific Peaks"
echo "==================================================="
# In Control but NOT in Fish
bedtools intersect -a ldm_control_consensus.bed -b ldm_fish_consensus.bed -v > ldm_control_specific.bed
echo "LDM Control-specific peaks: $(wc -l < ldm_control_specific.bed)"

# In Fish but NOT in Control
bedtools intersect -a ldm_fish_consensus.bed -b ldm_control_consensus.bed -v > ldm_fish_specific.bed
echo "LDM Fish-specific peaks: $(wc -l < ldm_fish_specific.bed)"

echo "==================================================="
echo "4. Overlap: ZORIMAR Consensus & LDM Control Consensus"
echo "==================================================="
bedtools intersect -a zorimar_consensus.bed -b ldm_control_consensus.bed -u > overlapping_zorimar_ldm_control.bed
echo "Peaks found in BOTH ZORIMAR and LDM Control: $(wc -l < overlapping_zorimar_ldm_control.bed)"

