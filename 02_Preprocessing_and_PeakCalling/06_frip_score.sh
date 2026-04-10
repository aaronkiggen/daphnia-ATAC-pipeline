#!/bin/bash
# Calculate FRiP using ZORIMAR consensus peaks for all samples

# SLURM options not strictly needed if running from interactive session, 
# but can be submitted. We'll run it sequentially.

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate atac_seq_fixed

BASE_DIR="/scratch/leuven/354/vsc35429/ATAC"
CONSENSUS="${BASE_DIR}/consensus_analysis/zorimar_consensus.bed"

OUTPUT_FILE="${BASE_DIR}/consensus_analysis/FRiP_summary_zorimar_consensus.txt"

echo "Calculating FRiP using ${CONSENSUS}"
echo -e "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP" > $OUTPUT_FILE

bams=$(find ${BASE_DIR} -name "*.shifted.sorted.bam" | grep concatenate)

for bam in $bams; do
    sample_name=$(basename $bam .shifted.sorted.bam)
    echo "Processing $sample_name..."
    
    # Mapped reads (using -F 260 to exclude unmapped and not primary alignments)
    total=$(samtools view -c -F 260 $bam)

    if [ "$total" -eq 0 ]; then
        frip="0.000"
        in_peaks=0
    else
        # Reads overlapping consensus peaks (using -L provides fast bed overlap)
        in_peaks=$(samtools view -c -L $CONSENSUS -F 260 $bam)
        
        # Calculate FRiP
        frip=$(echo "$in_peaks $total" | awk '{printf "%.4f", $1/$2}')
    fi
    
    echo -e "${sample_name}\t${total}\t${in_peaks}\t${frip}" >> $OUTPUT_FILE
    echo "  Total: $total | In Peaks: $in_peaks | FRiP: $frip"
done

echo "Done! Results saved to $OUTPUT_FILE"
