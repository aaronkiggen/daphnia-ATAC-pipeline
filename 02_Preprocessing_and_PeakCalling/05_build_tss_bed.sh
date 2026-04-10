#!/bin/bash
# Script to build a BED file containing Transcription Start Sites (TSS) from a GFF/GTF genome annotation.
# Usage: ./build_tss_bed.sh path/to/genome_annotation.gff3 [output_file.bed]

INPUT_GFF=$1
OUTPUT_BED=${2:-"/user/leuven/354/vsc35429/DATA/PhD/ATAC/Dmagna_TSS.bed"}

if [ -z "$INPUT_GFF" ]; then
    echo "Usage: bash build_tss_bed.sh <annotation.gff3> [output.bed]"
    echo ""
    echo "Please provide your Daphnia magna genome annotation (GFF/GFF3 or GTF file)."
    exit 1
fi

if [ ! -f "$INPUT_GFF" ]; then
    echo "Error: Cannot find file '$INPUT_GFF'."
    exit 1
fi

echo "Extracting TSS positions from '$INPUT_GFF'..."

# Parse the annotation file.
# 1. Filters for the feature type "gene" (ignores exons/mRNA/CDS).
# 2. For (+) strand genes: TSS = Start Coordinate
# 3. For (-) strand genes: TSS = End Coordinate
# 4. BED is 0-based, so we subtract 1 from the genomic coordinate.

zcat -f "$INPUT_GFF" | awk -F'\t' 'BEGIN {OFS="\t"} 
    !/^#/ && $3 == "gene" {
        chrom = $1
        start_gff = $4
        end_gff = $5
        strand = $7
        
        # Try to extract the gene ID or Name to label it in the BED file
        id = "."
        n=split($9, attr, ";")
        for(i=1; i<=n; i++) {
            if(attr[i] ~ /^ID=/) {
                id = substr(attr[i], 4)
                break
            } else if(attr[i] ~ /gene_id "/) {
                split(attr[i], arr, "\"")
                id = arr[2]
                break
            }
        }
        
        # Identify the precise 1bp TSS based on strand orientation
        if (strand == "+") {
            tss_start = start_gff - 1
            tss_end = start_gff
        } else if (strand == "-") {
            tss_start = end_gff - 1
            tss_end = end_gff
        } else {
            next
        }
        
        print chrom, tss_start, tss_end, id, ".", strand
    }' > "$OUTPUT_BED"

echo "Success! The TSS positions were saved to: $OUTPUT_BED"
echo "Total TSS regions found: $(wc -l < "$OUTPUT_BED")"
