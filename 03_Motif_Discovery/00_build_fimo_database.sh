#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=BUILD_DB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=04:00:00
#SBATCH -A lp_svbelleghem

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /data/leuven/354/vsc35429/miniconda3/envs/atac_motif_env

DB_DIR="/scratch/leuven/354/vsc35429/ATAC/motif_databases"
mkdir -p $DB_DIR
cd $DB_DIR

echo "==== 1. Downloading JASPAR 2026 FULL CORE ===="
wget -q "https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_meme.txt" -O JASPAR2026_CORE_all.meme
wget -q "https://mencius.uio.no/JASPAR/JASPAR_metadata/2026/ultimate_metadata_table_CORE.tsv" -O JASPAR2026_metadata.tsv

echo "==== 2. Downloading CIS-BP Daphnia pulex ===="
wget -q "https://cisbp.ccbr.utoronto.ca/tmp/Daphnia_pulex_2026_04_08_11:43_am.zip" -O cisbp_daphnia_pulex.zip
unzip -q -o cisbp_daphnia_pulex.zip -d cisbp_daphnia_raw/

echo "==== 3. Converting CIS-BP to MEME ===="
# cisbp2meme extracts PWMs and ties them to TF names
cisbp2meme -cisbp cisbp_daphnia_raw/pwms -info cisbp_daphnia_raw/TF_Information.txt > cisbp_daphnia.meme 2>/dev/null

echo "==== 4. Combining databases ===="
cat JASPAR2026_CORE_all.meme cisbp_daphnia.meme > combined_raw.meme
echo "Combined total motifs: $(grep -c "^MOTIF" combined_raw.meme)"

echo "==== 5. TOMTOM Redundancy Clustering ===="
tomtom -oc tomtom_selfcomp/ -thresh 0.05 -min-overlap 5 -dist pearson -p 8 combined_raw.meme combined_raw.meme 2> tomtom_selfcomp/tomtom.log

echo "==== 6. Parsing TOMTOM Redundancies (Python) ===="
python3 << 'PYEOF'
import pandas as pd
import os

try:
    tomtom = pd.read_csv("tomtom_selfcomp/tomtom.tsv", sep="\t", comment="#")
    tomtom = tomtom[tomtom["Query_ID"] != tomtom["Target_ID"]]
    similar = tomtom[tomtom["q-value"] < 0.05]
    
    seen = set()
    redundant = set()
    for _, row in similar.iterrows():
        q, t = row["Query_ID"], row["Target_ID"]
        if q not in redundant:
            seen.add(q)
            redundant.add(t)
            
    print(f"Redundant motifs to remove: {len(redundant)}")
    with open("redundant_motif_ids.txt", "w") as f:
        for m in redundant:
            f.write(str(m) + "\n")
except Exception as e:
    print(f"Python failed (check if pandas is installed): {e}")
    os.system("touch redundant_motif_ids.txt")
PYEOF

echo "==== 7. Extracting Non-Redundant Subset & HOMER Conversion (R) ===="
Rscript - << 'REOF'
# Make sure universalmotif is installed in your R environment!
# install.packages("BiocManager"); BiocManager::install("universalmotif")

if (!requireNamespace("universalmotif", quietly = TRUE)) {
    cat("ERROR: universalmotif package not installed in this R environment.\n")
    cat("Please install it to run the final filtration and HOMER conversion.\n")
    quit(status=1)
}
library(universalmotif)

motifs <- read_meme("combined_raw.meme")
if (file.exists("redundant_motif_ids.txt")) {
    redundant <- readLines("redundant_motif_ids.txt")
    if (length(redundant) > 0) {
        # Remove redundant motifs based on their MEME ID
        motifs <- motifs[!sapply(motifs, function(m) m@name %in% redundant)]
    }
}

write_meme(motifs, "combined_nonredundant.meme")
write_homer(motifs, "combined_nonredundant_homer.motif")
system("awk -F\\t\ \{if($0 ~ /^>/){ print $1 \"\\t\" $2 \"\\t\" $3 * 0.69314718 } else { print $0 }}\ combined_nonredundant_homer.motif > fixed_homer.motif")
cat("SUCCESS! Saved", length(motifs), "non-redundant motifs.\n")
REOF

echo "==== Done ===="
echo "Final MEME DB:  $DB_DIR/combined_nonredundant.meme"
echo "Final HOMER DB: $DB_DIR/combined_nonredundant_homer.motif"
