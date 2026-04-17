import sys

def mod_file(path, replacements):
    with open(path, "r") as f:
        content = f.read()
    for old, new in replacements:
        content = content.replace(old, new)
    with open(path, "w") as f:
        f.write(content)

p4_reps = [
    ('paste0("GRN_Biological_Edges_All_"', 'paste0("02_Filtered_Bio_Network_Edges_All_"'),
    ('paste0("Cytoscape_Collapsed_Edges_All_DEG_"', 'paste0("04_Cyto_Collapsed_Edges_All_DEGTFs_"'),
    ('paste0("Cytoscape_Collapsed_Nodes_All_DEG_"', 'paste0("04_Cyto_Collapsed_Nodes_All_DEGTFs_"')
]
mod_file("/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/08_grn_phase4_collapse.R", p4_reps)

p5_reps = [
    ('paste0("GRN_Biological_Edges_All_"', 'paste0("02_Filtered_Bio_Network_Edges_All_"'),
    ('paste0("Validation_Top20_TFs_"', 'paste0("05_Ortholog_Validation_Top20_TFs_"'),
    ('Wrote Validation_Top20_TFs_', 'Wrote 05_Ortholog_Validation_Top20_TFs_')
]
mod_file("/user/leuven/354/vsc35429/DATA/PhD/ATAC/ATAC_Pipeline_GitHub/04_GRN_Construction/09_grn_phase5_validation.R", p5_reps)

print("Phase 4 and 5 Updated!")
