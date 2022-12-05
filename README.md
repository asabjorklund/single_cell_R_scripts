## Collection of custom scripts used in Single Cell analyses

Source a script with devtools::source_url("https://raw.github.com/single_cell_R_scripts/<script_name>.R")




* `misc_seurat.R` - Contains functions for running seurat pipeline, DEG testing etc.
* `top_expressed_genes.R` - Boxplot with top expressed genes from a Seurat object
* `overlap_phyper_v2.R` - Heatmap with overlap of two lists, with coloring by phyper p-value.
* `vdj_analysis.R` - Parsing of cellranger VDJ data
* `clonal_diversity.R` - Some functions for calculating diversity scores for clonal distributions.
