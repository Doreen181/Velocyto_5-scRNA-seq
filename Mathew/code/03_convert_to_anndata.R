#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Mathew: convert files from RDS to h5ad
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(Seurat)
library(MuDataSeurat)

# ==============================================================================
# load processed dataset: cellranger

cellranger <- readRDS(file = "./02_filtering_and_preprocessing/cellranger/06_cluster/seurat_object.rds")
cellranger <- UpdateSeuratObject(object = cellranger)

# convert to anndata
MuDataSeurat::WriteH5AD(cellranger, "../data/processed_seurat_cellranger.h5ad", assay="RNA")

# ==============================================================================
# load processed dataset: salmon

salmon <- readRDS(file = "./02_filtering_and_preprocessing/salmon/06_cluster/seurat_object.rds")
salmon <- UpdateSeuratObject(object = salmon)

# convert to anndata
MuDataSeurat::WriteH5AD(salmon, "../data/processed_seurat_salmon.h5ad", assay="RNA")