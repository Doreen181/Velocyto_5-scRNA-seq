#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Mathew: convert files from RDS to h5ad
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(Seurat)
library(SeuratData)
library(SeuratDisk)

# ==============================================================================
# load processed dataset: cellranger

cellranger <- readRDS(file = "./02_filtering_and_preprocessing/cellranger/06_cluster/seurat_object.rds")
cellranger <- UpdateSeuratObject(object = cellranger)

# convert to anndata
SaveH5Seurat(cellranger, filename = "../data/processed_seurat_cellranger.h5Seurat")
Convert("../data/processed_seurat_cellranger.h5Seurat", dest = "../data/processed_seurat_cellranger.h5ad")


# ==============================================================================
# load processed dataset: salmon

salmon <- readRDS(file = "./02_filtering_and_preprocessing/salmon/06_cluster/seurat_object.rds")
salmon <- UpdateSeuratObject(object = salmon)

# convert to anndata
SaveH5Seurat(salmon, filename = "../data/processed_seurat_salmon.h5Seurat")
Convert("../data/processed_seurat_salmon.h5Seurat", dest = "../data/processed_seurat_salmon.h5ad")