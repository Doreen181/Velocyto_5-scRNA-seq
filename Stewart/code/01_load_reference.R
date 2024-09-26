#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# convert files from scRNAseq dataset by Stewart et al. (doi: 10.3389/fimmu.2021.602539)
# to load them into python
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(Seurat)
library(SeuratData)
library(SeuratDisk)
# ==============================================================================
# load reference dataset
seurat <- readRDS(file = "../data/scPure2_HB6_UMAP3D.rds")
seurat <- UpdateSeuratObject(object = seurat)

# convert to anndata
SaveH5Seurat(seurat, filename = "../data/stewart_reference.h5Seurat")
Convert("../data/stewart_reference.h5Seurat", dest = "../data/stewart_reference.h5ad")

