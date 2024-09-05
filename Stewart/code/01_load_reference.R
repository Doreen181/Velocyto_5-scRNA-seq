#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# convert files from scRNAseq dataset by Stewart et al. (doi: 10.3389/fimmu.2021.602539)
# to load them into python
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(Seurat)
library(MuDataSeurat)

# ==============================================================================
# load reference dataset
seurat <- readRDS(file = "../data/scPure2_AllIntegrated.rds")
seurat <- UpdateSeuratObject(object = seurat)

# convert to anndata
MuDataSeurat::WriteH5AD(seurat, "../data/stewart_reference.h5ad", assay="RNA")
