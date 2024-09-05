#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fu: convert files from RDS to h5ad
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(Seurat)
library(MuData)
library(Matrix)

# ==============================================================================
# load reference from script
ref <- readRDS(file = "../data/GSE252994_integrated8.rpca0919.rds")

# save raw count matrix
counts_matrix <- GetAssayData(ref, assay='RNA', layer='counts')
writeMM(counts_matrix, file='../data/ref_raw_counts.mtx')

# save normalized count matrix
norm_matrix <- GetAssayData(ref, assay='RNA')
writeMM(norm, file='../data/ref_norm_counts.mtx')

# save gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='../data/ref_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

# save variable features
write.csv(ref@assays[["integrated"]]@var.features, file='../data/ref_hvg.csv', quote=F, row.names=F)

# add barcodes and UMAP to metadata
ref[["barcode"]] <- colnames(ref)
ref[["UMAP_1"]] <- ref@reductions[["umap"]]@cell.embeddings[,1]
ref[["UMAP_2"]] <- ref@reductions[["umap"]]@cell.embeddings[,2]

# save meta data
write.csv(ref@meta.data, file='../data/ref_metadata.csv', quote=F, row.names=F)

# save PCA
write.csv(ref@reductions[["pca"]]@cell.embeddings, file='../data/ref_pca.csv', quote=F, row.names=F)

# ==============================================================================
# load processed dataset: cellranger

cellranger <- readRDS(file = "../data/cellranger_preprocessed.rds")

# save raw count matrix
cellranger_counts_matrix <- GetAssayData(cellranger, assay='RNA', layer='counts')
writeMM(cellranger_counts_matrix, file='../data/cellranger_raw_counts.mtx')

# save normalized count matrix
cellranger_norm_matrix <- GetAssayData(cellranger, assay='RNA')
writeMM(cellranger_norm_matrix, file='../data/cellranger_norm_counts.mtx')

# save gene names
write.table(
  data.frame('gene'=rownames(cellranger_counts_matrix)),file='../data/cellranger_gene_names.csv',
  quote=F,row.names=F,col.names=F
)


# add barcodes to metadata
cellranger[["barcode"]] <- colnames(cellranger)

# save meta data
write.csv(cellranger@meta.data, file='../data/cellranger_metadata.csv', quote=F, row.names=F)

# ==============================================================================
# load processed dataset: salmon

salmon <- readRDS(file = "../data/salmon_preprocessed.rds")

# save raw count matrix
salmon_counts_matrix <- GetAssayData(salmon, assay='RNA', layer='counts')
writeMM(salmon_counts_matrix, file='../data/salmon_raw_counts.mtx')

# save normalized count matrix
salmon_norm_matrix <- GetAssayData(salmon, assay='RNA')
writeMM(salmon_norm_matrix, file='../data/salmon_norm_counts.mtx')

# save gene names
write.table(
  data.frame('gene'=rownames(salmon_counts_matrix)),file='../data/salmon_gene_names.csv',
  quote=F,row.names=F,col.names=F
)


# add barcodes to metadata
salmon[["barcode"]] <- colnames(salmon)

# save meta data
write.csv(salmon@meta.data, file='../data/salmon_metadata.csv', quote=F, row.names=F)