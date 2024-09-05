#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nCount_RNA,nFeature_RNA,perc_mito,CC.Diff'  # regress CC.Diff (S - G2M) instead of S and G2M separately
main='~/Mathew/code/02_filtering_and_preprocessing/salmon'
script_path='~/Mathew/sauron/scripts'

cd $main

##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
conda activate Sauron.v1  # macOS (local)
# source activate Sauron.v1  # linux/unix (cluster)



#####################
### LOAD DATASETS ###
#####################
mkdir $main/'01_qc'
Rscript $script_path/00_load_data_salmon.R \
--input_path '../../salmon' \
--dataset_metadata_path '~/Mathew/data/metadata.csv' \
--assay 'RNA' \
--output_path $main/'01_qc' \
2>&1 | tee $main/'01_qc/00_load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $main/'01_qc/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes 'Ighd, Ighm, Ighg1, Ighg2c, Ighg2b, Ighg3, Igha, Ighe' \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--pct_mito_range '0,25' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $main/'01_qc' \
2>&1 | tee $main/'01_qc/01_qc_log.txt'



###################################################
### RUN INTEGRATION WITH PRE-SCALING/REGRESSION ###
###################################################
mkdir $main/'02_cluster'
Rscript '~/Mathew/data/preprocessing_scripts/02b_scale_integrate.R' \
--Seurat_object_path $main/'01_qc/filt_seurat_object.rds' \
--sauron_script_path $script_path \
--regress $var_to_regress \
--integration_method 'mnn,dataset' \
--var_genes 'seurat' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $main/'02_cluster' \
2>&1 | tee $main/'02_cluster/integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'02_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42;'\
'umap10, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42' \
--pre_dim_reduct 'mnn' \
--cluster_use 'all' \
--cluster_method 'louvain' \
--assay 'RNA' \
--output_path $main/'02_cluster' \
2>&1 | tee $main/'02_cluster/dr_and_cluster_log.txt'



###########################
## CELL TYPE PREDICTION ###
###########################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'02_cluster/seurat_object.rds' \
--marker_lists '~/Mathew/data/gene_lists/main_cell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_0.65' \
--output_path $main/'02_cluster/cell_type_prediction' \
2>&1 | tee $main/'02_cluster/cell_type_prediction_log.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
mkdir $main/'03_diff_expr'
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'02_cluster/seurat_object.rds' \
--clustering_use 'louvain_0.65' \
--metadata_use 'infection,organ' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'03_diff_expr' \
2>&1 | tee $main/'03_diff_expr/diff_expr_log.txt'



#######################################
## REMOVE NON-B-CELLS FROM THE DATA ###
#######################################
# Cells that are NOT predicted as B-cells will be removed from the data, and the pipeline
# re-run starting from the initial data loading step.
# The following code will extract the barcodes of cells that were in the specified clusters OR cells
# that were not predicted by the correlation analysis to be B-cells.
mkdir $main/'04_remove_cells'
Rscript '~/Mathew/data/preprocessing_scripts/remove_cells.R' \
--Seurat_object_path $main/'02_cluster/seurat_object.rds' \
--remove 'louvain_0.65,0,11,12,13,15' \
--keep 'cell_pred_correlation_main_cell_types,B_cell' \
--combine_method 'union' \
--output_type 'barcodes' \
--output_path $main/'04_remove_cells' \
2>&1 | tee $main/'04_remove_cells/remove_cells_log.txt'



###################################################
## RE-LOAD DATASETS, EXCLUDING THE NON-B-BCELLS ###
###################################################
mkdir $main/'05_qc'
Rscript $script_path/00_load_data_salmon.R \
--input_path '~/Mathew/salmon' \
--dataset_metadata_path '~/Mathew/data/metadata.csv' \
--assay 'RNA' \
--remove_cells $main/'04_remove_cells/remove_cell_barcodes.txt' \
--output_path $main/'05_qc' \
2>&1 | tee $main/'05_qc/load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $main/'05_qc/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes 'Ighd, Ighm, Ighg1, Ighg2c, Ighg2b, Ighg3, Igha, Ighe' \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--pct_mito_range '0,25' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $main/'05_qc' \
2>&1 | tee $main/'05_qc/qc_log.txt'



###################################################
### RUN INTEGRATION WITH PRE-SCALING/REGRESSION ###
###################################################
mkdir $main/'06_cluster'
Rscript '~/Mathew/data/preprocessing_scripts/02b_scale_integrate.R' \
--Seurat_object_path $main/'05_qc/filt_seurat_object.rds' \
--sauron_script_path $script_path \
--regress $var_to_regress \
--integration_method 'mnn,dataset' \
--var_genes 'seurat' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $main/'06_cluster' \
2>&1 | tee $main/'06_cluster/integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'06_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42;'\
'umap10, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42' \
--pre_dim_reduct 'mnn' \
--cluster_use 'all' \
--cluster_method 'louvain,HC' \
--assay 'RNA' \
--output_path $main/'06_cluster' \
2>&1 | tee $main/'06_cluster/dr_and_cluster_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'06_cluster/seurat_object.rds' \
--marker_lists '~/Mathew/data/gene_lists/main_cell_types.csv','~/data/Mathew/gene_lists/bcell_types.csv','~/data/Mathew/gene_lists/bcell_types_germsub.csv','~/data/Mathew/gene_lists/bcell_types_germsub_zonesub.csv' \
--clustering_use 'HC_16' \
--assay 'RNA' \
--output_path $main/'06_cluster/cell_type_prediction' \
2>&1 | tee $main/'06_cluster/cell_subtype_prediction_log.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'06_cluster/seurat_object.rds' \
--clustering_use 'HC_16' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'07_diff_expr' \
2>&1 | tee $main/'07_diff_expr/diff_expr_log.txt'




conda deactivate
