## Directory structure
```
|-- code
|-- sauron
|   |-- ...
|-- data
|   |-- preprocessing_scripts
|   |   |-- 02b_scale_integrate.R
|   |   |-- remove_cells.R
|   |-- gene_lists
|   |   |-- main_cell_types.csv
|   |   |-- bcell_types.csv
|   |   |-- bcell_types_germsub.csv
|   |   |-- bcell_types_germsub_zonesub.csv
|   |-- metadata.csv
|-- cellranger
|   |-- SampleID_1_14_mar19
|   |   |-- barcodes.tsv.gz
|   |   |-- features.tsv.gz
|   |   |-- matrix.mtx.gz
|   |-- ...
|-- salmon
|   |-- SampleID_1_14_mar19
|   |   |-- _quant_res
|   |-- ...
|-- velocyto
|   |-- SampleID_1_14_mar19.loom
|   |-- ...
```

## Download FASTQ files
Raw sequencing reads were downloaded from [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9491).

## Download Sauron workflow
Sauron pipeline for preprocessing was downloaded into ./sauron from[github](https://github.com/NBISweden/sauron/tree/master)

## Download Seurat scripts
Scripts for preprocessing were downloaded into ./data/preprocessing_scripts from [github](https://github.com/angelettilab/scMouseBcellFlu/tree/master/scripts/scRNAseq_pipeline).

## Download gene lists
Required gene lists for preprocessing was downloaded from [github](https://github.com/angelettilab/scMouseBcellFlu/tree/master/data/gene_lists).

## Download metadata
Required metadata for preprocessing was downloaded from [github](https://github.com/angelettilab/scMouseBcellFlu/blob/master/data/metadata.csv).

## Run Cell Ranger
cellranger count was run on each sample separately, folder structure compatible with published
code is shown as an example for sample "SampleID_1_14_mar19".

## Run salmon
salmon/alevin-fry was run on each sample saperately, folder structure compatible with published
code is shown as an example for sample "SampleID_1_14_mar19".

## Run velocyto
Velocyto was run on each sample separately, folder structure compatible with published
code is shown as an example for sample "SampleID_1_14_mar19".
