## Directory structure

|-- code
|-- data
|   |-- scPure2_AllIntegrated.rds
|-- cellranger
|   |-- naive
|       |-- outs
|   |   |   |-- filtered_feature_bc_matrix
|   |   |-- velocyto
|   |   |   |-- naive.loom
|   |-- ...
|-- salmon
|   |-- naive
|   |   |-- naive_quant_res
|   |-- ...


## Download FASTQ files
Raw sequencing reads were downloaded from [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9544?query=E-MTAB-9544).

## Download reference dataset
scPure2_AllIntegrated.rds (processed dataset by Stewart et al. was downloaded from
[ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9544?query=E-MTAB-9544).

## Run Cell Ranger
cellranger count was run on each sample seperately, folder structure compatible with published
code is shown as an example for sample "naive".

## Run salmon
salmon/alevin-fry was run on each sample seperately, folder structure compatible with published
code is shown as an example for sample "naive".

## Run velocyto
Velocyto was run on each sample seperately, inside directory ./cellranger.
Folder structure compatible with published code is shown as an example for sample "naive"