## Directory structure

|-- code \
|-- data \
|   |-- GSE252994_integrated8.rpca0919.rds \
|   |-- cell \
|   |   |-- mj001.csv \
|   |   |-- ... \
|   |-- metadata.csv \
|-- cellranger \
|   |-- MJ001 \
|   |   |-- barcodes.tsv.gz \
|   |   |-- features.tsv.gz \
|   |   |-- matrix.mtx.gz \
|   |-- ... \
|-- salmon \
|   |-- MJ001 \
|   |   |-- MJ001_quant_res \
|   |-- ... \
|-- velocyto \
|   |-- MJ001.loom \
|   |-- ... \


## Download BAM files
BAM files were downloaded from [SRA](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=SRP417027) and converted to FASTQ files by cellranger bamtofastq.

## Download Seurat scripts
Scripts for preprocessing were downloaded from [github](https://github.com/princello/scRNA-seq-TRM-paper/blob/main/SC_processing.R) and adjusted.

## Download reference dataset
GSE252994_integrated8.rpca0919.rds (processed dataset by Fu et al. was downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252994).

## Download metadata
Required metadata for preprocessing was downloaded into ./cell from [github](https://github.com/princello/scRNA-seq-TRM-paper/tree/main/cell).

## Run Cell Ranger
Cell Ranger count was run on each sample seperately, folder structure compatible with published
code is shown as an example for sample "MJ001".

## Run salmon
salmon/alevin-fry was run on each sample seperately, folder structure compatible with published
code is shown as an example for sample "MJ001".

## Run velocyto
Velocyto was run on each sample seperately, folder structure compatible with published
code is shown as an example for sample "MJ001".