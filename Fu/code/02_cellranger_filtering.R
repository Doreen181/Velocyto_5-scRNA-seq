#preprocessing script from Fu et al.
# (https://github.com/princello/scRNA-seq-TRM-paper/blob/main/SC_processing.R)

# get Seurat version 4
#remotes::install_version("Seurat", "4.0.0")

library(clustree)
library(dplyr)
library(Seurat)
library(DESeq2)
library(cluster)
suppressPackageStartupMessages({
    library(rlang)
})
library(ggplot2)
library(grid)
library(SignacX)
library(future)

#load data

mj001_data <- Read10X(data.dir = "../cellranger/MJ001") #1
mj002_data <- Read10X(data.dir = "../cellranger/MJ002") #2
mj003_data <- Read10X(data.dir = "../cellranger/MJ003") #3
mj005_data <- Read10X(data.dir = "../cellranger/MJ005") #4
mj006_data <- Read10X(data.dir = "../cellranger/MJ006") #5
mj007_data <- Read10X(data.dir = "../cellranger/MJ007") #12
mj008_data <- Read10X(data.dir = "../cellranger/MJ008") #6
mj009_data <- Read10X(data.dir = "../cellranger/MJ009") #7

mj016_data <- Read10X(data.dir = "../cellranger/fu/MJ016") #8
mj017_data <- Read10X(data.dir = "../cellranger/fu/MJ017") #9
mj018_data <- Read10X(data.dir = "../cellranger/fu/MJ018") #10
mj019_data <- Read10X(data.dir = "../cellranger/fu/MJ019") #11

#downsample size
#downsize  = 6400
# --> use cells from paper instead

# get reference from paper
ref <- readRDS(file = "../data/GSE252994_integrated8.rpca0919.rds")

######################################################################################

#processing individual data


drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj001_data,), value = FALSE)
mj001_count.data <- mj001_data[-drop.index, ]
mj001_me <- CreateSeuratObject(counts = mj001_count.data, project = "Pt15_POD1194")
mj001_me$pos <- "Pt15_POD1194"
mj001_me$mt <- PercentageFeatureSet(mj001_me, pattern = "^MT-")
mj001_me_Q1 = quantile(mj001_me$nFeature_RNA)[2]
mj001_me_Q3 = quantile(mj001_me$nFeature_RNA)[4]
mj001_me_interQ = (mj001_me_Q3-mj001_me_Q1)*1.5
mj001_me_upperb = mj001_me_Q3 + mj001_me_interQ
mj001_me_lowerb = mj001_me_Q1 - mj001_me_interQ

##########################################Mj002##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj002_data,), value = FALSE)
mj002_count.data <- mj002_data[-drop.index, ]
mj002_me <- CreateSeuratObject(counts = mj002_count.data, project = "Pt13_POD1032_IEL")
mj002_me$pos <- "Pt13_POD1032_IEL"
mj002_me$mt <- PercentageFeatureSet(mj002_me, pattern = "^MT-")
mj002_me_Q1 = quantile(mj002_me$nFeature_RNA)[2]
mj002_me_Q3 = quantile(mj002_me$nFeature_RNA)[4]
mj002_me_interQ = (mj002_me_Q3-mj002_me_Q1)*1.5
mj002_me_upperb = mj002_me_Q3 + mj002_me_interQ
mj002_me_lowerb = mj002_me_Q1 - mj002_me_interQ

##########################################Mj003##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj003_data,), value = FALSE)
mj003_count.data <- mj003_data[-drop.index, ]
mj003_me <- CreateSeuratObject(counts = mj003_count.data, project = "Pt13_POD1032_LPL")
mj003_me$pos <- "Pt13_POD1032_LPL"
mj003_me$mt <- PercentageFeatureSet(mj003_me, pattern = "^MT-")
mj003_me_Q1 = quantile(mj003_me$nFeature_RNA)[2]
mj003_me_Q3 = quantile(mj003_me$nFeature_RNA)[4]
mj003_me_interQ = (mj003_me_Q3-mj003_me_Q1)*1.5
mj003_me_upperb = mj003_me_Q3 + mj003_me_interQ
mj003_me_lowerb = mj003_me_Q1 - mj003_me_interQ

##########################################Mj005##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj005_data,), value = FALSE)
mj005_count.data <- mj005_data[-drop.index, ]
mj005_me <- CreateSeuratObject(counts = mj005_count.data, project = "Pt14_POD1764")
mj005_me$pos <- "Pt14_POD1764"
mj005_me$mt <- PercentageFeatureSet(mj005_me, pattern = "^MT-")
mj005_me_Q1 = quantile(mj005_me$nFeature_RNA)[2]
mj005_me_Q3 = quantile(mj005_me$nFeature_RNA)[4]
mj005_me_interQ = (mj005_me_Q3-mj005_me_Q1)*1.5
mj005_me_upperb = mj005_me_Q3 + mj005_me_interQ
mj005_me_lowerb = mj005_me_Q1 - mj005_me_interQ

##########################################Mj006##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj006_data,), value = FALSE)
mj006_count.data <- mj006_data[-drop.index, ]
mj006_me <- CreateSeuratObject(counts = mj006_count.data, project = "Pt21_POD626")
mj006_me$pos <- "Pt21_POD626"
mj006_me$mt <- PercentageFeatureSet(mj006_me, pattern = "^MT-")
mj006_me_Q1 = quantile(mj006_me$nFeature_RNA)[2]
mj006_me_Q3 = quantile(mj006_me$nFeature_RNA)[4]
mj006_me_interQ = (mj006_me_Q3-mj006_me_Q1)*1.5
mj006_me_upperb = mj006_me_Q3 + mj006_me_interQ
mj006_me_lowerb = mj006_me_Q1 - mj006_me_interQ

##########################################Mj007##########################################

mj007_me <- CreateSeuratObject(counts = mj007_data, project = "D251")
mj007_me <- RenameCells(mj007_me, new.names = paste0(colnames(x = mj007_me), "_12"))
mj007_me$clonotype <- 'Un'
mj007_me$pos <- "D251"
mj007_me$pre <- "Un"
mj007_me$post <- "Un"
mj007_me$pre_post <- "Un; Un"
mj007_me$mt <- PercentageFeatureSet(mj007_me, pattern = "^MT-")
mj007_me_Q1 = quantile(mj007_me$nFeature_RNA)[2]
mj007_me_Q3 = quantile(mj007_me$nFeature_RNA)[4]
mj007_me_interQ = (mj007_me_Q3-mj007_me_Q1)*1.5
mj007_me_upperb = mj007_me_Q3 + mj007_me_interQ
mj007_me_lowerb = mj007_me_Q1 - mj007_me_interQ
mj007_me <- subset(mj007_me, subset = nFeature_RNA > mj007_me_lowerb & nFeature_RNA < mj007_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "D251")))

##########################################Mj008##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj008_data,), value = FALSE)
mj008_count.data <- mj008_data[-drop.index, ]
mj008_me <- CreateSeuratObject(counts = mj008_count.data, project = "Pt04_POD1606_IEL")
mj008_me$pos <- "Pt04_POD1606_IEL"
mj008_me$mt <- PercentageFeatureSet(mj008_me, pattern = "^MT-")
mj008_me_Q1 = quantile(mj008_me$nFeature_RNA)[2]
mj008_me_Q3 = quantile(mj008_me$nFeature_RNA)[4]
mj008_me_interQ = (mj008_me_Q3-mj008_me_Q1)*1.5
mj008_me_upperb = mj008_me_Q3 + mj008_me_interQ
mj008_me_lowerb = mj008_me_Q1 - mj008_me_interQ

##########################################Mj009##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj009_data,), value = FALSE)
mj009_count.data <- mj009_data[-drop.index, ]
mj009_me <- CreateSeuratObject(counts = mj009_count.data, project = "Pt04_POD1606_LPL")
mj009_me$pos <- "Pt04_POD1606_LPL"
mj009_me$mt <- PercentageFeatureSet(mj009_me, pattern = "^MT-")
mj009_me_Q1 = quantile(mj009_me$nFeature_RNA)[2]
mj009_me_Q3 = quantile(mj009_me$nFeature_RNA)[4]
mj009_me_interQ = (mj009_me_Q3-mj009_me_Q1)*1.5
mj009_me_upperb = mj009_me_Q3 + mj009_me_interQ
mj009_me_lowerb = mj009_me_Q1 - mj009_me_interQ

##########################################Mj016##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj016_data,), value = FALSE)
mj016_count.data <- mj016_data[-drop.index, ]
mj016_me <- CreateSeuratObject(counts = mj016_count.data, project = "Pt16_POD1004_IEL")
mj016_me$pos <- "Pt16_POD1004_IEL"
mj016_me$mt <- PercentageFeatureSet(mj016_me, pattern = "^MT-")
mj016_me_Q1 = quantile(mj016_me$nFeature_RNA)[2]
mj016_me_Q3 = quantile(mj016_me$nFeature_RNA)[4]
mj016_me_interQ = (mj016_me_Q3-mj016_me_Q1)*1.5
mj016_me_upperb = mj016_me_Q3 + mj016_me_interQ
mj016_me_lowerb = mj016_me_Q1 - mj016_me_interQ

##########################################Mj017##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj017_data,), value = FALSE)
mj017_count.data <- mj017_data[-drop.index, ]
mj017_me <- CreateSeuratObject(counts = mj017_count.data, project = "Pt16_POD1004_LPL")
mj017_me$pos <- "Pt16_POD1004_LPL"
mj017_me$mt <- PercentageFeatureSet(mj017_me, pattern = "^MT-")
mj017_me_Q1 = quantile(mj017_me$nFeature_RNA)[2]
mj017_me_Q3 = quantile(mj017_me$nFeature_RNA)[4]
mj017_me_interQ = (mj017_me_Q3-mj017_me_Q1)*1.5
mj017_me_upperb = mj017_me_Q3 + mj017_me_interQ
mj017_me_lowerb = mj017_me_Q1 - mj017_me_interQ

##########################################Mj018##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj018_data,), value = FALSE)
mj018_count.data <- mj018_data[-drop.index, ]
mj018_me <- CreateSeuratObject(counts = mj018_count.data, project = "Pt21_POD1145_IEL")
mj018_me$pos <- "Pt21_POD1145_IEL"
mj018_me$mt <- PercentageFeatureSet(mj018_me, pattern = "^MT-")
mj018_me_Q1 = quantile(mj018_me$nFeature_RNA)[2]
mj018_me_Q3 = quantile(mj018_me$nFeature_RNA)[4]
mj018_me_interQ = (mj018_me_Q3-mj018_me_Q1)*1.5
mj018_me_upperb = mj018_me_Q3 + mj018_me_interQ
mj018_me_lowerb = mj018_me_Q1 - mj018_me_interQ

##########################################Mj019##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj019_data,), value = FALSE)
mj019_count.data <- mj019_data[-drop.index, ]
mj019_me <- CreateSeuratObject(counts = mj019_count.data, project = "Pt21_POD1145_LPL")
mj019_me$pos <- "Pt21_POD1145_LPL"
mj019_me$mt <- PercentageFeatureSet(mj019_me, pattern = "^MT-")
mj019_me_Q1 = quantile(mj019_me$nFeature_RNA)[2]
mj019_me_Q3 = quantile(mj019_me$nFeature_RNA)[4]
mj019_me_interQ = (mj019_me_Q3-mj019_me_Q1)*1.5
mj019_me_upperb = mj019_me_Q3 + mj019_me_interQ
mj019_me_lowerb = mj019_me_Q1 - mj019_me_interQ


#category meta data function
cat <- function(info) {
    categories <- info[,7]
    categories <- as.character(categories)
    categories[info[,5] == 'CD4_HvG' & info[,6] == "CD4_H'vG"] <- "Persistent HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "CD8_H'vG"] <- "Persistent HvG"
    categories[info[,5] == 'CD4_HvG' & info[,6] == "CD4_nonH'vG"] <- "Tolerant HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "CD8_nonH'vG"] <- "Tolerant HvG"
    categories[info[,5] == 'CD4_HvG' & info[,6] == "Un"] <- "Missing HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "Un"] <- "Missing HvG"
    categories[info[,5] == 'CD4_nonHvG' & info[,6] == "CD4_H'vG"] <- "Acquired H'vG"
    categories[info[,5] == 'CD8_nonHvG' & info[,6] == "CD8_H'vG"] <- "Acquired H'vG"
    categories[info[,5] == 'Un' & info[,6] == "CD4_H'vG"] <- "De novo H'vG"
    categories[info[,5] == 'Un' & info[,6] == "CD8_H'vG"] <- "De novo H'vG"
    categories[info[,5] == 'CD4_nonHvG' & info[,6] == "CD4_nonH'vG"] <- "Persistent nonHvG"
    categories[info[,5] == 'CD8_nonHvG' & info[,6] == "CD8_nonH'vG"] <- "Persistent nonHvG"
    categories[categories == 'Un; Un'] <- "Others"
    categories <- as.factor(categories)
    return(categories)
}

# rename barcodes to match reference
mj001_me <- RenameCells(mj001_me, new.names = paste0(colnames(x = mj001_me), "_1"))
mj002_me <- RenameCells(mj002_me, new.names = paste0(colnames(x = mj002_me), "_2"))
mj003_me <- RenameCells(mj003_me, new.names = paste0(colnames(x = mj003_me), "_3"))
mj005_me <- RenameCells(mj005_me, new.names = paste0(colnames(x = mj005_me), "_4"))
mj006_me <- RenameCells(mj006_me, new.names = paste0(colnames(x = mj006_me), "_5"))
mj008_me <- RenameCells(mj008_me, new.names = paste0(colnames(x = mj008_me), "_6"))
mj009_me <- RenameCells(mj009_me, new.names = paste0(colnames(x = mj009_me), "_7"))
mj016_me <- RenameCells(mj016_me, new.names = paste0(colnames(x = mj016_me), "_8"))
mj017_me <- RenameCells(mj017_me, new.names = paste0(colnames(x = mj017_me), "_9"))
mj018_me <- RenameCells(mj018_me, new.names = paste0(colnames(x = mj018_me), "_10"))
mj019_me <- RenameCells(mj019_me, new.names = paste0(colnames(x = mj019_me), "_11"))

##########################################DownSample##########################################

info = read.csv('../data/cell/mj005.csv')
mj005_me$clonotype <- info[,4]
mj005_me$pre <- info[,5]
mj005_me$post <- info[,6]
mj005_me$pre_post <- info[,7]
mj005_me$categories <- cat(info)

mj005_me <- subset(mj005_me, subset = nFeature_RNA > mj005_me_lowerb & nFeature_RNA < mj005_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt14_POD1764")))
Idents(mj005_me) <- 'categories'
mj005_me_hvg <- subset(mj005_me,idents = c("Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
mj005_me_un <- subset(mj005_me,idents = 'Others')
mj005_me_un_sub <- subset(mj005_me_un)
mj005_me <- merge(mj005_me_hvg, mj005_me_un_sub,  project = "Pt14_POD1764")
Idents(mj005_me) <- 'orig.ident'



#Idents(mj005_me) <- 'categories'
#mj005_me <- subset(mj005_me, subset = nFeature_RNA > mj005_me_lowerb & nFeature_RNA < mj005_me_upperb  & mt < 15)
#mj005_me_hvg <- subset(mj005_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj005_me_un <- subset(mj005_me,idents = 'Others')
#mj005_me_un_sub <- subset(mj005_me_un, downsample = downsize-length(mj005_me_hvg$orig.ident))
#mj005_me <- merge(mj005_me_hvg, mj005_me_un_sub,  project = "Pt14_POD1764")
#Idents(mj005_me) <- 'orig.ident'
#mj005_me <- NormalizeData(mj005_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj005_me <- FindVariableFeatures(mj005_me, selection.method = "vst", nfeatures = 5000)
#mj005_me <- SCTransform(mj005_me, verbose = FALSE)



info = read.csv('../data/cell/mj006.csv')
mj006_me$clonotype <- info[,4]
mj006_me$pre <- info[,5]
mj006_me$post <- info[,6]
mj006_me$pre_post <- info[,7]
mj006_me$categories <- cat(info)
mj006_me <- subset(mj006_me, subset = nFeature_RNA > mj006_me_lowerb & nFeature_RNA < mj006_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt21_POD626")))

Idents(mj006_me) <- 'categories'
mj006_me_hvg <- subset(mj006_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
mj006_me_un <- subset(mj006_me,idents = 'Others')
mj006_me_un_sub <- subset(mj006_me_un)
mj006_me <- merge(mj006_me_hvg, mj006_me_un_sub,  project = "Pt21_POD626")
Idents(mj006_me) <- 'orig.ident'



#mj006_me <- subset(mj006_me, subset = nFeature_RNA > mj006_me_lowerb & nFeature_RNA < mj006_me_upperb  & mt < 15)
#mj006_me_hvg <- subset(mj006_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj006_me_un <- subset(mj006_me,idents = 'Others')
#mj006_me_un_sub <- subset(mj006_me_un, downsample = downsize-length(mj006_me_hvg$orig.ident))
#mj006_me <- merge(mj006_me_hvg, mj006_me_un_sub,  project = "Pt21_POD626")
#Idents(mj006_me) <- integrated8.rpca
#mj006_me <- NormalizeData(mj006_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj006_me <- FindVariableFeatures(mj006_me, selection.method = "vst", nfeatures = 5000)
#mj006_me <- SCTransform(mj006_me, verbose = FALSE)



info = read.csv('../data/cell/mj001.csv')
mj001_me$clonotype <- info[,4]
mj001_me$pre <- info[,5]
mj001_me$post <- info[,6]
mj001_me$pre_post <- info[,7]
mj001_me$categories <- cat(info)
mj001_me <- subset(mj001_me, subset = nFeature_RNA > mj001_me_lowerb & nFeature_RNA < mj001_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt15_POD1194")))

#mj001_me <- NormalizeData(mj001_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj001_me <- FindVariableFeatures(mj001_me, selection.method = "vst", nfeatures = 5000)
#mj001_me <- SCTransform(mj001_me, verbose = TRUE)

info = read.csv('../data/cell/mj002.csv')
mj002_me$clonotype <- info[,4]
mj002_me$pre <- info[,5]
mj002_me$post <- info[,6]
mj002_me$pre_post <- info[,7]
mj002_me$categories <- cat(info)
mj002_me <- subset(mj002_me, subset = nFeature_RNA > mj002_me_lowerb & nFeature_RNA < mj002_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt13_POD1032_IEL")))

#mj002_me <- NormalizeData(mj002_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj002_me <- FindVariableFeatures(mj002_me, selection.method = "vst", nfeatures = 5000)
#mj002_me <- SCTransform(mj002_me, verbose = FALSE)

info = read.csv('../data/cell/mj003.csv')
mj003_me$clonotype <- info[,4]
mj003_me$pre <- info[,5]
mj003_me$post <- info[,6]
mj003_me$pre_post <- info[,7]
mj003_me$categories <- cat(info)
mj003_me <- subset(mj003_me, subset = nFeature_RNA > mj003_me_lowerb & nFeature_RNA < mj003_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt13_POD1032_LPL")))


#mj003_me <- NormalizeData(mj003_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj003_me <- FindVariableFeatures(mj003_me, selection.method = "vst", nfeatures = 5000)
#mj003_me <- SCTransform(mj003_me, verbose = FALSE)

info = read.csv('../data/cell/mj008.csv')
mj008_me$clonotype <- info[,4]
mj008_me$pre <- info[,5]
mj008_me$post <- info[,6]
mj008_me$pre_post <- info[,7]
mj008_me$categories <- cat(info)
mj008_me <- subset(mj008_me, subset = nFeature_RNA > mj008_me_lowerb & nFeature_RNA < mj008_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt04_POD1606_IEL")))


#mj008_me <- NormalizeData(mj008_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj008_me <- FindVariableFeatures(mj008_me, selection.method = "vst", nfeatures = 5000)

#mj008_me <- SCTransform(mj008_me, verbose = FALSE)

info = read.csv('../data/cell/mj009.csv')
mj009_me$clonotype <- info[,4]
mj009_me$pre <- info[,5]
mj009_me$post <- info[,6]
mj009_me$pre_post <- info[,7]
mj009_me$categories <- cat(info)
mj009_me <- subset(mj009_me, subset = nFeature_RNA > mj009_me_lowerb & nFeature_RNA < mj009_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt04_POD1606_LPL")))




#Idents(mj009_me) <- 'categories'
#mj009_me <- subset(mj009_me, subset = nFeature_RNA > mj009_me_lowerb & nFeature_RNA < mj009_me_upperb  & mt < 15)
#mj009_me_hvg <- subset(mj009_me,idents = c("Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj009_me_un <- subset(mj009_me,idents = 'Others')
#mj009_me_un_sub <- subset(mj009_me_un, downsample = downsize-length(mj009_me_hvg$orig.ident))
#mj009_me <- merge(mj009_me_hvg, mj009_me_un_sub,  project = "Pt04_POD1606_LPL")
#Idents(mj009_me) <- 'orig.ident'
#mj009_me <- NormalizeData(mj009_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj009_me <- FindVariableFeatures(mj009_me, selection.method = "vst", nfeatures = 5000)
#mj009_me <- SCTransform(mj009_me, verbose = FALSE)



info = read.csv('../data/cell/mj016.csv')
mj016_me$clonotype <- info[,4]
mj016_me$pre <- info[,5]
mj016_me$post <- info[,6]
mj016_me$pre_post <- info[,7]
mj016_me$categories <- cat(info)
mj016_me <- subset(mj016_me, subset = nFeature_RNA > mj016_me_lowerb & nFeature_RNA < mj016_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt16_POD1004_IEL")))


info = read.csv('../data/cell/mj017.csv')
mj017_me$clonotype <- info[,4]
mj017_me$pre <- info[,5]
mj017_me$post <- info[,6]
mj017_me$pre_post <- info[,7]
mj017_me$categories <- cat(info)
mj017_me <- subset(mj017_me, subset = nFeature_RNA > mj017_me_lowerb & nFeature_RNA < mj017_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt16_POD1004_LPL")))


info = read.csv('../data/cell/mj018.csv')
mj018_me$clonotype <- info[,4]
mj018_me$pre <- info[,5]
mj018_me$post <- info[,6]
mj018_me$pre_post <- info[,7]
mj018_me$categories <- cat(info)
mj018_me <- subset(mj018_me, subset = nFeature_RNA > mj018_me_lowerb & nFeature_RNA < mj018_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt21_POD1145_IEL")))



info = read.csv('../data/cell/mj019.csv')
mj019_me$clonotype <- info[,4]
mj019_me$pre <- info[,5]
mj019_me$post <- info[,6]
mj019_me$pre_post <- info[,7]
mj019_me$categories <- cat(info)
mj019_me <- subset(mj019_me, subset = nFeature_RNA > mj019_me_lowerb & nFeature_RNA < mj019_me_upperb  & mt < 15,
                   cells = colnames(subset(ref, subset = orig.ident == "Pt21_POD1145_LPL")))

##########################################Integrating##########################################

samples <- list( mj001_me,mj002_me,mj003_me,mj005_me,mj006_me,mj008_me,mj009_me,mj016_me,mj017_me,mj018_me,mj019_me,mj007_me)

mj.list <- lapply(X = samples, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6400)
})

features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 6400)

mj.list <- lapply(X = mj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors.rpca <- FindIntegrationAnchors(object.list = mj.list, reduction = "rpca", anchor.features =features)


integrated8.rpca <- IntegrateData(anchorset = anchors.rpca,dims = 1:64)
rm(anchors.rpca)

saveRDS(integrated8.rpca,'../data/fu_cellranger_preprocessed.rds')