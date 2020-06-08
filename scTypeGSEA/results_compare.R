library(Seurat)
library(fgsea)
library(scTypeGSEA)

## data set 1

dta <- readRDS("dataset/pbmc1/pbmc1.rds")
dta_celltype <- readRDS("dataset/pbmc1/pbmc1_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype[, -1])
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
celltype_list <- readRDS("scTypeGSEA/celltype_pbmc1.rds")

## data set 2

dta <- readRDS("dataset/pbmc2/pbmc2.rds")
dta_celltype <- readRDS("dataset/pbmc2/pbmc2_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype[, -1])
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
celltype_list <- readRDS("scTypeGSEA/celltype_pbmc2.rds")
  
## data set 3

dta <- readRDS("dataset/GSE81861/GSE81861.rds")
dta_celltype <- readRDS("dataset/GSE81861/GSE81861_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
celltype_list <- readRDS("scTypeGSEA/celltype_GSE81861.rds")



## data set 4

dta <- readRDS("dataset/GSE72056/GSE72056.rds")
dta_celltype <- readRDS("dataset/GSE72056/GSE72056_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
celltype_list <- readRDS("scTypeGSEA/celltype_GSE72056.rds")



## data set 5

dta <- readRDS("dataset/E-MTAB-6149/E_MTAB-6149.rds")
dta_celltype <- readRDS("dataset/E-MTAB-6149/E_MTAB-6149_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
celltype_list <- readRDS("scTypeGSEA/celltype_E_MTAB-6149.rds")
