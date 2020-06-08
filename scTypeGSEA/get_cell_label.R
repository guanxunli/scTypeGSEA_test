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
cluster_cell <- rep(NA, length(dta_seurat$celltype))
for (i in 1:length(celltype)){
  cluster_cell[which(dta_seurat$celltype == celltype[i])] = i - 1
}

dta_seurat <- doClustering(dta_seurat, cluster_cell = cluster_cell)
cluster_list <- getFC(dta_seurat)
cluster_celltype <- doGSEA(cluster_list = cluster_list)

write.table(x = data.frame("TRUE" = celltype, "method" = cluster_celltype), file = "scTypeGSEA_pbmc1.txt")

## data set 2

dta <- readRDS("dataset/pbmc2/pbmc2.rds")
dta_celltype <- readRDS("dataset/pbmc2/pbmc2_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype[, -1])
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
cluster_cell <- rep(NA, length(dta_seurat$celltype))
for (i in 1:length(celltype)){
  cluster_cell[which(dta_seurat$celltype == celltype[i])] = i - 1
}

dta_seurat <- doClustering(dta_seurat, cluster_cell = cluster_cell)
cluster_list <- getFC(dta_seurat)
cluster_celltype <- doGSEA(cluster_list = cluster_list)

write.table(x = data.frame("TRUE" = celltype, "method" = cluster_celltype), file = "scTypeGSEA_pbmc2.txt")

## data set 3

dta <- readRDS("dataset/GSE81861/GSE81861.rds")
dta_celltype <- readRDS("dataset/GSE81861/GSE81861_celltype.rds")
dta <- dta[, -133]
dta_celltype <- dta_celltype[-133]
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype, min.cells = 10, min.features = 1000)

celltype <- unique(dta_seurat$celltype)
cluster_cell <- rep(NA, length(dta_seurat$celltype))
for (i in 1:length(celltype)){
  cluster_cell[which(dta_seurat$celltype == celltype[i])] = i - 1
}

dta_seurat <- doClustering(dta_seurat, cluster_cell = cluster_cell)
cluster_list <- getFC(dta_seurat)
cluster_celltype <- doGSEA(cluster_list = cluster_list)

saveRDS(cluster_celltype, "scTypeGSEA/celltype_GSE81861.rds")


## data set 4

dta <- readRDS("dataset/GSE72056/GSE72056.rds")
dta_celltype <- readRDS("dataset/GSE72056/GSE72056_celltype.rds")
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
cluster_cell <- rep(NA, length(dta_seurat$celltype))
for (i in 1:length(celltype)){
  cluster_cell[which(dta_seurat$celltype == celltype[i])] = i - 1
}

dta_seurat <- doClustering(dta_seurat, cluster_cell = cluster_cell)
cluster_list <- getFC(dta_seurat)
cluster_celltype <- doGSEA(cluster_list = cluster_list)

write.table(x = data.frame("TRUE" = celltype, "method" = cluster_celltype), file = "scTypeGSEA_GSE72056.txt")

## data set 5

dta <- readRDS("dataset/E-MTAB-6149/E_MTAB-6149.rds.rds")
dta_celltype <- readRDS("dataset/E-MTAB-6149/E_MTAB-6149_celltype.rds.rds")
dta_celltype <- as.data.frame(dta_celltype)
rownames(dta_celltype) <- colnames(dta)
colnames(dta_celltype) <- "celltype"

dta_seurat <- scqc(dta, metadata = dta_celltype)

celltype <- unique(dta_seurat$celltype)
cluster_cell <- rep(NA, length(dta_seurat$celltype))
for (i in 1:length(celltype)){
  cluster_cell[which(dta_seurat$celltype == celltype[i])] = i - 1
}

dta_seurat <- doClustering(dta_seurat, cluster_cell = cluster_cell)
cluster_list <- getFC(dta_seurat)
cluster_celltype <- doGSEA(cluster_list = cluster_list)

write.table(x = data.frame("TRUE" = celltype, "method" = cluster_celltype), file = "scTypeGSEA_E-MTAB-6149.txt")
