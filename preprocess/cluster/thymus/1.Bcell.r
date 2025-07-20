library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
thymus_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/thymus_step_1_processed.rds')

library(harmony)

Bcell <- thymus_step_1[, Idents(thymus_step_1) %in% c('8', '32')]
Bcell <- NormalizeData(Bcell, verbose = T)
Bcell <- FindVariableFeatures(Bcell, selection.method = "vst", nfeatures = 2000)
Bcell <- ScaleData(Bcell)
Bcell <- RunPCA(Bcell, npcs = 50, verBcellose = FALSE)
ElbowPlot(Bcell, ndims = 50)

Bcell <- RunTSNE(Bcell, reduction = "pca", dims = 1:20)
Bcell <- RunUMAP(Bcell, reduction = "pca", dims = 1:20)
Bcell <- FindNeighbors(Bcell, reduction = "pca", dims = 1:20)



Bcell <- FindClusters(Bcell, resolution = 0.4)
DimPlot(Bcell, reduction = 'tsne', label = T)
DimPlot(Bcell, reduction = 'tsne', label = T, group.by = 'batch_2')
FeaturePlot(Bcell, reduction = 'tsne', pt.size = 0.5,
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

options(repr.plot.width=12, repr.plot.height=12)
FeaturePlot(Bcell, reduction = 'tsne', pt.size = 0.5,
            features = c('CD3E', 'CD3D', 'MS4A1','IGKC'))

###removed T-cell doublets
Bcell_2 <- Bcell[, !Idents(Bcell) %in% c('6', '8')]
Bcell_2 <- NormalizeData(Bcell_2, verbose = T)
Bcell_2 <- FindVariableFeatures(Bcell_2, selection.method = "vst", nfeatures = 2000)
Bcell_2 <- ScaleData(Bcell_2, vars.to.regress = c("batch_2"))
Bcell_2 <- RunPCA(Bcell_2, npcs = 50, verBcell_2ose = FALSE)
ElbowPlot(Bcell_2, ndims = 50)

Bcell_2 <- RunTSNE(Bcell_2, reduction = "pca", dims = 1:20)
Bcell_2 <- RunUMAP(Bcell_2, reduction = "pca", dims = 1:20)



Bcell_2 <- FindNeighbors(Bcell_2, reduction = "pca", dims = 1:20)

Bcell_2 <- FindClusters(Bcell_2, resolution = 0.26)
DimPlot(Bcell_2, reduction = 'tsne', label = T)
DimPlot(Bcell_2, reduction = 'tsne', label = T, group.by = 'old_identity')
DimPlot(Bcell_2, reduction = 'tsne', label = T, group.by = 'batch_2')
DimPlot(Bcell_2, reduction = 'tsne', label = T, group.by = 'Phase')


FeaturePlot(Bcell_2 ,reduction='tsne', cols = c('lightgrey', 'red'),
            features = c('CD27', 'SDC1', 'MKI67', 'TCL1A'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(Bcell_2 ,reduction='tsne', cols = c('lightgrey', 'red'),
            features = c('CD27'))

Idents(Bcell_2) <- Bcell_2$seurat_clusters
Bcell_2 <- RenameIdents(Bcell_2, '0' = 'MemoryB', '1' = 'MemoryB', 
                        '2' = 'MemoryB', '3' = 'NaiveB', '4' = 'Plasma', 
                        '5' = 'MemoryB', '6' = 'MemoryB')
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Bcell_2, reduction = 'tsne', label = T)

save(Bcell, Bcell_2, file = 'Bcell.RData')

