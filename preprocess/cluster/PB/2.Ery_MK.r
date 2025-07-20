library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)

library(harmony)
library(openxlsx)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')

PB_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/PB_step_1_processed.rds')

Idents(PB_step_1) <- PB_step_1$RNA_snn_res.0.1

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(PB_step_1, label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(PB_step_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('PF4', 'PPBP', 'HBA1', 'GYPA'))####RBC+MK

Ery_MK <- PB_step_1[, Idents(PB_step_1) %in% c('5', '7')]



Ery_MK <- NormalizeData(Ery_MK)
Ery_MK <- FindVariableFeatures(Ery_MK, selection.method = "vst", nfeatures = 2000)
Ery_MK <- ScaleData(Ery_MK)
Ery_MK <- RunPCA(Ery_MK, npcs = 50)

Ery_MK <- RunHarmony(Ery_MK, 'batch')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(Ery_MK, ndims = 50, reduction = 'harmony')

Ery_MK <- RunUMAP(Ery_MK, reduction = "harmony", dims = 1:20)
Ery_MK <- FindNeighbors(Ery_MK, reduction = "harmony", dims = 1:20)

Ery_MK <- FindClusters(Ery_MK, resolution = 0.15)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Ery_MK, reduction = 'umap',label = T)
DimPlot(Ery_MK, reduction = 'umap',group.by = 'batch')
DimPlot(Ery_MK, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Ery_MK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)#############checking doublets
FeaturePlot(Ery_MK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MS4A1', 'SPI1', 'CD33','HBA2', 'PF4', "CD3E", 'TYROBP', 'NCR1', 'KLRD1'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)#checking doublets
FeaturePlot(Ery_MK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('PF4', 'HBA1', 'PTPRC', 'CD3E'))

options(repr.plot.width=16, repr.plot.height=16)#############checking doublets
FeaturePlot(Ery_MK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TUBB1', 'HBB', 'GP9','CLU', 'GATA2', "ITGA2B", 
                         'HBA2', 'HBA1'), pt.size = 0.1)

###extracted pure cells
Ery_MK_2 <- Ery_MK[, Idents(Ery_MK) %in% c('0', '1', '3')]
Ery_MK_2 <- NormalizeData(Ery_MK_2)
Ery_MK_2 <- FindVariableFeatures(Ery_MK_2, selection.method = "vst", nfeatures = 2000)
Ery_MK_2 <- ScaleData(Ery_MK_2)
Ery_MK_2 <- RunPCA(Ery_MK_2, npcs = 50)

Ery_MK_2 <- RunHarmony(Ery_MK_2, 'batch')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(Ery_MK_2, ndims = 50, reduction = 'harmony')

Ery_MK_2 <- RunUMAP(Ery_MK_2, reduction = "harmony", dims = 1:20)
Ery_MK_2 <- FindNeighbors(Ery_MK_2, reduction = "harmony", dims = 1:20)

options(repr.plot.width=7, repr.plot.height=7)
Ery_MK_2 <- FindClusters(Ery_MK_2, resolution = 0.05)
DimPlot(Ery_MK_2, reduction = 'umap',label = T)
DimPlot(Ery_MK_2, reduction = 'umap',group.by = 'batch')
DimPlot(Ery_MK_2, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Ery_MK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('PF4', 'HBB', 'TUBB1', 'HBA2'))

options(repr.plot.width=8, repr.plot.height=8)
Idents(Ery_MK_2) <- Ery_MK_2$seurat_clusters
Ery_MK_2 <- RenameIdents(Ery_MK_2, '0' = 'MK', '1' = 'RBC')
DimPlot(Ery_MK_2, reduction = 'umap',label = T)

save(Ery_MK, Ery_MK_2, file = 'Ery_MK.RData')
