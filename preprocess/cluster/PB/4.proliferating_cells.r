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
            features = c('CDK1', 'MKI67', 'CCNA2'))####proliferation

proliferating <- PB_step_1[, Idents(PB_step_1) %in% c('9')]



proliferating <- NormalizeData(proliferating)
proliferating <- FindVariableFeatures(proliferating, selection.method = "vst", nfeatures = 2000)
proliferating <- ScaleData(proliferating, vars.to.regress = c('S.Score', 'G2M.Score', 'batch'))
proliferating <- RunPCA(proliferating, npcs = 50)

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(proliferating, ndims = 50, reduction = 'pca')

proliferating <- RunUMAP(proliferating, reduction = "pca", dims = 1:10)
proliferating <- FindNeighbors(proliferating, reduction = "pca", dims = 1:10)

proliferating <- FindClusters(proliferating, resolution = 0.15)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(proliferating, reduction = 'umap',label = T)
DimPlot(proliferating, reduction = 'umap',group.by = 'batch')
DimPlot(proliferating, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(proliferating, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)######checking doublets
FeaturePlot(proliferating, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MS4A1', 'SPI1', 'CD33','HBA2', 'PF4', "CD3E", 'TYROBP', 'NCR1', 'KLRD1'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)#############
FeaturePlot(proliferating, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3D', 'CD247', 'CD3G'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(proliferating, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('PF4', 'HBA1', 'PTPRC', 'CD3E'))

#removed doublets
proliferating_2 <- proliferating[, !Idents(proliferating) %in% c('2')]
proliferating_2 <- NormalizeData(proliferating_2)
proliferating_2 <- FindVariableFeatures(proliferating_2, selection.method = "vst", nfeatures = 2000)
proliferating_2 <- ScaleData(proliferating_2, vars.to.regress = c('S.Score', 'G2M.Score', 'batch'))
proliferating_2 <- RunPCA(proliferating_2, npcs = 50)

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(proliferating_2, ndims = 50, reduction = 'pca')

proliferating_2 <- RunUMAP(proliferating_2, reduction = "pca", dims = 1:10)
proliferating_2 <- FindNeighbors(proliferating_2, reduction = "pca", dims = 1:10)

options(repr.plot.width=7, repr.plot.height=7)
proliferating_2 <- FindClusters(proliferating_2, resolution = 0.2)
DimPlot(proliferating_2, reduction = 'umap',label = T)
DimPlot(proliferating_2, reduction = 'umap',group.by = 'batch')
DimPlot(proliferating_2, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)#############
FeaturePlot(proliferating_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MS4A1', 'SPI1', 'CD33','HBA2', 'PF4', "CD3E", 'TYROBP', 'NCR1', 'KLRD1'), pt.size = 0.1)

options(repr.plot.width=8, repr.plot.height=8)
Idents(proliferating_2) <- proliferating_2$seurat_clusters
proliferating_2 <- RenameIdents(proliferating_2, '0' = 'proliferating_NK', '3' = 'proliferating_B', 
                               '1'='proliferating_T', '2'='proliferating_T')
DimPlot(proliferating_2, reduction = 'umap',label = T)

save(proliferating, proliferating_2, file = 'proliferating.RData')
