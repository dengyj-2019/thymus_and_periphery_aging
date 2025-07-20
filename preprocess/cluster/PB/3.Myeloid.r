library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
library(harmony)
library(openxlsx)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
library(pheatmap)

PB_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/PB_step_1_processed.rds')

Idents(PB_step_1) <- PB_step_1$RNA_snn_res.0.1

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(PB_step_1, label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(PB_step_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3D', 'CD3E', 'CD3G'))###

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(PB_step_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('SPI1', 'CD33', 'TYROBP'))###

Myeloid <- PB_step_1[, Idents(PB_step_1) %in% c('2', '6')]
Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 2000)
Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid, npcs = 50)


Myeloid <- RunHarmony(Myeloid, 'batch')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(Myeloid, ndims = 50, reduction = 'harmony')

Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:20)
Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:20)

Myeloid <- FindClusters(Myeloid, resolution = 0.85)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Myeloid, reduction = 'umap',label=T)

options(repr.plot.width=16, repr.plot.height=16)
DimPlot(Myeloid, reduction = 'umap',group.by = 'batch')
DimPlot(Myeloid, reduction = 'umap',group.by = 'Phase')
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))   


options(repr.plot.width=16, repr.plot.height=16)###checking doublets
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3E', 'CD3D', 'CD3G','PF4', 'TUBB1', 'PPBP',  'MS4A1', 'CD79A', 'CD79B','CD69', 
                        'NR4A1', 'ID3', 'EGR1', 'SPI1', 'FCGR3A', 'CD14'))   

options(repr.plot.width=16, repr.plot.height=16)###pDC
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('IL3RA', 'TCF4', 'LILRA4','CLEC4C'))

options(repr.plot.width=16, repr.plot.height=16)###HSPC
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD34', 'KIT', 'SPINK2', 'HOXA9'))

options(repr.plot.width=8, repr.plot.height=8)###cDC2
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CLEC10A'))

options(repr.plot.width=8, repr.plot.height=8)###cDC1
FeaturePlot(Myeloid, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CLEC9A'))

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Myeloid, reduction = 'umap',group.by = 'orig_idents')

###removed doublets and cells of low quality
Myeloid_2 <- Myeloid[, !Idents(Myeloid) %in% c('5', '12','13', '14','7', '16', '15')]
Myeloid_2 <- NormalizeData(Myeloid_2)
Myeloid_2 <- FindVariableFeatures(Myeloid_2, selection.method = "vst", nfeatures = 2000)
Myeloid_2 <- ScaleData(Myeloid_2)
Myeloid_2 <- RunPCA(Myeloid_2, npcs = 50)                  

Myeloid_2 <- RunHarmony(Myeloid_2, 'batch') 

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(Myeloid_2, ndims = 50, reduction = 'harmony')

Myeloid_2 <- RunUMAP(Myeloid_2, reduction = "harmony", dims = 1:20)
Myeloid_2 <- FindNeighbors(Myeloid_2, reduction = "harmony", dims = 1:20)

Myeloid_2 <- FindClusters(Myeloid_2,resolution =  0.4)


options(repr.plot.width=9, repr.plot.height=9)
DimPlot(Myeloid_2, reduction = 'umap',label = T)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(Myeloid_2, reduction = 'umap',group.by = 'batch')
DimPlot(Myeloid_2, reduction = 'umap',group.by = 'Phase')


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Myeloid_2, reduction = 'umap',group.by = 'orig_idents')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CXCR2', 'CXCR1', 'FCGR3B')) 

options(repr.plot.width=16, repr.plot.height=16)###checking doublets
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3G','PF4', 'TUBB1', 'PPBP',  'MS4A1', 'CD79A', 'CD79B','CD69', 
                        'NR4A1', 'ID3', 'EGR1', 'SPI1', 'CLEC9A', 'CD14', 'CLEC10A'))   

options(repr.plot.width=16, repr.plot.height=16)###看看
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'CLEC9A', 'CLEC10A', 'CD14','FCGR3A'))   

options(repr.plot.width=16, repr.plot.height=16)###看看
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'CD34', 'KIT', 'SPINK2'))   

options(repr.plot.width=16, repr.plot.height=16)###pDC
FeaturePlot(Myeloid_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('IL3RA', 'TCF4', 'LILRA4','CLEC4C'))

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(Myeloid_2, reduction = 'umap',label = T)

Idents(Myeloid_2) <- Myeloid_2$seurat_clusters
Myeloid_2 <- RenameIdents(Myeloid_2, '0' = 'Mono1', '1' = 'Mono1', '2' = 'Mono1', '4' = 'Mono1', '5' = 'cDC2', 
                          '6' = 'pDC', '7'='HSPC', '8' = 'Mono3','9' ='cDC1','3' = 'Mono2'
                         )
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Myeloid_2, reduction = 'umap', label =T)

save(Myeloid_2, Myeloid, file = 'Myeloid.RData')
