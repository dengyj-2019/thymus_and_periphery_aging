library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)

load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
thymus_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/thymus_step_1_processed.rds')

library(openxlsx)

library(harmony)

Myeloid <- thymus_step_1[, Idents(thymus_step_1) %in% c('24', '37', '38')]
Myeloid <- NormalizeData(Myeloid, verbose = T)
Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 2000)
Myeloid <- ScaleData(Myeloid)#, vars.to.regress = c("batch_2"))
Myeloid <- RunPCA(Myeloid, npcs = 50, verbose = FALSE)
ElbowPlot(Myeloid, ndims = 50)

Myeloid <- RunTSNE(Myeloid, reduction = "pca", dims = 1:30)
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:30)



Myeloid <- FindNeighbors(Myeloid, reduction = "pca", dims = 1:30)

Myeloid <- FindClusters(Myeloid, resolution = 3)
DimPlot(Myeloid, reduction = 'tsne', label = T)
DimPlot(Myeloid, reduction = 'tsne', label =T, group.by = 'old_identity')
DimPlot(Myeloid, reduction = 'tsne', label = T, group.by = 'batch_2')
FeaturePlot(Myeloid, reduction = 'tsne', 
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
DimPlot(Myeloid, reduction = 'tsne', label = T, group.by = 'Phase')
DimPlot(Myeloid, reduction = 'tsne', label = T, group.by = 'batch')
options(repr.plot.width=9, repr.plot.height=9)
FeaturePlot(Myeloid, reduction = 'tsne', 
            features = c('CD3D', 'CLEC9A', 'MKI67', 
                         'S100A8', 'CLEC10A', 'CLEC4A', 
                         'CD19', 'CD4', 'CD3E'), pt.size = 0.5)
FeaturePlot(Myeloid, reduction = 'tsne', 
            features = c('CD8A', 'CD8B', 'CD3E', 'CD3G'), pt.size = 0.5)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Myeloid, reduction = 'tsne', 
            features = c('SPI1', 'IL3RA', 'TYROBP', 'NCR1'), pt.size = 0.5)

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(Myeloid, reduction = 'tsne', 
            features = c('CPA3'), pt.size = 0.5)

###removed T-cell doublets
Myeloid_2 <- Myeloid[, !Idents(Myeloid) %in% c('0', '9','13', '16', '20')]

DefaultAssay(Myeloid_2) <- 'RNA'
Myeloid_2 <- NormalizeData(Myeloid_2, verbose = T)
Myeloid_2 <- FindVariableFeatures(Myeloid_2, selection.method = "vst", nfeatures = 2000)
Myeloid_2 <- ScaleData(Myeloid_2,verbose = T, vars.to.regress = c('nCount_RNA'))
Myeloid_2 <- RunPCA(Myeloid_2, npcs = 50, verbose = FALSE)


ElbowPlot(Myeloid_2, ndims = 50)
Myeloid_2 <- RunHarmony(Myeloid_2, "batch", plot_convergence = TRUE)

Myeloid_2 <- RunTSNE(Myeloid_2, reduction = "harmony", dims = 1:30)
Myeloid_2 <- RunUMAP(Myeloid_2, reduction = "harmony", dims = 1:30)
Myeloid_2 <- FindNeighbors(Myeloid_2, reduction = "harmony", dims = 1:30)


DefaultAssay(Myeloid_2) <- 'RNA'
Myeloid_2 <- FindClusters(Myeloid_2, resolution = 0.2)
DimPlot(Myeloid_2, reduction = 'tsne', label = T)
DimPlot(Myeloid_2, reduction = 'tsne', label=T, group.by = 'old_identity')
DimPlot(Myeloid_2, reduction = 'tsne',  group.by = 'Phase')
DimPlot(Myeloid_2, reduction = 'tsne',  group.by = 'batch_2')

FeaturePlot(Myeloid_2, reduction = 'tsne', cols = c('lightgrey', 'red'), 
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

options(repr.plot.width=12, repr.plot.height=15)
FeaturePlot(Myeloid_2, reduction = 'tsne', cols = c('lightgrey', 'red'), 
            features = c('CLEC9A','CPVL','LYZ',
'C1orf54','DNASE1L3','LGALS2','IRF8','RGCC','BATF3','HLA-DPA1',
'CST3','S100A10','C1orf162','S100B','CPNE3','HLA-DPB1','CADM1','LMNA','SNX3','HLA-DQB1'))
options(repr.plot.width=12, repr.plot.height=12)


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Myeloid_2, reduction = 'tsne', cols = c('lightgrey', 'red'), 
            features = c('CLEC10A','CD1C','LAMP3','CCL19','GZMB','LILRA4','C1QC',
                         'C1QB','CDC20','MKI67','CLEC9A','XCR1','VCAN','FCN1','ELANE','AZU1'))




DimPlot(Myeloid_2, reduction = 'tsne', label = T)
DimPlot(Myeloid_2, reduction = 'tsne', label=T, group.by = 'old_identity')

Myeloid_2 <- RenameIdents(Myeloid_2, '0' = 'DC2', '2' = 'aDC', '1' = 'pDC', 
                         '4' = 'Macrophage', '5' = 'DC1P', '3' = 'DC1', 
                          '6' = 'Monocyte', '7' = 'Neutrophil')

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Myeloid_2, reduction = 'tsne',label =T)

save(Myeloid, Myeloid_2, file = 'Myeloid.RData')
