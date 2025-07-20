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
            features = c('MS4A1', 'PAX5', 'IGKC', 'IGHM'))###

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(PB_step_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('JCHAIN', 'MZB1', 'DERL3', 'TNFRSF17'))###

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(PB_step_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3D', 'CD3E', 'CD3G'))###

Bcell <- PB_step_1[, Idents(PB_step_1) %in% c('4', '8')]

Bcell <- NormalizeData(Bcell)
Bcell <- FindVariableFeatures(Bcell, selection.method = "vst", nfeatures = 2000)
Bcell <- ScaleData(Bcell)
Bcell <- RunPCA(Bcell, npcs = 50)


Bcell <- RunHarmony(Bcell, 'batch')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(Bcell, ndims = 50, reduction = 'harmony')


Bcell <- RunUMAP(Bcell, reduction = "harmony", dims = 1:20)
Bcell <- FindNeighbors(Bcell, reduction = "harmony", dims = 1:20)


Bcell <- FindClusters(Bcell, resolution = 0.2)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Bcell, reduction = 'umap',label=T)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Bcell, reduction = 'umap',group.by = 'batch')
DimPlot(Bcell, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)#############checking doublets
FeaturePlot(Bcell, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3E', 'CD3G', 'CD3D', 'CD19', 'MS4A1',  'IGKC', 'CD79A', 'CD79B',
                        'PF4', 'CLU', 'GP9','TUBB1', 'CCL5', 'PF4V1', 'ITGB3','ITGA2B'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)#checking doublets
FeaturePlot(Bcell, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CEBPA', 'SPI1', 'CD33','HBA2'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)############
FeaturePlot(Bcell, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TBX21', 'ITGAX', 'MS4A1'), pt.size = 0.1)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Bcell, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Bcell, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MS4A1', 'PAX5', 'IGKC', 'IGHM', 'JCHAIN', 'MZB1', 'DERL3', 'TNFRSF17'))###

###removed doublets and cells of low quality
Bcell_2 <- Bcell[, !Idents(Bcell) %in% c('3','5', '6')]
Bcell_2 <- NormalizeData(Bcell_2)
Bcell_2 <- FindVariableFeatures(Bcell_2, selection.method = "vst", nfeatures = 2000)
Bcell_2 <- ScaleData(Bcell_2)#, vars.to.regress = c('nFeature_RNA'))
Bcell_2 <- RunPCA(Bcell_2, npcs = 50)

Bcell_2 <- RunHarmony(Bcell_2, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Bcell_2, ndims = 50, reduction = 'harmony')

Bcell_2 <- RunUMAP(Bcell_2, reduction = "harmony", dims = 1:20)
Bcell_2 <- FindNeighbors(Bcell_2, reduction = "harmony", dims = 1:20)


Bcell_2 <- FindClusters(Bcell_2, resolution = 0.6)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Bcell_2, reduction = 'umap',label =T)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Bcell_2, reduction = 'umap',group.by = 'batch')
DimPlot(Bcell_2, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Bcell_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)#############checking doublets
FeaturePlot(Bcell_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3E', 'CD3G', 'CD3D', 'CD19', 'MS4A1',  'IGKC', 'CD79A', 'CD79B',
                        'PF4', 'CLU', 'GP9','TUBB1', 'CCL5', 'PF4V1', 'ITGB3','ITGA2B'), pt.size = 0.1)

Bcell_3 <- Bcell_2[, !Idents(Bcell_2) %in% c('10')]
Bcell_3 <- NormalizeData(Bcell_3)
Bcell_3 <- FindVariableFeatures(Bcell_3, selection.method = "vst", nfeatures = 2000)
Bcell_3 <- ScaleData(Bcell_3)#, vars.to.regress = c('nFeature_RNA'))
Bcell_3 <- RunPCA(Bcell_3, npcs = 50)


Bcell_3 <- RunHarmony(Bcell_3, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Bcell_3, ndims = 50, reduction = 'harmony')

Bcell_3 <- RunUMAP(Bcell_3, reduction = "harmony", dims = 1:10)
Bcell_3 <- FindNeighbors(Bcell_3, reduction = "harmony", dims = 1:10)

Bcell_3 <- FindClusters(Bcell_3, resolution = 0.17)


options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Bcell_3, reduction = 'umap',label =T)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Bcell_3, reduction = 'umap',group.by = 'batch')
DimPlot(Bcell_3, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)#############checking doublets
FeaturePlot(Bcell_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3E', 'CD3G', 'CD3D', 'CD19', 'MS4A1',  'IGKC', 'CD79A', 'CD79B',
                        'PF4', 'CLU', 'GP9','TUBB1', 'CCL5', 'PF4V1', 'ITGB3','ITGA2B'), pt.size = 0.1)



options(future.globals.maxSize = 10 * 1024^3)
DefaultAssay(Bcell_3) <- "RNA"
plan('multisession', workers = 10)
markers_Bcell_3 <- FindAllMarkers(Bcell_3, pseudocount.use = 0)
write.xlsx(markers_Bcell_3, file = 'markers_Bcell_3.xlsx', row.names = T)
plan('sequential')

for(i in unique(markers_Bcell_3$cluster)){
    df <- markers_Bcell_3[markers_Bcell_3$cluster == i, ]
    df <- df[df$avg_logFC>0 & df$p_val_adj< 0.05, ]
    features <- df$gene[1:16]
    features <- features[!is.na(features)]
    p <- FeaturePlot(Bcell_3, features = features, cols = c('lightgrey','red'))
    print(p)
    print(i)
}

FeaturePlot(Bcell_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('IGKC', 'IGHA1', 'IGHA2', 'JCHAIN', 'IGHG4', 'MZB1', 'SDC1', 
                        'TNFRSF17', 'IGLL5'))####plasma

options(repr.plot.width=7.5, repr.plot.height=7.5)
#############MKI67 were expressed by plasma
FeaturePlot(Bcell_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MKI67'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Bcell_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD27', 'IGHG1','CD24', 'BACH2', 'TCL1A', 'IGHD', 'IGHM', 'CD38', 
                         'LAIR1','TBX21', 'ITGAX', 'CD19', 'CD72', 'FCRL5', 'FCRL2', 
                         'FCRL3'))####Naive and memoryB and ABC(CD19hi and CD11c+)

options(repr.plot.width=8, repr.plot.height=8)

FeaturePlot(Bcell_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('DERL3'))

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Bcell_3, label = T)

Idents(Bcell_3) <- Bcell_3$seurat_clusters
Bcell_3 <- RenameIdents(Bcell_3, 
                        '0' = 'Naive_B', '1' = 'Memory_B',
                        '3' = 'ABC', '2' = 'Plasma')
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(Bcell_3, label = T)

save(Bcell, Bcell_2,  Bcell_3,file = 'PB_Bcell.RData')
