library(Seurat)
library(dplyr)
library(future)
library(ggplot2)

library(data.table)

library(harmony)
library(openxlsx)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')

library(grid)
library(rlang)

Idents_proj <- function(
    obj1, 
    obj2
){
    new_identity <- rbind(data.frame(b = as.character(Idents(obj1)), row.names = names(Idents(obj1))),
                    data.frame(b = as.character(Idents(obj2)), 
                    row.names = names(Idents(obj2)))[setdiff(colnames(obj2), colnames(obj1)),,drop=F])
    obj2 <- AddMetaData(obj2, metadata = new_identity, col.name ='new_identity')            
    Idents(obj2) <- obj2$new_identity
    return(obj2)
}

PB_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/PB_step_1_processed.rds')

Idents(PB_step_1) <- PB_step_1$RNA_snn_res.0.1

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(PB_step_1, label = T)

T_NK <- PB_step_1[, Idents(PB_step_1) %in% c('0','1')]

T_NK <- NormalizeData(T_NK)
T_NK <- FindVariableFeatures(T_NK, selection.method = "vst", nfeatures = 2000)
T_NK <- ScaleData(T_NK)#, vars.to.regress = c('S.Score', 'G2M.Score'))
T_NK <- RunPCA(T_NK, npcs = 50)


T_NK <- RunHarmony(T_NK, 'batch')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(T_NK, ndims = 50, reduction = 'harmony')

T_NK <- RunUMAP(T_NK, reduction = "harmony", dims = 1:20)
T_NK <- FindNeighbors(T_NK, reduction = "harmony", dims = 1:20)

T_NK <- FindClusters(T_NK, resolution = seq(0.1, 1, 0.05))


options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(T_NK, reduction = 'umap',group.by = 'batch')
DimPlot(T_NK, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA')) 

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(T_NK, reduction = 'umap',group.by = 'orig_idents')+
scale_color_manual(values=cluster_colors)

options(repr.plot.width=16, repr.plot.height=16)###checking doublets
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('HBB', 'S100A8', 'S100A9', 'CD68', 'FCER1G', 'HLA-DRA', 
                        'LYZ', 'SPI1', 'S100A12','PF4', 'GP9', 'TUBB1', 'ITGA2B', 'NCR1', 'CD3D', 'TYROBP'))

options(repr.plot.width=16, repr.plot.height=16)###checking doublets
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MS4A1', 'CD19', 'PAX5', 'IGHM'))

Idents(T_NK) <- T_NK$RNA_snn_res.0.1
T_NK_2 <- T_NK[, !Idents(T_NK) %in% c('4', '6')]#####removed doublets
T_NK_2 <- NormalizeData(T_NK_2)
T_NK_2 <- FindVariableFeatures(T_NK_2, selection.method = "vst", nfeatures = 2000)

T_NK_2 <- ScaleData(T_NK_2)

T_NK_2 <- RunPCA(T_NK_2, npcs = 50)

T_NK_2 <- RunHarmony(T_NK_2, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2, ndims = 50, reduction = 'harmony')

T_NK_2 <- RunUMAP(T_NK_2, reduction = "harmony", dims = 1:20)
T_NK_2 <- FindNeighbors(T_NK_2, reduction = "harmony", dims = 1:20)

T_NK_2 <- FindClusters(T_NK_2, resolution = seq(0.1, 2, 0.05))

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(T_NK_2, reduction = 'umap',group.by = 'batch')
DimPlot(T_NK_2, reduction = 'umap',group.by = 'Phase')
  


options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 2, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK_2, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
}    

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(T_NK_2, reduction = 'umap',group.by = 'orig_idents')+
scale_color_manual(values = cluster_colors)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('CD4', 'CD40LG', 'CD8A', 'CD8B'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('CD3D', 'CD3E', 'CD3G', 'CD247'), cols = c('lightgrey', 'red'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(T_NK_2, features = c('CD38'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('IFIT1', 'OAS1', 'IFI44', 'USP18'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TYROBP', 'CD3E', 'TRDV2', 'TRGV9'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('IL2RA'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('NKG7', 'CST7', 'PLAC8', 'GZMK'))

options(repr.plot.width=16, repr.plot.height=16)##############在CD4 Naive中，也有一个
FeaturePlot(T_NK_2, features = c('CCR7', 'IL7R', 'TCF7', 'LEF1', 'SELL', 
                             'S1PR1', 'ACTN1', 'CD55', 'KLF2'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('IL2RB', 'FAS', 'CRTAM', 'ANXA1', 'ANXA2', 'PRDM1', 
                            'CD58', 'AHNAK', 'S100A4'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('PASK', 'ITGB1', 'GPR183', 'PTGER2'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('GZMH', 'CD27', 'GZMB', 'GZMK'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, features = c('HAVCR2', 'LAG3', 'TIGIT', 'IL10'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)##############Treg
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('FOXP3', 'IL2RA', 'TNFRSF9'))   

options(repr.plot.width=16, repr.plot.height=16)####gdT
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRGV9', 'TRDV2', 'TRDV1', 'TRDC')) 

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('KLRB1')) 

options(repr.plot.width=16, repr.plot.height=16)###NK
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('NCR1', 'NCAM1', 'TYROBP', 'FCGR3A')) 

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('MKI67')) 

T_NK_2_c1 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('1')]
T_NK_2_c1 <- NormalizeData(T_NK_2_c1)
T_NK_2_c1 <- FindVariableFeatures(T_NK_2_c1, selection.method = "vst", nfeatures = 2000)
T_NK_2_c1 <- ScaleData(T_NK_2_c1)
T_NK_2_c1 <- RunPCA(T_NK_2_c1, npcs = 50)

T_NK_2_c1 <- RunHarmony(T_NK_2_c1, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c1, ndims = 50, reduction = 'harmony')

options(repr.plot.width=16, repr.plot.height=16)
for(i in 1:10){
    p <- FeaturePlot(T_NK_2_c1, reduction = 'harmony', cols = c('lightgrey', 'red'),
                 dims = c(2*i-1, 2*i), features = c('CD4', 'CD8A', 'CD8B'))
    print(p)
}

###PC1 represents batch effect
T_NK_2_c1 <- RunUMAP(T_NK_2_c1, reduction = "harmony", dims = 2:20)
T_NK_2_c1 <- FindNeighbors(T_NK_2_c1, reduction = "harmony", dims = 2:20)

T_NK_2_c1 <- FindClusters(T_NK_2_c1, resolution = seq(0.1, 1.5, 0.05))

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1.5, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK_2_c1, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 

options(repr.plot.width=7.5, repr.plot.height=7.5)

DimPlot(T_NK_2_c1, reduction = 'umap',group.by = 'batch')
DimPlot(T_NK_2_c1, reduction = 'umap',group.by = 'Phase')

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(T_NK_2_c1, reduction = 'umap',group.by = 'RNA_snn_res.0.9', label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD4', 'CD40LG','CD8A', 'CD8B')) 


T_NK_2_c1$tmp_idents <- as.character(T_NK_2_c1$`RNA_snn_res.0.9`)
T_NK_2_c1$tmp_idents[T_NK_2_c1$RNA_snn_res.0.9 %in% c('0','1','3','4','11')] <- 'CD4+CTL'
T_NK_2_c1$tmp_idents[!T_NK_2_c1$RNA_snn_res.0.9 %in% c('0','1','3','4','11')] <- 'CD8+TEM-B'
Idents(T_NK_2_c1) <- T_NK_2_c1$tmp_idents
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2_c1, reduction = 'umap',label = T)



T_NK_2_c2 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('2')]


T_NK_2_c2 <- NormalizeData(T_NK_2_c2)
T_NK_2_c2 <- FindVariableFeatures(T_NK_2_c2, selection.method = "vst", nfeatures = 2000)
T_NK_2_c2 <- ScaleData(T_NK_2_c2)
T_NK_2_c2 <- RunPCA(T_NK_2_c2, npcs = 50)

T_NK_2_c2 <- RunHarmony(T_NK_2_c2, c('batch'))

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c2, ndims = 50, reduction = 'harmony')

options(repr.plot.width=16, repr.plot.height=16)
for(i in 1:10){
    p <- FeaturePlot(T_NK_2_c2, reduction = 'harmony', cols = c('lightgrey', 'red'),
                 dims = c(2*i-1, 2*i), features = c('CD4', 'CD8A', 'CD8B'))
    print(p)
}

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c2, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}

###PC1 represents batch effect
T_NK_2_c2 <- RunUMAP(T_NK_2_c2, reduction = "harmony", dims = c(2:10))
T_NK_2_c2 <- FindNeighbors(T_NK_2_c2, reduction = "harmony", dims = c(2:10))

options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c2 <- FindClusters(T_NK_2_c2, resolution = seq(0.1, 1.5, 0.05))
DimPlot(T_NK_2_c2, reduction = 'umap', label =T)

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1.5, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK_2_c2, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 



options(repr.plot.width=9, repr.plot.height=7.5)

DimPlot(T_NK_2_c2, reduction = 'umap',group.by = 'batch', label =T)
DimPlot(T_NK_2_c2, reduction = 'umap',group.by = 'Phase', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)

DimPlot(T_NK_2_c2, reduction = 'umap',group.by = 'orig_idents', label =T)+
scale_color_manual(values = cluster_colors)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRDC', 'TRDV2', 'TRGV9', 'TRDV1')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD3E', 'CD3G', 'CD3D', 'TYROBP')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD8A', 'CD8B', 'CXCR6', 'GZMK')) 

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(T_NK_2_c2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD27', 'LAG3')) 

Idents(T_NK_2_c2) <- T_NK_2_c2$`RNA_snn_res.0.1` 
options(repr.plot.width=8, repr.plot.height=8)
T_NK_2_c2 <- RenameIdents(T_NK_2_c2, '1'= 'CD27lo_Vδ1T', '0' = 'LAG3+NK')
DimPlot(T_NK_2_c2, reduction = 'umap',label = T)



T_NK_2_c3 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('3')]
T_NK_2_c3 <- NormalizeData(T_NK_2_c3)
T_NK_2_c3 <- FindVariableFeatures(T_NK_2_c3, selection.method = "vst", nfeatures = 2000)
T_NK_2_c3 <- ScaleData(T_NK_2_c3)
T_NK_2_c3 <- RunPCA(T_NK_2_c3, npcs = 50)
T_NK_2_c3 <- RunHarmony(T_NK_2_c3, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c3, ndims = 50, reduction = 'harmony')

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c3, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}

###PC1 represents batch effect
T_NK_2_c3 <- RunUMAP(T_NK_2_c3, reduction = "harmony", dims = 2:10)
T_NK_2_c3 <- FindNeighbors(T_NK_2_c3, reduction = "harmony", dims = 2:10)

options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c3 <- FindClusters(T_NK_2_c3, resolution = 0.3)
DimPlot(T_NK_2_c3, reduction = 'umap', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)

DimPlot(T_NK_2_c3, reduction = 'umap',group.by = 'batch', label =T)
DimPlot(T_NK_2_c3, reduction = 'umap',group.by = 'Phase', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)

DimPlot(T_NK_2_c3, reduction = 'umap',group.by = 'orig_idents', label =T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRDC', 'TRDV2', 'TRGV9', 'TRDV1')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('IL7R', 'CD27','CCR7', 'TCF7')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CCL5', 'GZMK', 'PRF1', 'GNLY')) 

Idents(T_NK_2_c3) <- T_NK_2_c3$seurat_clusters
options(repr.plot.width=8, repr.plot.height=8)
T_NK_2_c3 <- RenameIdents(T_NK_2_c3, '3'= 'CD27+Vδ1T', '2' = 'CD8+TCM', '0' = 'CD8+TEM-K', '1'='CD8+TEM-K')
DimPlot(T_NK_2_c3, reduction = 'umap',label = T)

T_NK_2_c6 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('6')]
T_NK_2_c6 <- NormalizeData(T_NK_2_c6)
T_NK_2_c6 <- FindVariableFeatures(T_NK_2_c6, selection.method = "vst", nfeatures = 2000)
T_NK_2_c6 <- ScaleData(T_NK_2_c6)
T_NK_2_c6 <- RunPCA(T_NK_2_c6, npcs = 50)
T_NK_2_c6 <- RunHarmony(T_NK_2_c6, 'batch')


options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c6, ndims = 50, reduction = 'harmony')

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c6, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}


###PC1 represents batch effect
T_NK_2_c6 <- RunUMAP(T_NK_2_c6, reduction = "harmony", dims = 2:10)
T_NK_2_c6 <- FindNeighbors(T_NK_2_c6, reduction = "harmony", dims = 2:10)

options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c6 <- FindClusters(T_NK_2_c6, resolution = 0.2)
DimPlot(T_NK_2_c6, reduction = 'umap', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)
DimPlot(T_NK_2_c6, reduction = 'umap',group.by = 'batch', label =T)
DimPlot(T_NK_2_c6, reduction = 'umap',group.by = 'Phase', label =T)


options(repr.plot.width=9, repr.plot.height=7.5)
DimPlot(T_NK_2_c6, reduction = 'umap',group.by = 'orig_idents', label =T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c6, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('GZMK', 'RORC', 'GATA3', 'CXCR3')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c6, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TBX21', 'IFNG', 'CCR6', 'CCR4')) 

T_NK_2_c6$tmp_idents <- as.character(T_NK_2_c6$seurat_clusters)
T_NK_2_c6$tmp_idents[T_NK_2_c6$seurat_clusters %in% c('0')] <- 'CD4+TEM'
T_NK_2_c6$tmp_idents[!T_NK_2_c6$seurat_clusters %in% c('0')] <- 'CD4+TEFF'###expressed Th2/Th17 related transcriptional factor
Idents(T_NK_2_c6) <- T_NK_2_c6$tmp_idents
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2_c6, reduction = 'umap',label = T)



T_NK_2_c0 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('0')]
T_NK_2_c0 <- NormalizeData(T_NK_2_c0)
T_NK_2_c0 <- FindVariableFeatures(T_NK_2_c0, selection.method = "vst", nfeatures = 2000)
T_NK_2_c0 <- ScaleData(T_NK_2_c0)#, vars.to.regress = 'batch')
T_NK_2_c0 <- RunPCA(T_NK_2_c0, npcs = 50)
T_NK_2_c0 <- RunHarmony(T_NK_2_c0, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c0, ndims = 50, reduction = 'harmony')

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c0, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}

###PC1 represents batch effect
T_NK_2_c0 <- RunUMAP(T_NK_2_c0, reduction = "harmony", dims = 2:10)
T_NK_2_c0 <- FindNeighbors(T_NK_2_c0, reduction = "harmony", dims = 2:10)


options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c0 <- FindClusters(T_NK_2_c0, resolution = 0.6)
DimPlot(T_NK_2_c0, reduction = 'umap', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)
DimPlot(T_NK_2_c0, reduction = 'umap',group.by = 'batch', label =T)
DimPlot(T_NK_2_c0, reduction = 'umap',group.by = 'Phase', label =T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c0, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('SOX4', 'PECAM1', 'LRRN3', 'CD38')) 

T_NK_2_c0$tmp_idents <- as.character(T_NK_2_c0$seurat_clusters)
T_NK_2_c0$tmp_idents[T_NK_2_c0$seurat_clusters %in% c('4')] <- 'CD4+RTE'
T_NK_2_c0$tmp_idents[!T_NK_2_c0$seurat_clusters %in% c('4')] <- 'CD4+TN'
Idents(T_NK_2_c0) <- T_NK_2_c0$tmp_idents
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2_c0, reduction = 'umap',label = T)



T_NK_2_c4 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('4')]
T_NK_2_c4 <- NormalizeData(T_NK_2_c4)
T_NK_2_c4 <- FindVariableFeatures(T_NK_2_c4, selection.method = "vst", nfeatures = 2000)
T_NK_2_c4 <- ScaleData(T_NK_2_c4)#, vars.to.regress = 'batch')
T_NK_2_c4 <- RunPCA(T_NK_2_c4, npcs = 50)
T_NK_2_c4 <- RunHarmony(T_NK_2_c4, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c4, ndims = 50, reduction = 'harmony')

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c4, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}

###PC1 represents batch effect
T_NK_2_c4 <- RunUMAP(T_NK_2_c4, reduction = "harmony", dims = 2:15)
T_NK_2_c4 <- FindNeighbors(T_NK_2_c4, reduction = "harmony", dims = 2:15)

options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c4 <- FindClusters(T_NK_2_c4, resolution = 0.2)
DimPlot(T_NK_2_c4, reduction = 'umap', label =T)

options(repr.plot.width=9, repr.plot.height=7.5)
DimPlot(T_NK_2_c4, reduction = 'umap',group.by = 'batch', label =T)
DimPlot(T_NK_2_c4, reduction = 'umap',group.by = 'Phase', label =T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c4, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('SOX4', 'TOX', 'LRRN3', 'CD38')) 

T_NK_2_c4$tmp_idents <- as.character(T_NK_2_c4$seurat_clusters)
T_NK_2_c4$tmp_idents[T_NK_2_c4$seurat_clusters %in% c('1')] <- 'CD8+RTE'
T_NK_2_c4$tmp_idents[!T_NK_2_c4$seurat_clusters %in% c('1')] <- 'CD8+TN'
Idents(T_NK_2_c4) <- T_NK_2_c4$tmp_idents
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2_c4, reduction = 'umap',label = T)

T_NK_2_c8 <- T_NK_2[, T_NK_2$RNA_snn_res.0.65 %in% c('8')]
T_NK_2_c8 <- NormalizeData(T_NK_2_c8)
T_NK_2_c8 <- FindVariableFeatures(T_NK_2_c8, selection.method = "vst", nfeatures = 2000)
T_NK_2_c8 <- ScaleData(T_NK_2_c8)#, vars.to.regress = 'batch')
T_NK_2_c8 <- RunPCA(T_NK_2_c8, npcs = 50)
T_NK_2_c8 <- RunHarmony(T_NK_2_c8, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(T_NK_2_c8, ndims = 50, reduction = 'harmony')

options(repr.plot.width=8, repr.plot.height=8)
for(i in 1:10){
    p <- DimPlot(T_NK_2_c8, reduction = 'harmony', group.by='batch',
                 dims = c(2*i-1, 2*i))
    print(p)
}


###PC1 represents batch effect
T_NK_2_c8 <- RunUMAP(T_NK_2_c8, reduction = "harmony", dims = 2:10)
T_NK_2_c8 <- FindNeighbors(T_NK_2_c8, reduction = "harmony", dims = 2:10)

options(repr.plot.width=7.5, repr.plot.height=7.5)
T_NK_2_c8 <- FindClusters(T_NK_2_c8, resolution = 0.2)
DimPlot(T_NK_2_c8, reduction = 'umap', label =T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c8, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRDV2', 'TRGV9', 'NCR3', 'TRAV1-2')) 

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(T_NK_2_c8, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRDC')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2_c8, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('RORC', 'RORA', 'KLRB1')) 

T_NK_2_c8$tmp_idents <- as.character(T_NK_2_c8$seurat_clusters)
T_NK_2_c8$tmp_idents[T_NK_2_c8$seurat_clusters %in% c('0')] <- 'MAIT'
T_NK_2_c8$tmp_idents[!T_NK_2_c8$seurat_clusters %in% c('0')] <- 'Vδ2Vγ9T'
Idents(T_NK_2_c8) <- T_NK_2_c8$tmp_idents
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2_c8, reduction = 'umap',label = T)







options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('GATA3', 'CCR4', "CCR8", 'CXCR3', 'HAVCR2', 'PDCD1', 'ICOS', 'IFNG' )) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('GZMK', 'GZMB', 'CD27', 'KLRG1', 'GZMH', 'TBX21', 'FCGR3A', 'CX3CR1', 'EOMES'
                        )) 

options(repr.plot.width=9, repr.plot.height=9)
Idents(T_NK_2) <- T_NK_2$`RNA_snn_res.0.65`
DimPlot(T_NK_2, label = T)

options(repr.plot.width=10, repr.plot.height=8)
Idents(T_NK_2) <- T_NK_2$`RNA_snn_res.0.65`
T_NK_2 <- RenameIdents(T_NK_2, '9' = 'Treg', '13' = 'Tex', '5' = 'CD4+TCM', 
                       '11'='IFN_T','10' = 'T_un',
                        '7' = 'CD56_Dim_NK', '12' = 'CD56_Bright_NK')
T_NK_2 <- Idents_proj(T_NK_2_c0, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c1, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c2, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c3, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c4, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c6, T_NK_2)   
T_NK_2 <- Idents_proj(T_NK_2_c8, T_NK_2)   
DimPlot(T_NK_2, label = T)

T_NK_2$age_3 <- ''
T_NK_2$age_3[T_NK_2$age <= 12] <- '<=12'
T_NK_2$age_3[T_NK_2$age >= 13 & T_NK_2$age <= 39] <- '13-39'
T_NK_2$age_3[T_NK_2$age >= 40 & T_NK_2$age<99] <- '40-99'
T_NK_2$age_3[T_NK_2$age >= 100] <- '>=100'
T_NK_2$age_3 <- factor(T_NK_2$age_3, levels = c('<=12', '13-39', '40-99', '>=100'))

save(T_NK, T_NK_2, T_NK_2_c0, T_NK_2_c1,T_NK_2_c2,T_NK_2_c3,
     T_NK_2_c4,T_NK_2_c6,T_NK_2_c8, file = 'T_NK.RData')

saveRDS(T_NK_2, file = 'T_NK_2.rds')

T_NK_2_idents <- as.character(Idents(T_NK_2))
names(T_NK_2_idents) <- colnames(T_NK_2)
saveRDS(T_NK_2_idents, file = 'T_NK_2_idents.rds')

saveRDS(Idents(T_NK_2), file = '/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_tcr_idents.rds')






Tcell <- T_NK_2[, !grepl('NK', Idents(T_NK_2))]
Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell)#, vars.to.regress = c('percent.mt'))
Tcell <- RunPCA(Tcell, npcs = 50)

Tcell <- RunHarmony(Tcell, 'batch')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(Tcell, ndims = 50, reduction = 'harmony')

Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:20)

options(repr.plot.width=8, repr.plot.height=8)
p <- DimPlot(Tcell,  label = F,pt.size = 0.3)+
scale_color_manual(values = cluster_colors)
p

saveRDS(Tcell, file = 'Tcell.rds')
