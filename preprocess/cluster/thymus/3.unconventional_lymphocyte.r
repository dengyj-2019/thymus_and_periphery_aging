library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
thymus_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/thymus_step_1_processed.rds')

library(openxlsx)
library(harmony)

unconventional <- 
    thymus_step_1[, Idents(thymus_step_1) %in% c('21', '28', '33')]
    
unconventional <- NormalizeData(unconventional)
unconventional <- FindVariableFeatures(unconventional, 
                                         selection.method = 'vst', nfeatures = 2000)
unconventional <- ScaleData(unconventional)
unconventional <- RunPCA(unconventional, npcs = 50, verbose = FALSE)
ElbowPlot(unconventional, ndims = 50)



unconventional <- RunUMAP(unconventional, reduction = "pca", dims = 1:30)
unconventional <- RunTSNE(unconventional, reduction = "pca", dims = 1:30)
unconventional <- FindNeighbors(unconventional, reduction = "pca", dims = 1:30)



unconventional <- FindClusters(unconventional, resolution = 0.4)
options(repr.plot.width=12, repr.plot.height=12)
DimPlot(unconventional, reduction = 'tsne', label = T)
DimPlot(unconventional, reduction = 'tsne', label = T, group.by = 'old_identity')
DimPlot(unconventional, reduction = 'tsne', label = T, group.by = 'batch')
FeaturePlot(unconventional, cols = c('lightgrey', 'red'), reduction = 'tsne',
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
VlnPlot(unconventional,  pt.size=0,
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

FeaturePlot(unconventional, reduction = 'tsne', cols = c('lightgrey', 'red'),
            features = c('CD40LG', 'NKG7', 'CD8B', 'CD8A', 'ZNF683', 'PTCRA', 'RAG1', 'IL7R', 'ANXA1', 'ELOVL4', 
                        'CD3E', 'CD3D', 'CD3G', 'CD19' ,"MS4A1"))

FeaturePlot(unconventional, reduction = 'tsne', cols = c('lightgrey', 'red'),
            features = c('PTCRA', 'RAG1','ELOVL4'
                        ))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(unconventional, reduction = 'tsne', cols = c('lightgrey', 'red'),
            features = c('TYROBP','NCR1', 'NCAM1', 'TRDC'
                        ))

###removed doublets and DP/SP contamination
unconventional_2 <- unconventional[, !Idents(unconventional) %in% 
                                         c('4', '7', '8', '11')]
unconventional_2 <- NormalizeData(unconventional_2)
unconventional_2 <- FindVariableFeatures(unconventional_2, 
                                           selection.method = 'vst', nfeatures = 2000)
unconventional_2 <- ScaleData(unconventional_2)
unconventional_2 <- RunPCA(unconventional_2, npcs = 50, verbose = FALSE)

ElbowPlot(unconventional_2, ndims = 50)

unconventional_2 <- RunUMAP(unconventional_2, reduction = "pca", dims = 1:30)
unconventional_2 <- RunTSNE(unconventional_2, reduction = "pca", dims = 1:30)
unconventional_2 <- FindNeighbors(unconventional_2, reduction = "pca", dims = 1:30)






unconventional_2 <- FindClusters(unconventional_2, resolution = 0.8)
DimPlot(unconventional_2, reduction = 'tsne', label = T)
DimPlot(unconventional_2, reduction = 'tsne', label = T, group.by = 'old_identity')
DimPlot(unconventional_2, reduction = 'tsne', group.by = 'batch_2')
DimPlot(unconventional_2, reduction = 'tsne', group.by = 'age')

FeaturePlot(unconventional_2, cols = c('lightgrey', 'red'), reduction = 'tsne',
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
VlnPlot(unconventional_2,  pt.size=0,
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))





###remove low quality cells
unconventional_3 <- unconventional_2[, !Idents(unconventional_2) %in% c('8')]
unconventional_3 <- NormalizeData(unconventional_3, verbose = T)
unconventional_3 <- FindVariableFeatures(unconventional_3, selection.method = "vst", nfeatures = 2000)
unconventional_3 <- ScaleData(unconventional_3, verbose = FALSE, vars.to.regress = 'batch')
unconventional_3 <- RunPCA(unconventional_3, npcs = 50, verbose = FALSE)

ElbowPlot(unconventional_3, ndims = 50)

unconventional_3 <- RunUMAP(unconventional_3, reduction = "pca", dims = 1:30)
unconventional_3 <- RunTSNE(unconventional_3, reduction = "pca", dims = 1:30)
unconventional_3 <- FindNeighbors(unconventional_3, reduction = "pca", dims = 1:30)



unconventional_3 <- FindClusters(unconventional_3, resolution = 0.3)
DimPlot(unconventional_3, reduction = 'umap', label = T)
DimPlot(unconventional_3, reduction = 'umap', label = T, group.by = 'old_identity')
DimPlot(unconventional_3, reduction = 'umap', group.by = 'batch_2')
DimPlot(unconventional_3, reduction = 'umap', group.by = 'Phase')
DimPlot(unconventional_3, reduction = 'umap', group.by = 'age')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(unconventional_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('CD8A', 'CD8B','CD3D', 'TCF7', 
                         'NCR1', 'CD4', 'CD40LG', 'CPA3'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(unconventional_3, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('NCAM1', 'FCGR3A', 'TYROBP', 'NCR1'))

sub_1 <- unconventional_3[, Idents(unconventional_3) %in% c('0')]
DefaultAssay(sub_1) <- 'RNA'
sub_1 <- NormalizeData(sub_1)
sub_1 <- FindVariableFeatures(sub_1, selection.method = 'vst', nfeatures = 2000)
sub_1 <- ScaleData(sub_1, vars.to.regress = 'nCount_RNA')
sub_1 <- RunPCA(sub_1, npcs = 50, verbose = FALSE)


ElbowPlot(sub_1, ndims = 50)

sub_1 <- RunUMAP(sub_1, reduction = "pca", dims = 1:30)
sub_1 <- RunTSNE(sub_1, reduction = "pca", dims = 1:30)
sub_1 <- FindNeighbors(sub_1, reduction = "pca", dims = 1:30)

options(repr.plot.width=8, repr.plot.height=6)
sub_1 <- FindClusters(sub_1, resolution = 1)
DimPlot(sub_1, reduction = 'umap', label = T)
DimPlot(sub_1, reduction = 'umap', group.by = 'old_identity', label=T)
DimPlot(sub_1, reduction = 'umap', group.by = 'batch_2')
DimPlot(sub_1, reduction = 'umap', group.by = 'age')
options(repr.plot.width=12, repr.plot.height=12)

FeaturePlot(sub_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('TRDC', 'CD8A', 'TRDV1', 'CD8B'))

FeaturePlot(sub_1, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('ZNF683', 'GNG4', 'ISG15'))

Idents(sub_1) <- sub_1$seurat_clusters
sub_1 <- RenameIdents(sub_1, '0' = 'γδT_3', '2' =  'γδT_2', '1' = 'γδT_1', '3' = 'GNG4+CD8ααT', 
                     '4' = 'ZNF683+CD8ααT', '5' = 'IFN_CD8ααT')
DimPlot(sub_1, reduction = 'umap', label = T)

sub_2 <- unconventional_3[, Idents(unconventional_3) %in% c('2')]
DefaultAssay(sub_2) <- 'RNA'
sub_2 <- NormalizeData(sub_2)
sub_2 <- FindVariableFeatures(sub_2, selection.method = 'vst', nfeatures = 2000)
sub_2 <- ScaleData(sub_2, vars.to.regress = 'nCount_RNA')
sub_2 <- RunPCA(sub_2, npcs = 50, verbose = FALSE)
ElbowPlot(sub_2, ndims = 50)


sub_2 <- RunUMAP(sub_2, reduction = "pca", dims = 1:20)
sub_2 <- RunTSNE(sub_2, reduction = "pca", dims = 1:20)
sub_2 <- FindNeighbors(sub_2, reduction = "pca", dims = 1:20)



sub_2 <- FindClusters(sub_2, resolution = 1.2)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(sub_2, reduction = 'umap', label = T)
DimPlot(sub_2, reduction = 'umap', group.by = 'old_identity', label=T)
DimPlot(sub_2, reduction = 'umap', group.by = 'batch_2')
DimPlot(sub_2, reduction = 'umap', group.by = 'age')
options(repr.plot.width=12, repr.plot.height=12)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(sub_2, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('NCAM1', 'FCGR3A', 'SELL', 'SPINK2', 'S1PR1', 'CXCR6','CCL5', 'KIT', 'TNFSF11'))

Idents(sub_2) <- sub_2$seurat_clusters
sub_2 <- RenameIdents(sub_2, '0' = 'MemoryNK', '3' = 'Lti', '1' = 'ImmatureNK', '2' = 'Lti')
DimPlot(sub_2, reduction = 'umap', label = T)

Idents(unconventional_3) <- unconventional_3$seurat_clusters
DimPlot(unconventional_3, reduction = 'umap', label = T)
DimPlot(unconventional_3, reduction = 'umap', label = T, group.by = 'old_identity')

Idents(unconventional_3) <- unconventional_3$seurat_clusters
unconventional_3 <- RenameIdents(unconventional_3,  '3' = 'NKT_like_CD8ααT',
                                 '4' = 'Mast', '1' = 'MatureNK')
combined_sub <- merge(sub_1, sub_2)
new_identity <- rbind(data.frame(b = as.character(Idents(combined_sub)), row.names = names(Idents(combined_sub))),
                  data.frame(b = as.character(Idents(unconventional_3)), row.names = names(Idents(unconventional_3)))[setdiff(colnames(unconventional_3), colnames(combined_sub)),,drop=F])
unconventional_3 <- AddMetaData(unconventional_3, metadata = new_identity, col.name ='new_identity') 
Idents(unconventional_3) <- unconventional_3$new_identity
DimPlot(unconventional_3, reduction = 'umap', label = T)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(unconventional_3, reduction = 'umap', label = T) + 
scale_color_manual(values = my_colors)
options(repr.plot.width=12, repr.plot.height=12)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(unconventional_3, reduction = 'umap',  group.by ='age_2') + 
scale_color_manual(values = c('<12'="#27468D", '12-24'="#268CC1", '25-39'="#50BF56",
                             '40-59'="#F6CB1E", '>=60'="#FD9303"))
DimPlot(unconventional_3, reduction = 'umap',  group.by ='Phase')
options(repr.plot.width=12, repr.plot.height=12)

save(unconventional, unconventional_2, unconventional_3, sub_1, sub_2, 
     file = 'unconventional.RData')
