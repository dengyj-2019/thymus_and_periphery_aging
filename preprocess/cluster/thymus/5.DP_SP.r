library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
library(openxlsx)
library(ggpubr)
library(harmony)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
thymus_step_1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/thymus_step_1_processed.rds')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/unconventional/unconventional.RData')

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

T_dev_identity <- c(0,1,2,3,4,5,6,7, 10,11,12,
                    13,14,15, 16, 17, 20,22, 
                    25,27, 29, 31,35,36)
T_dev <- thymus_step_1[, Idents(thymus_step_1) %in% T_dev_identity]
T_dev <- NormalizeData(T_dev, verbose = T)
T_dev <- FindVariableFeatures(T_dev, selection.method = "vst", nfeatures = 2000)
T_dev <- ScaleData(T_dev)
T_dev <- RunPCA(T_dev, npcs = 50, verbose = FALSE)
T_dev <- RunHarmony(T_dev, 'batch')

ElbowPlot(T_dev, ndims = 50, reduction = 'harmony')

T_dev <- RunUMAP(T_dev, reduction = "harmony", dims = 1:20)
T_dev <- FindNeighbors(T_dev, reduction = "harmony", dims = 1:20)



options(repr.plot.width=12, repr.plot.height=12)
T_dev <- FindClusters(T_dev, resolution = 0.5)


DimPlot(T_dev, reduction = 'umap', label = T)
DimPlot(T_dev, reduction = 'umap', group.by = 'Phase')
DimPlot(T_dev, reduction = 'umap', group.by = 'batch_2')
DimPlot(T_dev, reduction = 'umap', group.by = 'batch')
DimPlot(T_dev, reduction = 'umap', group.by = 'old_identity')
FeaturePlot(T_dev, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('nCount_RNA', 'percent.mt'))
VlnPlot(T_dev, pt.size = 0, features = c('nCount_RNA', 'percent.mt'))

DimPlot(unconventional, reduction = 'tsne', label = T)

###removed cells of low quality and combined mixed DP/SP cells from clustering analysis of unconventional lymphocyte 
T_dev_2 <- merge(T_dev[, !Idents(T_dev) %in% c('4', '8', '9', '15')], 
                 unconventional[, Idents(unconventional) %in% c('4', '8')])
T_dev_2 <- NormalizeData(T_dev_2, verbose = T)
T_dev_2 <- FindVariableFeatures(T_dev_2, selection.method = "vst", nfeatures = 2000)
T_dev_2 <- ScaleData(T_dev_2)
T_dev_2 <- RunPCA(T_dev_2, npcs = 50, verbose = FALSE)
T_dev_2 <- RunHarmony(T_dev_2, 'batch')

ElbowPlot(T_dev_2, ndims = 50, reduction = 'harmony')

T_dev_2 <- RunUMAP(T_dev_2, reduction = "harmony", dims = 1:20)
T_dev_2 <- FindNeighbors(T_dev_2, reduction = "harmony", dims = 1:20)


T_dev_2 <- FindClusters(T_dev_2, resolution = 0.1)
DimPlot(T_dev_2, reduction = 'umap', label = T)
DimPlot(T_dev_2, reduction = 'umap', group.by = 'Phase')
DimPlot(T_dev_2, reduction = 'umap', group.by = 'batch_2')
DimPlot(T_dev_2, reduction = 'umap', group.by = 'batch')
DimPlot(T_dev_2, reduction = 'umap', group.by = 'old_identity')


FeaturePlot(T_dev_2, cols = c('lightgrey', 'red'), 
           features = c('CD40LG', 'TOX2', 'CCR7', 'SELL', 'TCF7', 'CD4', 'CD8A', 'BACH2', 
                       'CD8B'))


FeaturePlot(T_dev_2, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('nCount_RNA', 'nFeature_RNA','percent.mt'))
VlnPlot(T_dev_2, pt.size = 0, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))#, pt.size=0)

FeaturePlot(T_dev_2, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('CD4', 'CD8B','CD8A', 'CD40LG', 'BACH2', 'TOX2', 'CCR7', 'SELL', 'TCF7'))

SP <- T_dev_2[, Idents(T_dev_2) %in% c('1', '3')]
DefaultAssay(SP) <- 'RNA'
SP <- NormalizeData(SP, verbose = T)
SP <- FindVariableFeatures(SP, selection.method = "vst", nfeatures = 2000)
SP <- ScaleData(SP)
SP <- RunPCA(SP, npcs = 50, verbose = FALSE)
SP <- RunHarmony(SP, 'batch')


ElbowPlot(SP, ndims = 50, reduction = 'harmony')

SP <- RunUMAP(SP, reduction = "harmony", dims = 1:20)
SP <- FindNeighbors(SP, reduction = "harmony", dims = 1:20)

SP <- FindClusters(SP, resolution = seq(0.1, 1, 0.025))

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.025)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(SP, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
}  

options(repr.plot.width=12, repr.plot.height=12)
SP <- FindClusters(SP, resolution = 0.45)#0.5
DimPlot(SP, reduction = 'umap', label = T)
DimPlot(SP, reduction = 'umap', group.by = 'Phase')
DimPlot(SP, reduction = 'umap', group.by = 'batch')
DimPlot(SP, reduction = 'umap', group.by = 'old_identity', label = T)

FeaturePlot(SP, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('ISG15', 'OAS1','IFIT1', 'MX2'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('CD27', 'PDCD1', 'GNG4', 'MIR155HG'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('CD4', 'CD40LG','CD8A', 'CD8B'))


options(repr.plot.width=16, repr.plot.height=16)##############
FeaturePlot(SP, features = c('IL2RB', 'FAS', 'CRTAM', 'ANXA1', 'ANXA2', 'PRDM1', 
                            'CD58', 'AHNAK', 'S100A4', 'ICAM2', 'PRF1', 
                                'GPR183', 'CCL5', 'GZMA', 'PTGER2'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)##########
FeaturePlot(SP, features = c('CCR7', 'IL7R', 'TCF7', 'LEF1', 'SELL', 
                             'S1PR1', 'ACTN1', 'CD55', 'KLF2'), cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(SP, features = c('EOMES', 'ICOS', 'TBX21', 
                            'CD244', 'CD44', 'ZEB2', 'IFNG', 
                            'CD27', 'FOXP1'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)##############在CD4 Naive中，也有一个
FeaturePlot(SP, features = c('CST7', 'NKG7', 'GZMK', 'GZMB'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)##############在CD4 Naive中，也有一个
FeaturePlot(SP, features = c('HAVCR2', 'LAG3', 'CTLA4', 'TIGIT'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)#############
FeaturePlot(SP, features = c('CD28', 'CD160', 'CXCR5', 'ALCAM', 'PDCD1', 'KLRK1', 'HIF1A', 'CD101', 
                            'BCL6', 'MYB', 'TRIM68', 'LY9', 'SLAMF6'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)##############在CD4 Naive中，也有一个
FeaturePlot(SP, features = c('IL2RA', 'CTLA4', 'FOXP3', 'CD40LG'), cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(SP, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('PASK', 'TNFRSF4'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('RAG1', 'AQP3', 'ELOVL4', 'RAG2'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap', cols  = c('lightgrey', 'red'),pt.size = 0.1,
            c('BCL11A', 'SOX4', 'CCR9', 'TOX2', 'ID3', 'CD38', 'SIRPG', 'CD28', 'BACH2', 
             "SATB1", 'STAT5A', 'TOX', 'TRAV8-2', 'CD2', 'SPINK2'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(SP, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c( 
                        "DHRS3", 'RAG1', 'DNTT', 'TNFRSF9', 'CCR7', 'CTLA4', 'GNG4', 'RORC'
                       ))
                        

pre_T <- SP[, Idents(SP) %in% c('0')]
pre_T <- NormalizeData(pre_T, verbose = T)
pre_T <- FindVariableFeatures(pre_T, selection.method = "vst", nfeatures = 2000)
pre_T <- ScaleData(pre_T)#, vars.to.regress = c('batch'))
pre_T <- RunPCA(pre_T, npcs = 50, verbose = FALSE)#, features = setdiff(pre_T$RNA@var.features, gogenes))
pre_T <- RunHarmony(pre_T, 'batch')


ElbowPlot(pre_T, ndims = 50, reduction = 'harmony')

pre_T <- RunUMAP(pre_T, reduction = "harmony", dims = 1:20)
pre_T <- FindNeighbors(pre_T, reduction = "harmony", dims = 1:20)


options(repr.plot.width=7.5, repr.plot.height=7.5)
pre_T <- FindClusters(pre_T, resolution = 0.3)
DimPlot(pre_T, reduction = 'umap', label = T)


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(pre_T, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c( 'CD40LG', 'CD4',  'CD8A','CD8B',  'CD27','PDCD1', 'TOX2', 'RAG1', 'CD1E'
                       ))

options(repr.plot.width=7, repr.plot.height=7)
Idents(pre_T) <- pre_T$seurat_clusters
pre_T <- RenameIdents(pre_T, '1' = 'αβ_Entry', '0' = 'CD4+pre_T',
                      '2' = 'CD8+pre_T', '3' = 'CD4+pre_T')
DimPlot(pre_T, reduction = 'umap', label = T)

Idents(SP) <- SP$seurat_clusters
CD4_Naive <- SP[, Idents(SP) %in% c('1')]
CD4_Naive <- NormalizeData(CD4_Naive, verbose = T)
CD4_Naive <- FindVariableFeatures(CD4_Naive, selection.method = "vst", nfeatures = 2000)
CD4_Naive <- ScaleData(CD4_Naive)#, vars.to.regress = c('batch'))
CD4_Naive <- RunPCA(CD4_Naive, npcs = 50, verbose = FALSE)#, features = setdiff(CD4_Naive$RNA@var.features, gogenes))
CD4_Naive <- RunHarmony(CD4_Naive, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD4_Naive, ndims = 50, reduction = 'harmony')

CD4_Naive <- RunUMAP(CD4_Naive, reduction = "harmony", dims = 1:10)
CD4_Naive <- FindNeighbors(CD4_Naive, reduction = "harmony", dims = 1:10)

CD4_Naive <- RunTSNE(CD4_Naive, reduction = "harmony", dims = 1:10)

mycolors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
'#968175')

options(repr.plot.width=7.5, repr.plot.height=7.5)
CD4_Naive <- FindClusters(CD4_Naive, resolution = 0.2)
DimPlot(CD4_Naive, reduction = 'umap', label = T)
DimPlot(CD4_Naive, reduction = 'umap', group.by = 'batch')+
scale_color_manual(values = mycolors)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(CD4_Naive, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('CD4', 'CD40LG', 'CD8A', 'CD8B', 'SOX4', 'CCR9', 'TOX', 'TOX2'))

sub_1 <- CD4_Naive[, Idents(CD4_Naive) %in% c('2')]
sub_1 <- NormalizeData(sub_1, verbose = T)
sub_1 <- FindVariableFeatures(sub_1, selection.method = "vst", nfeatures = 2000)
sub_1 <- ScaleData(sub_1, vars.to.regress = c('batch'))
sub_1 <- RunPCA(sub_1, npcs = 50, verbose = FALSE)#, features = setdiff(sub_1$RNA@var.features, gogenes))


options(repr.plot.width=5, repr.plot.height=5)
ElbowPlot(sub_1, ndims = 50, reduction = 'pca')

sub_1 <- RunUMAP(sub_1, reduction = "pca", dims = 1:10)
sub_1 <- FindNeighbors(sub_1, reduction = "pca", dims = 1:10)

options(repr.plot.width=6, repr.plot.height=6)
sub_1 <- FindClusters(sub_1, resolution = 1.3)
DimPlot(sub_1, reduction = 'umap', label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(sub_1, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c(  'CD4', 'CD40LG', 'CD8A', 'CD8B'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(sub_1, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c(  'SOX4', 'CCR9', 'TOX', 'TOX2'))

options(repr.plot.width=7, repr.plot.height=7)
Idents(sub_1) <- sub_1$seurat_clusters
sub_1 <- RenameIdents(sub_1, '0' = 'CD4+SOX4+Naive', '1' = 'CD4+SOX4+Naive','6' = 'CD4+SOX4+Naive',
                      '2' = 'CD8+SOX4+Naive','3' = 'CD8+SOX4+Naive','4' = 'CD8+SOX4+Naive','5' = 'CD8+SOX4+Naive'
                         )
DimPlot(sub_1, reduction = 'umap', label = T)

options(repr.plot.width=6, repr.plot.height=6)
Idents(CD4_Naive) <- CD4_Naive$seurat_clusters
CD4_Naive <- RenameIdents(CD4_Naive, '0' = 'CD4+SOX4+Naive', 
                         '1' = 'CD4+SOX4-Naive')
new_identity <- rbind(data.frame(b = as.character(Idents(sub_1)), row.names = names(Idents(sub_1))),
                  data.frame(b = as.character(Idents(CD4_Naive)), row.names = names(Idents(CD4_Naive)))[setdiff(colnames(CD4_Naive), colnames(sub_1)),,drop=F])
CD4_Naive <- AddMetaData(CD4_Naive, metadata = new_identity, col.name ='new_identity')            
Idents(CD4_Naive) <- CD4_Naive$new_identity
DimPlot(CD4_Naive, reduction = 'umap', label = T)



Idents(SP) <- SP$seurat_clusters
CD8_Naive <- SP[, Idents(SP) %in% c('6')]
CD8_Naive <- NormalizeData(CD8_Naive, verbose = T)
CD8_Naive <- FindVariableFeatures(CD8_Naive, selection.method = "vst", nfeatures = 2000)
CD8_Naive <- ScaleData(CD8_Naive)#, vars.to.regress = c('batch'))
CD8_Naive <- RunPCA(CD8_Naive, npcs = 50, verbose = FALSE)#, features = setdiff(CD8_Naive$RNA@var.features, gogenes))
CD8_Naive <- RunHarmony(CD8_Naive, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD8_Naive, ndims = 50, reduction = 'harmony')

CD8_Naive <- RunUMAP(CD8_Naive, reduction = "harmony", dims = 1:10)
CD8_Naive <- FindNeighbors(CD8_Naive, reduction = "harmony", dims = 1:10)

options(repr.plot.width=7.5, repr.plot.height=7.5)
CD8_Naive <- FindClusters(CD8_Naive, resolution = 0.1)
DimPlot(CD8_Naive, reduction = 'umap', label = T)
DimPlot(CD8_Naive, reduction = 'umap', group.by = 'batch')+
scale_color_manual(values = mycolors)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(CD8_Naive, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('CD4', 'CD40LG', 'CD8A', 'CD8B', 'SOX4', 'CCR9', 'TOX', 'TOX2'))

options(repr.plot.width=6, repr.plot.height=6)
Idents(CD8_Naive) <- CD8_Naive$seurat_clusters
CD8_Naive <- RenameIdents(CD8_Naive, '0' = 'CD8+SOX4+Naive', 
                         '1' = 'CD8+SOX4-Naive')
DimPlot(CD8_Naive, reduction = 'umap', label = T)





Treg_like <- SP[, Idents(SP) %in% c('7', '8')]
Treg_like <- NormalizeData(Treg_like, verbose = T)
Treg_like <- FindVariableFeatures(Treg_like, selection.method = "vst", nfeatures = 2000)
Treg_like <- ScaleData(Treg_like, vars.to.regress = c('batch'))
Treg_like <- RunPCA(Treg_like, npcs = 50, verbose = FALSE)#, features = setdiff(Treg_like$RNA@var.features, gogenes))
#Treg_like <- RunHarmony(Treg_like, 'batch')


options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(Treg_like, ndims = 50, reduction = 'pca')

Treg_like <- RunUMAP(Treg_like, reduction = "pca", dims = 1:20)#20
Treg_like <- FindNeighbors(Treg_like, reduction = "pca", dims = 1:20)

options(repr.plot.width=7.5, repr.plot.height=7.5)
Treg_like <- FindClusters(Treg_like, resolution = 0.5)#0.5
DimPlot(Treg_like, reduction = 'umap', label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Treg_like, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c( 'CD40LG', 'CD4',  'CD8A','CD8B', 'FOXP3', 'IL2RA', 'IL2RB', 'PDCD1', 'MIR155HG', 
                        "DHRS3", 'RAG1', 'DNTT', 'TNFRSF9', 'CCR7', 'CTLA4', 'GNG4'
                       ))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Treg_like, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c('RAG1', 'RAG2', 'DNTT', 'ELOVL4'
                       ))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Treg_like, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c( 'TCF7', 'SELL', 'LEF1', 'CCR7'
                       ))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Treg_like, reduction = 'umap',cols = c('lightgrey', 'red'), 
           features = c( 
                        'MIF', 'TNFRSF9', 'ITM2A','NME1', 'ATP5MC1', 'ERH', 
               'IL2RA','RANBP1', 'MIR155HG',  'CD40LG', 'NCL', 'BTG1', 'HCST', 'GIMAP4', 
               'SMC4', 'ARHGEF1'
                       ))

options(repr.plot.width=7, repr.plot.height=7)
Idents(Treg_like) <- Treg_like$seurat_clusters
Treg_like <- RenameIdents(Treg_like, '0' = 'pre_Treg1', '3' = 'pre_Treg1',
                      '2' = 'pre_Treg2', '1' = 'Treg', '4' ='doublets'
                         )
DimPlot(Treg_like, reduction = 'umap', label = T)

Idents(SP) <- SP$seurat_clusters
DimPlot(SP, reduction = 'umap', label = T)

options(repr.plot.width=9, repr.plot.height=9)
Idents(SP) <- SP$seurat_clusters
SP <- RenameIdents(SP, '3' = 'αβ_Entry', '9' = 'Agonist_2', '2' = 'Agonist_1','5'= 'CD8+Memory_like',
                     #'6' = 'CD8+Naive', 
                     '4'= 'CD4+Memory_like', '10' = 'IFN_response', '11' = 'doublets')
SP <- Idents_proj(CD4_Naive, SP)
SP <- Idents_proj(CD8_Naive, SP)
SP <- Idents_proj(Treg_like, SP)
SP <- Idents_proj(pre_T, SP)

DimPlot(SP, reduction = 'umap', label = T)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(SP, reduction = 'umap', group.by= 'batch')

DimPlot(T_dev_2, label = T)

Idents(T_dev_2) <- T_dev_2$seurat_clusters
T_dev_2 <- RenameIdents(T_dev_2, '4' = 'HSP_DP', '2' = 'DP_P',
                      '0' = 'DP_Q')
DimPlot(T_dev_2, reduction = 'umap', label = T)

T_dev_3 <- merge(T_dev_2[, Idents(T_dev_2) %in% c('HSP_DP','DP_P','DP_Q')], 

                     SP[, !Idents(SP) %in% c('doublets')] 
                )

T_dev_3 <- NormalizeData(T_dev_3, verbose = T)
T_dev_3 <- FindVariableFeatures(T_dev_3, selection.method = "vst", nfeatures = 2000)
T_dev_3 <- ScaleData(T_dev_3)
T_dev_3 <- RunPCA(T_dev_3, npcs = 50, verbose = FALSE)
T_dev_3 <- RunHarmony(T_dev_3, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(T_dev_3, ndims = 50, reduction = 'harmony')

T_dev_3 <- RunUMAP(T_dev_3, reduction = "harmony", dims = 1:20)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(T_dev_3, reduction = 'umap', label = T)



matureSP <- SP[, grepl('Naive',Idents(SP))]

options(repr.plot.width=5, repr.plot.height=5)
ElbowPlot(matureSP, ndims = 50, reduction = 'harmony')

matureSP <- RunUMAP(matureSP, reduction = 'harmony', dims = 1:20, verbose = F)

lbls_func <- function(x){
    x <- gsub('CD8\\+', 'CD8 ', x)
    x <- gsub('CD4\\+', 'CD4 ', x)
    x <- gsub('Naive', 'SP', x)
}

options(repr.plot.width=4, repr.plot.height=4)
p <- DimPlot(matureSP, reduction = 'umap', label = T)+
guides(color = F)+
scale_color_manual(values = cluster_colors, label = lbls_func)+
d_theme_w(size = 13.75)+
theme(axis.text = element_blank(), axis.ticks = element_blank())

p$layers[[2]]$data$ident <- lbls_func(p$layers[[2]]$data$ident)

p$layers[[2]]$geom_params$size <- 13.75/.pt

p

matureSP$gsea_idents <- as.character(Idents(matureSP))

matureSP$gsea_idents <- gsub('SOX4\\+', 'SOX4_pos_', matureSP$gsea_idents)
matureSP$gsea_idents <- gsub('SOX4\\-', 'SOX4_neg_', matureSP$gsea_idents)
matureSP$gsea_idents <- gsub('\\+', '_', matureSP$gsea_idents)


setwd('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/')
get_gsea_file2(matureSP$RNA@data, group_var = matureSP$gsea_idents, exp_name = 'matureSP_matrix.txt', 
              cls_name = 'matureSP_cls.cls',gseapy = F, add_prefix = F)
setwd('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev')

#Linux environment
#!/bin/bash
###
gsea-cli.sh GSEA -res matureSP_matrix.txt -cls matureSP_cls.cls\
#CD4_SOX4_pos_Naive_versus_CD4_SOX4_neg_Naive \
-gmx /data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/GO_KEGG_REACTOME_C7.gmt \
-collapse No_Collapse \
-rpt_label CD4_SOX4_pos_Naive_versus_CD4_SOX4_neg_Naive


#Linux environment
#!/bin/bash
###
gsea-cli.sh GSEA -res matureSP_matrix.txt -cls matureSP_cls.cls\
#CD8_SOX4_pos_Naive_versus_CD8_SOX4_neg_Naive \
-gmx /data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/GO_KEGG_REACTOME_C7.gmt \
-collapse No_Collapse \
-rpt_label CD8_SOX4_pos_Naive_versus_CD8_SOX4_neg_Naive





T_dev_3$id_lv1 <- as.character(Idents(T_dev_3))
T_dev_3$id_lv2 <- T_dev_3$id_lv1
T_dev_3$id_lv2[grep('CD4.*Naive', T_dev_3$id_lv2)] <- 'CD4+Naive'
T_dev_3$id_lv2[grep('CD8.*Naive', T_dev_3$id_lv2)] <- 'CD8+Naive'
Idents(T_dev_3) <- T_dev_3$id_lv2

T_dev_3$age_3 <- ''
T_dev_3$age_3[T_dev_3$age <=12] <- '<=12'
T_dev_3$age_3[T_dev_3$age >= 13 & T_dev_3$age <= 39 ] <- '13-39'
T_dev_3$age_3[T_dev_3$age >=40] <- '>=40'



IFN_response <- T_dev_3[, Idents(T_dev_3) %in% c('IFN_response')]

IFN_response <- FindVariableFeatures(IFN_response, selection.method = "vst", nfeatures = 2000,verbose = F)
IFN_response <- ScaleData(IFN_response,verbose = F)#, vars.to.regress = c('batch'))
IFN_response <- RunPCA(IFN_response, npcs = 50, verbose = F)#, 
IFN_response <- RunHarmony(IFN_response, 'batch', verbose = F)

IFN_response <- RunUMAP(IFN_response, reduction = "harmony", dims = 1:20, verbose = F)
IFN_response <- FindNeighbors(IFN_response, reduction = "harmony", dims = 1:20, verbose = F)

IFN_response <- FindClusters(IFN_response, resolution = seq(0.1, 1, 0.1), verbose = F)

options(repr.plot.width=5, repr.plot.height=5)
DimPlot(IFN_response, reduction = 'umap',label = T, group.by = 'RNA_snn_res.0.6')

options(repr.plot.width=9, repr.plot.height=9)
FeaturePlot(IFN_response, features = c('CD8A', 'CD8B', 'CD4', 'CD40LG'))

options(repr.plot.width=5, repr.plot.height=5)
Idents(IFN_response) <- IFN_response$`RNA_snn_res.0.6`
IFN_response <- RenameIdents(IFN_response, '4' = 'CD8+IFN_response', 
                             '0' = 'CD4+IFN_response', 
                            '1' = 'CD4+IFN_response',
                            '2' = 'CD4+IFN_response', 
                            '3' = 'CD4+IFN_response')

DimPlot(IFN_response, reduction = 'umap', label = F)                     

saveRDS(IFN_response, file = 'IFN_response.rds')









saveRDS(T_dev_3, file = 'T_dev_3.rds')

save(T_dev, T_dev_2,T_dev_3,Treg_like, SP,
     pre_T, CD4_Naive, sub_1, file = 'T_dev.RData')

T_dev_3$tcr_idents <- as.character(Idents(T_dev_3))
T_dev_3$idents <- as.character(Idents(T_dev_3))
saveRDS(T_dev_3$tcr_idents, file = 'T_dev_tcr_idents.rds')
saveRDS(T_dev_3$idents, file = 'T_dev_idents.rds')
saveRDS(matureSP, file = 'matureSP.rds')
