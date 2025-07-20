library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(harmony)
library(SCopeLoomR)



library(openxlsx)

load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')

library(randomForest) 


library(rlang)
library(ggforce)



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

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

color_used <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
                "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
                "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
                "#40E0D0","#5F9EA0","#FF1493",
                "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
                "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

Parekh_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Parekh_immunity/Parekh_step_2.rds')

Tom_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Tom_immunity/Tom_step_2.rds')

DN_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus_multiomics/DN_2.rds')

Parekh_step_2$source <- 'Parekh'
Tom_step_2$source <- 'Tom'
DN_2$source <- 'multiomics'


Parekh_step_2$idents <- as.character(Idents(Parekh_step_2))
Tom_step_2$idents <- as.character(Idents(Tom_step_2))
DN_2$idents <- as.character(Idents(DN_2))




Hs_earlyT <- merge(Parekh_step_2[, !Idents(Parekh_step_2) %in% c('other')], 
                   list(Tom_step_2[, !Idents(Tom_step_2) %in% c('nonT', 'CD8T')], DN_2))

Hs_earlyT <- NormalizeData(Hs_earlyT)
Hs_earlyT <- FindVariableFeatures(Hs_earlyT, selection.method = 'vst', nfeatures = 2000)
Hs_earlyT <- CellCycleScoring(Hs_earlyT, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)


plan('multisession', workers = 16)
options(future.globals.maxSize = 20 * 1024^3)
Hs_earlyT <- ScaleData(Hs_earlyT, vars.to.regress = c('S.Score', 'G2M.Score'))
plan('sequential')

Hs_earlyT <- RunPCA(Hs_earlyT, npcs = 50)

Hs_earlyT <- RunHarmony(Hs_earlyT, c('batch', 'source')) #))#, 'experiment', 

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT, ndims = 50, reduction = 'harmony')

Hs_earlyT <- RunUMAP(Hs_earlyT, reduction = "harmony", dims = 1:10, verbose = F)


options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT) <- Hs_earlyT$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT, reduction = 'umap', label=T)
DimPlot(Hs_earlyT, reduction = 'umap', label=T, group.by = 'Phase')

options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT) <- Hs_earlyT$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT, reduction = 'umap', label=T)
DimPlot(Hs_earlyT, reduction = 'umap', label=T, group.by = 'source')

plot_df <- cbind(Embeddings(Hs_earlyT, 'umap'), FetchData(Hs_earlyT, vars = c('idents', 'source', 'batch')))

options(repr.plot.width=16, repr.plot.height=16)
ggplot(plot_df, aes(UMAP_1, UMAP_2, color = idents))+
geom_point()+
scale_color_manual(values = color_used)+
facet_wrap(~source, ncol  = 2)

cluster_colors <- 
c(
    'CD8+Naive'='#53A85F','ETP_1' = '#297A6B','ETP_2'='#C05EA5','ETP'='#a7aad5',
    'ETP_3'='#e6ba77',
                             'CD4+Naive'='#D6E7A3','ETP-TP'='#a7aad5',
 'TP1'='#C05EA5','TP2'='#e6ba77','DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA',
    'DN4'='#ff927c','γδT_P'='#ee587d','pre_Treg1'='#E5D2DD',
                    'CD4+Memory_like'='#808040','CD8+Memory_like'='#F3B1A0',
'IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3',
    'Agonist_2'='#a83117',
    'Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','CD8+pre_T'='#625D9E', 
                    'Agonist_1' = '#509cbe',
                    'HSP_SP' = '#d98ed6', 
                    'pre_Treg2' = '#9566ac','Treg_diff_2'='#9566ac', 

    'CD8ααT' = '#e1812b',
   'MemoryB' = '#e0b139', 
    'NaiveB' = '#698c68', 
    'Plasma'='#ff7765',
   'MemoryNK'='#7ed670','NKT_like_CD8ααT'='#786085',
    'MatureNK'= '#b06fa2',
    'ImmatureNK'='#5065ad',
    'Lti'='#F0E442',
    'ZNF683+CD8ααT'='#de243a',
    'GNG4+CD8ααT'='#926c9c',
    'γδT_1'='#009E73','γδT_2'='#c39c62','γδT_3'='#9ba207',
    'IFN_CD8ααT'='#10c7ec',
    'Mast'='#e28f89', 
    "DC2"='#6e9a8d',#"#a8aad6",
    "aDC"="#73BDC6",
    "pDC"='#08843e',
    "Macrophage"="#E28F89",
    "DC1P"="#a3b9de","DC1"='#6aa6ec',#"#61a3f1",
    "Monocyte"="#8f6478","Neutrophil"="#51a0c6",#"#C19FC4"

    'Fb_1' ='#E5D2DD','cTEC_hi' ='#70a887', 'MKI67+Fb'='#E59CC4', 
                   'Edo'='#F3B1A0', 'Fb_2'='#D6E7A3', 'TEC(neuron)'='#57C3F3', 'mTEC_lo'='#5a829d',#'#476D87',
'VSMCs'='#E95C59', 'mTEC_hi'='#F1BB72','Ionocyte'='#cdb979',
                    'Lymph'='#AB3282', 'TEC(myo)'='#7db0bb',#'#91D0BE', 
                    'post_AIRE_mTEC'='#92a5a5', 'Tuft'='#c370ac', #'#3c97c4',
                    'Ionocyte'='#cdb979', 'Tuft/Ionocyte' = '#ff001c',
                    'Ciliated'='#cc5f91', 'cTEC_lo'='#c4c04f', 'MKI67+cTEC'='#ff7e00', 'Immature_TEC'='#9dabd5', 
                    'Mesothelial'='#d73e4b', 'Myelin'='#ff8ae0',
    "Endo"="#ec6960","Epi"="#ae80dc",
                              'Fibro'="#72bd90","Hema"="#9cbfe7",
    'thymus_CD8+Naive'='#53A85F',
                              
                             'thymus_CD4+Naive'='#D6E7A3','PB_CD8+RTE'='#276d6d',
'PB_CD4+TN'='#7bacd9','PB_CD8+TN'='#DDC4DA',
'PB_CD4+RTE'='#66c5b9'
   )

Hs_earlyT <- FindNeighbors(Hs_earlyT, reduction = "harmony", dims = 1:10, verbose = F)


plan('multisession', workers = 16)
options(future.globals.maxSize = 20 * 1024^3)
Hs_earlyT <- FindClusters(Hs_earlyT, resolution = seq(0.1,1, 0.05))
plan('sequential')

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(Hs_earlyT, reduction = 'umap',label = T, group.by = 'Phase')

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


Hs_earlyT




options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('PTPRC', 'CD7', 'CD1A','CD34'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c( 'CD7', 'CD1A','CD2','CD44'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c( 'MPO'),
            cols = c('lightgrey', 'red'))


unique(Hs_earlyT$source)



Hs_earlyT_multiomics <- Hs_earlyT[, Hs_earlyT$source == 'multiomics']

DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
Hs_earlyT_multiomics <- NormalizeData(Hs_earlyT_multiomics, normalization.method = 'CLR')


options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics, min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD2-TotalSeqB'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics, min.cutoff = 'q90',
            #max.cutoff = 'q95',
            features = c('CD44-TotalSeqB'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics, min.cutoff = 'q5', max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))



Hs_earlyT_S0 <- Hs_earlyT[, Hs_earlyT$RNA_snn_res.0.2 == '1']

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S0, ndims = 50, reduction = 'harmony')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S0) <- 'RNA'
FeaturePlot(Hs_earlyT_S0, features = c('PTCRA', 'BCL11B', 'RAG1', 'CD34', 'CD44', 'CD1A', 'CD2', 'CD7', 'DTX1'), 
            cols = c('lightgrey', 'red'))&
theme(plot.title=element_text(size = 20, face ='italic'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S0), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S0) <- 'RNA'
FeaturePlot(Hs_earlyT_S0,
            features = c('CD7', 'CD1A', 'CD2', 'CD44'))


DefaultAssay(Hs_earlyT_S0) <- 'RNA'
Hs_earlyT_S0 <- FindNeighbors(Hs_earlyT_S0, dims= 1:10, reduction = 'harmony',verbose = F)


options(future.globals.maxSize = 20 * 1024^3)
plan('multisession', workers = 16)
Hs_earlyT_S0 <- FindClusters(Hs_earlyT_S0, resolution = seq(0.05, 1, 0.05), verbose = F)
plan('sequential')

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S0, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


Hs_earlyT_S0_S1 <- Hs_earlyT_S0[, Hs_earlyT_S0$RNA_snn_res.0.4 %in% c('1','7', '3','4')]

DefaultAssay(Hs_earlyT_S0_S1) <- 'RNA'
Hs_earlyT_S0_S1 <- FindNeighbors(Hs_earlyT_S0_S1, dims= 1:10, reduction = 'harmony',verbose = F)


options(future.globals.maxSize = 20 * 1024^3)
plan('multisession', workers = 16)
Hs_earlyT_S0_S1 <- FindClusters(Hs_earlyT_S0_S1, resolution = seq(0.05, 1, 0.05), verbose = F)
plan('sequential')

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S0_S1, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 




Hs_earlyT_S0_S1_S1 <- Hs_earlyT_S0_S1[, Hs_earlyT_S0_S1$RNA_snn_res.0.7 %in% c('7')]


DefaultAssay(Hs_earlyT_S0_S1_S1) <- 'RNA'
Hs_earlyT_S0_S1_S1 <- FindNeighbors(Hs_earlyT_S0_S1_S1, dims= 1:10, reduction = 'harmony',verbose = F)


# options(future.globals.maxSize = 20 * 1024^3)
# plan('multisession', workers = 16)
Hs_earlyT_S0_S1_S1 <- FindClusters(Hs_earlyT_S0_S1_S1, resolution = seq(0.05, 1, 0.05), verbose = F)
# plan('sequential')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S0_S1_S1), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))


options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S0_S1_S1, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


Hs_earlyT_S0_S1_S1$Identity <- as.character(Hs_earlyT_S0_S1_S1$`RNA_snn_res.0.1`)
Hs_earlyT_S0_S1_S1$Identity[Hs_earlyT_S0_S1_S1$Identity %in% c('1')]<- 'Thy2a'
Hs_earlyT_S0_S1_S1$Identity[Hs_earlyT_S0_S1_S1$Identity %in% c('0')]<- 'Thy2b'

Idents(Hs_earlyT_S0_S1_S1) <- Hs_earlyT_S0_S1_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S0_S1_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))



Hs_earlyT_S0_S1_S2 <- Hs_earlyT_S0_S1[, Hs_earlyT_S0_S1$RNA_snn_res.0.7 %in% c('5')]


DefaultAssay(Hs_earlyT_S0_S1_S2) <- 'RNA'
Hs_earlyT_S0_S1_S2 <- FindNeighbors(Hs_earlyT_S0_S1_S2, dims= 1:10, reduction = 'harmony',verbose = F)


# options(future.globals.maxSize = 20 * 1024^3)
# plan('multisession', workers = 16)
Hs_earlyT_S0_S1_S2 <- FindClusters(Hs_earlyT_S0_S1_S2, resolution = seq(0.05, 1, 0.05), verbose = F)
# plan('sequential')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S0_S1_S2), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S0_S1_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S0_S1_S2, #min.cutoff = 'q5',
#            max.cutoff = 'q95',
            features = c('CD7', 'CD1A', 'CD2', 'CD44'))


options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S0_S1_S2, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 

Hs_earlyT_S0_S1_S2$Identity <- as.character(Hs_earlyT_S0_S1_S2$`RNA_snn_res.0.5`)
Hs_earlyT_S0_S1_S2$Identity[!Hs_earlyT_S0_S1_S2$Identity %in% c('2')]<- 'Thy2b'
Hs_earlyT_S0_S1_S2$Identity[Hs_earlyT_S0_S1_S2$Identity %in% c('2')]<- 'Thy2a'

Idents(Hs_earlyT_S0_S1_S2) <- Hs_earlyT_S0_S1_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S0_S1_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))

Hs_earlyT_S0_S1_S3 <- Hs_earlyT_S0_S1[, Hs_earlyT_S0_S1$RNA_snn_res.0.7 %in% c('0')]


DefaultAssay(Hs_earlyT_S0_S1_S3) <- 'RNA'
Hs_earlyT_S0_S1_S3 <- FindNeighbors(Hs_earlyT_S0_S1_S3, dims= 1:10, reduction = 'harmony',verbose = F)


# options(future.globals.maxSize = 20 * 1024^3)
# plan('multisession', workers = 16)
Hs_earlyT_S0_S1_S3 <- FindClusters(Hs_earlyT_S0_S1_S3, resolution = seq(0.05, 1, 0.05), verbose = F)
# plan('sequential')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S0_S1_S3), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S0_S1_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S0_S1_S3, #min.cutoff = 'q5',
#            max.cutoff = 'q95',
            features = c('CD7', 'CD1A', 'CD2', 'CD44'))


options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S0_S1_S3, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


Hs_earlyT_S0_S1_S3$Identity <- as.character(Hs_earlyT_S0_S1_S3$`RNA_snn_res.0.55`)
Hs_earlyT_S0_S1_S3$Identity[!Hs_earlyT_S0_S1_S3$Identity %in% c('0','3')]<- 'Thy2b'
Hs_earlyT_S0_S1_S3$Identity[Hs_earlyT_S0_S1_S3$Identity %in% c('0','3')]<- 'Thy2a'

Idents(Hs_earlyT_S0_S1_S3) <- Hs_earlyT_S0_S1_S3$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S0_S1_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))





Hs_earlyT_S0_S1$Identity <- as.character(Hs_earlyT_S0_S1$`RNA_snn_res.0.7`)
Hs_earlyT_S0_S1$Identity[Hs_earlyT_S0_S1$Identity %in% c('2','4','8')]<- 'Thy2b'
Hs_earlyT_S0_S1$Identity[Hs_earlyT_S0_S1$Identity %in% c('1','3','6')]<- 'Thy2a'

Idents(Hs_earlyT_S0_S1) <- Hs_earlyT_S0_S1$Identity
Hs_earlyT_S0_S1 <- Idents_proj(Hs_earlyT_S0_S1_S1, Hs_earlyT_S0_S1)
Hs_earlyT_S0_S1 <- Idents_proj(Hs_earlyT_S0_S1_S2, Hs_earlyT_S0_S1)
Hs_earlyT_S0_S1 <- Idents_proj(Hs_earlyT_S0_S1_S3, Hs_earlyT_S0_S1)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S0_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))

Hs_earlyT_S0$Identity <- as.character(Hs_earlyT_S0$`RNA_snn_res.0.4`)
Hs_earlyT_S0$Identity[Hs_earlyT_S0$Identity %in% c('6')]<- 'Thy1'
# Hs_earlyT_S0$Identity[Hs_earlyT_S0$Identity %in% c('1','7')]<- 'Thy2a'
# Hs_earlyT_S0$Identity[Hs_earlyT_S0$Identity %in% c('3','4')]<- 'Thy2b'
Hs_earlyT_S0$Identity[Hs_earlyT_S0$Identity %in% c('0','2', '5')]<- 'Thy2c'
Idents(Hs_earlyT_S0) <- Hs_earlyT_S0$Identity
Hs_earlyT_S0 <- Idents_proj(Hs_earlyT_S0_S1, Hs_earlyT_S0)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S0, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))






Hs_earlyT_S1 <- Hs_earlyT[, Hs_earlyT$RNA_snn_res.0.2 == '4']

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('PTCRA', 'BCL11B', 'RAG1', 'CD34', 'CD44', 'CD1A', 'CD2', 'CD7', 'DTX1'), 
            cols = c('lightgrey', 'red'))&
theme(plot.title=element_text(size = 20, face ='italic'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S1), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))


DefaultAssay(Hs_earlyT_S1) <- 'RNA'
Hs_earlyT_S1 <- FindNeighbors(Hs_earlyT_S1, dims= 1:10, reduction = 'harmony',verbose = F)


Hs_earlyT_S1 <- FindClusters(Hs_earlyT_S1, resolution = seq(0.05, 1, 0.05), verbose = F)

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S1, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


Hs_earlyT_S1$Identity <- as.character(Hs_earlyT_S1$`RNA_snn_res.0.3`)
Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('3')]<- 'Thy2a'
Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('2')]<- 'Thy2c'
Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('0','1')]<- 'Thy3'
Idents(Hs_earlyT_S1) <- Hs_earlyT_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))



Hs_earlyT_S2 <- Hs_earlyT[, Hs_earlyT$`RNA_snn_res.0.2`%in% c('2', '3')]

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S2, features = c('PTCRA', 'BCL11B', 'RAG1', 'CD34', 'CD44', 'CD1A', 'CD2', 'CD7', 'DTX1'), 
            cols = c('lightgrey', 'red'))&
theme(plot.title=element_text(size = 20, face ='italic'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
FeaturePlot(Hs_earlyT_multiomics[, intersect(colnames(Hs_earlyT_S2), 
                                             colnames(Hs_earlyT_multiomics))], min.cutoff = 'q5',
            max.cutoff = 'q95',
            features = c('CD7-TotalSeqB', 'CD1a-TotalSeqB', 'CD2-TotalSeqB', 'CD44-TotalSeqB'))


DefaultAssay(Hs_earlyT_S2) <- 'RNA'
Hs_earlyT_S2 <- FindNeighbors(Hs_earlyT_S2, dims= 1:10, reduction = 'harmony',verbose = F)


plan('multisession', workers = 16)
options(future.globals.maxSize = 20 * 1024^3)
Hs_earlyT_S2 <- FindClusters(Hs_earlyT_S2, resolution = seq(0.05, 1, 0.05), verbose = F)
plan('sequential')



options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S2, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 


Hs_earlyT_S2$Identity <- as.character(Hs_earlyT_S2$`RNA_snn_res.0.75`)
Hs_earlyT_S2$Identity[!Hs_earlyT_S2$Identity %in% c('6','8')]<- 'Thy3'
Hs_earlyT_S2$Identity[Hs_earlyT_S2$Identity %in% c('6','8')]<- 'Thy2c'
Idents(Hs_earlyT_S2) <- Hs_earlyT_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))



Hs_earlyT$Identity <- as.character(Hs_earlyT$`RNA_snn_res.0.2`)
Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('0','2','3')]<- 'Thy3'
# Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('1','4')]<- 'C2'
# Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('2')]<- 'C5'
Idents(Hs_earlyT) <- Hs_earlyT$Identity
Hs_earlyT <- Idents_proj(Hs_earlyT_S0, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S1, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S2, Hs_earlyT)

levels(Hs_earlyT) <- c('Thy1', 'Thy2a', 'Thy2b', 'Thy2c', 'Thy3')

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))

saveRDS(Idents(Hs_earlyT), file = 'Hs_earlyT_immunophenotypical_ID.rds')



Hs_earlyT_multiomics <- AddMetaData(Hs_earlyT_multiomics, metadata = Idents(Hs_earlyT), col.name = 'idents')

tmp_cluster_colors <- c('Thy1' = '#DC143C', 
                        'Thy2b' = '#808000', 
                        'Thy2a' = '#20B2AA',
                        'Thy2c' = '#FFA500',
                       'Thy3' = '#800080')

saveRDS(Hs_earlyT, file = 'Hs_earlyT_query.rds')

save(Hs_earlyT_S0, Hs_earlyT_S0_S1, Hs_earlyT_S0_S1_S1, Hs_earlyT_S0_S1_S2, Hs_earlyT_S0_S1_S3,
     Hs_earlyT_S1, Hs_earlyT_S2, file = 'Hs_earlyT_subsets.RData')

saveRDS(Hs_earlyT_multiomics, file = 'Hs_earlyT_multiomics.rds')



options(repr.plot.width=3.7, repr.plot.height=3.7)
p <- DimPlot(Hs_earlyT, reduction = 'umap',label = F, label.size = 13.75/.pt)+
scale_color_manual(values = tmp_cluster_colors)+
d_theme_w(size = 13.75)+
theme(plot.margin = ggplot2::margin(rep(0)),
      legend.text = element_text(size = 13.75, margin = ggplot2::margin(l=-0.3,unit = 'cm')),
      legend.position = c(0.01,0.19))+
guides(color = guide_legend(override.aes = list(size = 5)))+
scale_y_continuous(breaks = c(-4, 0, 4))
p

load('Thy1_3_df.RData')

options(repr.plot.width=3.7, repr.plot.height=3.7)
p + geom_mark_hull(linetype = 'dashed',color = 'black',
                    mapping = aes(x = UMAP_1, y = UMAP_2),data = Thy1_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                     radius = unit(1, "mm"),
                        con.type = 'straight',
                    con.border= 'none')+
 geom_mark_hull( linetype = 'dashed',color = 'black',
                    mapping = aes(x = UMAP_1, y = UMAP_2),data = Thy2a_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
 geom_mark_hull( linetype = 'dashed',color = 'black',
                    mapping = aes(x = UMAP_1, y = UMAP_2),data = Thy2b_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
 geom_mark_hull( linetype = 'dashed',color = 'black',
                    mapping = aes(x = UMAP_1, y = UMAP_2),data = Thy2c_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
 geom_mark_hull( linetype = 'dashed',color = 'black',
                    mapping = aes(x = UMAP_1, y = UMAP_2),data = Thy3_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')




Thy_df <- rbind(Thy1_df, Thy2a_df, Thy2b_df, Thy2c_df, Thy3_df)

DefaultAssay(Hs_earlyT_multiomics) <- 'ADT'
plot_df <- cbind(FetchData(Hs_earlyT_multiomics, rownames(Hs_earlyT_multiomics$ADT)), 
                Embeddings(Hs_earlyT_multiomics, reduction = 'umap'))

options(repr.plot.width=7.5, repr.plot.height=6)
tmp_SM <- 'CD7-TotalSeqB'
tmp_plot_df <- plot_df
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.95)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.05)
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 20)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = tmp_SM,hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))


options(repr.plot.width=7.5, repr.plot.height=6)
tmp_SM <- 'CD1a-TotalSeqB'
tmp_plot_df <- plot_df
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.95)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.05)
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 20)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = tmp_SM,hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))


options(repr.plot.width=7.5, repr.plot.height=6)
tmp_SM <- 'CD2-TotalSeqB'
tmp_plot_df <- plot_df
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.7)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.05)
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

p <- ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 20)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = tmp_SM,hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))
p$data <- p$data[order(p$data$`CD2-TotalSeqB`), ]
p


options(repr.plot.width=7.5, repr.plot.height=6)
tmp_SM <- 'CD44-TotalSeqB'
tmp_plot_df <- plot_df
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.99)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.9)
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 20)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = tmp_SM,hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))






plan('multicore', workers = 16)
options(future.globals.maxSize = 50 * 1024^3)
markers_Hs_earlyT <- FindAllMarkers(Hs_earlyT, pseudocount.use = 0, logfc.threshold = log(1.5))
plan('sequential')
saveRDS(markers_Hs_earlyT, file = 'markers_Hs_earlyT.rds')

openxlsx::write.xlsx(markers_Hs_earlyT, file = 'markers_Hs_earlyT.xlsx')



options(repr.plot.width=6, repr.plot.height=4)
DotPlot(Hs_earlyT, features = rev(c('CD34', 'BAALC', 'MEIS1', 'MEF2C', 'IRF8', 'BLNK', 'TYROBP',
                                'DTX1','LYL1', 'MPO','IL1B', 'IGFBP7',
                                'CD7', 'CD1A', 'CD44', 'CD2', 
                               'PTCRA','BCL11B','TCF7', 'RAG1')), 
       cols=  c('lightgrey', 'red'))+
scale_y_discrete(limits=rev)+
labs(x = NULL, y = NULL)+
d_theme_w(size  = 13.75)+
scale_size_continuous(range = c(1, 6))+
theme(axis.text.x = element_text(angle = 90,vjust =0.5, 
                                 hjust=1, face = 'italic'))+
guides(color = guide_colorbar(title = 'Expression'), size = guide_legend(title = '% Cells'))








