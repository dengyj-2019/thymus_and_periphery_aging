library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(harmony)
library(pheatmap)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')


library(plyr)

Hs_earlyT <- readRDS('old/250404/Hs_earlyT.rds')

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


color_used <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
                "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
                "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
                "#40E0D0","#5F9EA0","#FF1493",
                "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
                "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

earlyT_and_ILC <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/NA_first_revised_additional_data/earlyT_and_ILC.rds')

Parekh_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Parekh_immunity/Parekh_step_2.rds')

Tom_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Tom_immunity/Tom_step_2.rds')

DN_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus_multiomics/DN_2.rds')

ETP_DN_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN_2.rds')

NA_first_revised_additional_ETP_DN <- earlyT_and_ILC[, Idents(earlyT_and_ILC) %in% c('ETP_DN')]



NA_first_revised_additional_ETP_DN$batch_2 <- 'sixth'

DN_2$batch_2 <- 'fifth'

Parekh_step_2$batch_2 <- Parekh_step_2$batch

Tom_step_2$batch_2 <- Tom_step_2$batch

ETP_DN_2$batch_3 <- ETP_DN_2$batch_2

NA_first_revised_additional_ETP_DN$batch_3 <- NA_first_revised_additional_ETP_DN$batch
DN_2$batch_3 <- DN_2$batch

Parekh_step_2$batch_3 <- Parekh_step_2$batch

Tom_step_2$batch_3 <- Tom_step_2$batch

ETP_DN_2$source <- 'Hu'
DN_2$source <- 'Hu'
NA_first_revised_additional_ETP_DN$source <- 'Hu'

Parekh_step_2$source <- 'Parekh'

Tom_step_2$source <- 'Tom'

Parekh_step_2$idents <- as.character(Idents(Parekh_step_2))
Tom_step_2$idents <- as.character(Idents(Tom_step_2))
DN_2$idents <- as.character(Idents(DN_2))
ETP_DN_2$idents <- as.character(Idents(ETP_DN_2))
NA_first_revised_additional_ETP_DN$idents <- as.character(NA_first_revised_additional_ETP_DN$RNA_snn_res.0.9)

Hs_earlyT <- merge(Parekh_step_2[, !Idents(Parekh_step_2) %in% c('other')], 
                   list(Tom_step_2[, !Idents(Tom_step_2) %in% c('nonT', 'CD8T')], DN_2, ETP_DN_2, 
                       NA_first_revised_additional_ETP_DN))

Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111101M-10X5'] <- 3
Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111401M-10X5'] <- 0.1
Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111701F-10X5'] <- 3
Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111702M-10X5'] <- 0.03
Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111801F-10X5'] <- 3
Hs_earlyT$age[Hs_earlyT$batch %in% 'PT2020111901F-10X5'] <- 0.17
Hs_earlyT$age[Hs_earlyT$batch %in% 'GSM4127993'] <- 1.6
Hs_earlyT$age[Hs_earlyT$batch %in% 'GSM4505165'] <- 1.9
Hs_earlyT$age[Hs_earlyT$batch %in% 'GSM4505166'] <- 0.02
Hs_earlyT$age[Hs_earlyT$batch %in% 'GSM4505168'] <- 0.035
Hs_earlyT$age[Hs_earlyT$batch %in% 'GSM4505169'] <- 2
Hs_earlyT$age[Hs_earlyT$source %in% 'Tom'] <- 2


Hs_earlyT$Age_group <- ''
Hs_earlyT$Age_group[Hs_earlyT$age <=12] <- 'Prepuberal'
Hs_earlyT$Age_group[Hs_earlyT$age >= 13 & Hs_earlyT$age <= 39 ] <- 'Adult'
Hs_earlyT$Age_group[Hs_earlyT$age >=40] <- 'Aged'

Hs_earlyT$age_3 <- ''
Hs_earlyT$age_3[Hs_earlyT$age <=12] <- '<=12'
Hs_earlyT$age_3[Hs_earlyT$age >= 13 & Hs_earlyT$age <= 39 ] <- '13-39'
Hs_earlyT$age_3[Hs_earlyT$age >=40] <- '>=40'

Hs_earlyT_list <- SplitObject(Hs_earlyT, split.by = 'source')
Hs_earlyT_list <- lapply(Hs_earlyT_list, function(x){
    x <- NormalizeData(x, verbose = F)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
})
anchor_features <- union(SelectIntegrationFeatures(Hs_earlyT_list, 2000), 
                         intersect(c(s.genes, g2m.genes),rownames(Hs_earlyT)))
Hs_earlyT_list <- lapply(Hs_earlyT_list, function(x){
    x <- ScaleData(x, features = anchor_features, verbose = FALSE)
    x <- RunPCA(x, features = anchor_features, verbose = FALSE)
})


plan('multicore', workers = 10)
options(future.globals.maxSize = 100 * 1024^3)

anchors <- FindIntegrationAnchors(object.list = Hs_earlyT_list, dims = 1:20, reduction = "rpca",
                                  anchor.features = anchor_features,verbose = F)
Hs_earlyT <- IntegrateData(anchorset = anchors, dims = 1:20,verbose = F)

plan('sequential')

DefaultAssay(Hs_earlyT) <- "integrated"
Hs_earlyT <- CellCycleScoring(Hs_earlyT, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
# plan('multicore', workers = 16)
# options(future.globals.maxSize = 20 * 1024^3)
Hs_earlyT <- ScaleData(Hs_earlyT)#, vars.to.regress = c('S.Score', 'G2M.Score'))
# plan('sequential')

Hs_earlyT <- RunPCA(Hs_earlyT, npcs = 50, verbose =F)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT, ndims = 50, reduction = 'pca')

Hs_earlyT <- RunUMAP(Hs_earlyT, reduction = "pca", dims = 1:15, verbose = F)


options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT) <- Hs_earlyT$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT, reduction = 'umap', label=T)
DimPlot(Hs_earlyT, reduction = 'umap', label=T, group.by = 'Phase')

options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT) <- Hs_earlyT$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT, reduction = 'umap', label=T)
DimPlot(Hs_earlyT, reduction = 'umap', label=T, group.by = 'source')




options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )


options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('CD79B'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('LYL1'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT, features =c ('KLRB1', 'IGLL1', 'PRSS57', 'IGFBP7'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('RORC'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('IL7R'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('CD4', 'CD40LG', 'CD8A', 'CD8B'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('PTPRC', 'CD7', 'CD1A','CD34'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('FOS'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('NR4A1'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('MME', 'HIVEP3', 'TRBV7-2', 'TRBC2'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('RAG1', 'RAG2', 'TCF7', 'BCL11B'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('CD34'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('HOXA9'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT) <- 'RNA'
FeaturePlot(Hs_earlyT, features = c('TRDV1'), 
            cols = c('lightgrey', 'red'))






DefaultAssay(Hs_earlyT) <- 'integrated'
Hs_earlyT <- FindNeighbors(Hs_earlyT, reduction = "pca", dims = 1:15, verbose = F)


plan('multicore', workers = 16)
options(future.globals.maxSize = 50 * 1024^3)
Hs_earlyT <- FindClusters(Hs_earlyT, resolution = seq(0.1,1, 0.05))
plan('sequential')

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 







Hs_earlyT_S1 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('3')]

DefaultAssay(Hs_earlyT_S1) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S1, ndims = 50, reduction = 'pca')


options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT_S1) <- Hs_earlyT_S1$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT_S1, reduction = 'umap', label=T)
DimPlot(Hs_earlyT_S1, reduction = 'umap', label=T, group.by = 'Phase')

options(repr.plot.width=8, repr.plot.height=8)
#Idents(Hs_earlyT_S1) <- Hs_earlyT_S1$`RNA_snn_res.0.5`
#DimPlot(Hs_earlyT_S1, reduction = 'umap', label=T)
DimPlot(Hs_earlyT_S1, reduction = 'umap', label=T, group.by = 'source')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('IGLL1', 'IGFBP7', 'LYL1', 'KLRB1'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('RAG1', 'RAG2', 'TCF7', 'BCL11B'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=8)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('MCM4', 'E2F1'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('HOXA9'), 
            cols = c('lightgrey', 'red'))





options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))

DefaultAssay(Hs_earlyT_S1) <- 'integrated'
Hs_earlyT_S1 <- FindNeighbors(Hs_earlyT_S1, reduction = 'pca', dims = 1:10, verbose = F)

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S1 <- FindClusters(Hs_earlyT_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')


###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S1, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}



Hs_earlyT_S1_S1 <- Hs_earlyT_S1[,Hs_earlyT_S1$`integrated_snn_res.1` %in% c('2')]

DefaultAssay(Hs_earlyT_S1_S1) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S1_S1, ndims = 50, reduction = 'pca')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S1, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S1, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S1, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S1, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S1, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))

DefaultAssay(Hs_earlyT_S1_S1) <- 'integrated'
Hs_earlyT_S1_S1 <- FindNeighbors(Hs_earlyT_S1_S1, reduction = 'pca', dims = 1:10, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S1_S1 <- FindClusters(Hs_earlyT_S1_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')


###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S1_S1, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S1_S1$Identity <- as.character(Hs_earlyT_S1_S1$`integrated_snn_res.0.8`)
# Hs_earlyT_S1_S1$Identity[Hs_earlyT_S1_S1$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S1_S1$Identity[!Hs_earlyT_S1_S1$Identity %in% c('0','1','2','6')]<- 'TP2'#'TP1'
Hs_earlyT_S1_S1$Identity[Hs_earlyT_S1_S1$Identity %in% c('0','1','2','6')]<- 'TP1'#'TP2'
# Hs_earlyT_S1_S1$Identity[Hs_earlyT_S1_S1$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S1_S1$Identity[Hs_earlyT_S1_S1$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S1_S1) <- Hs_earlyT_S1_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S1_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))




Hs_earlyT_S1_S2 <- Hs_earlyT_S1[,Hs_earlyT_S1$`integrated_snn_res.1` %in% c('7')]

DefaultAssay(Hs_earlyT_S1_S2) <- 'integrated'


options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S1_S2, ndims = 50, reduction = 'pca')



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S2, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S2, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S2, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S2, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S2, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S1_S2) <- 'integrated'
Hs_earlyT_S1_S2 <- FindNeighbors(Hs_earlyT_S1_S2, reduction = 'pca', dims = 1:10, verbose = F)



# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S1_S2 <- FindClusters(Hs_earlyT_S1_S2, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')



###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S1_S2, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S1_S2$Identity <- as.character(Hs_earlyT_S1_S2$`integrated_snn_res.0.85`)
# Hs_earlyT_S1_S2$Identity[Hs_earlyT_S1_S2$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S1_S2$Identity[!Hs_earlyT_S1_S2$Identity %in% c('0','3')]<- 'TP1'#'TP1'
Hs_earlyT_S1_S2$Identity[Hs_earlyT_S1_S2$Identity %in% c('0','3')]<- 'TP2'#'TP2'
# Hs_earlyT_S1_S2$Identity[Hs_earlyT_S1_S2$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S1_S2$Identity[Hs_earlyT_S1_S2$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S1_S2) <- Hs_earlyT_S1_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S1_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))




Hs_earlyT_S1_S3 <- Hs_earlyT_S1[,Hs_earlyT_S1$`integrated_snn_res.1` %in% c('3')]

DefaultAssay(Hs_earlyT_S1_S3) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S1_S3, ndims = 50, reduction = 'pca')



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S3, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S3, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S3, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S1_S3, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S1_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S1_S3, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))



DefaultAssay(Hs_earlyT_S1_S3) <- 'integrated'
Hs_earlyT_S1_S3 <- FindNeighbors(Hs_earlyT_S1_S3, reduction = 'pca', dims = 1:10, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S1_S3 <- FindClusters(Hs_earlyT_S1_S3, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')


###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S1_S3, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S1_S3$Identity <- as.character(Hs_earlyT_S1_S3$`integrated_snn_res.1`)
# Hs_earlyT_S1_S3$Identity[Hs_earlyT_S1_S3$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S1_S3$Identity[!Hs_earlyT_S1_S3$Identity %in% c('2','3','4')]<- 'TP1'#'TP1'
Hs_earlyT_S1_S3$Identity[Hs_earlyT_S1_S3$Identity %in% c('2','3','4')]<- 'TP2'#'TP2'
# Hs_earlyT_S1_S3$Identity[Hs_earlyT_S1_S3$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S1_S3$Identity[Hs_earlyT_S1_S3$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S1_S3) <- Hs_earlyT_S1_S3$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S1_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))




Hs_earlyT_S1$Identity <- as.character(Hs_earlyT_S1$`integrated_snn_res.1`)
# Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S1$Identity[!Hs_earlyT_S1$Identity %in% c('2','3','6','7','8')]<- 'TP2'#'TP1'
Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('6','8')]<- 'TP1'#'TP2'
# Hs_earlyT_S1$Identity[Hs_earlyT_S1$Identity %in% c('3')]<- 'DN_P'


Idents(Hs_earlyT_S1) <- Hs_earlyT_S1$Identity
Hs_earlyT_S1 <- Idents_proj(Hs_earlyT_S1_S1, Hs_earlyT_S1)
Hs_earlyT_S1 <- Idents_proj(Hs_earlyT_S1_S2, Hs_earlyT_S1)
Hs_earlyT_S1 <- Idents_proj(Hs_earlyT_S1_S3, Hs_earlyT_S1)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))






Hs_earlyT_S2 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('2')]

DefaultAssay(Hs_earlyT_S2) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S2, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S2, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S2, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S2) <- 'integrated'
Hs_earlyT_S2 <- FindNeighbors(Hs_earlyT_S2, reduction = 'pca', dims = 1:10, verbose = F)

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S2 <- FindClusters(Hs_earlyT_S2, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S2, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S2_S1 <- Hs_earlyT_S2[, Hs_earlyT_S2$integrated_snn_res.0.65%in% c('5')]
DefaultAssay(Hs_earlyT_S2_S1) <- 'integrated'
Hs_earlyT_S2_S1 <- FindNeighbors(Hs_earlyT_S2_S1, reduction = 'pca', verbose = F)

# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S2_S1 <- FindClusters(Hs_earlyT_S2_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S2_S1, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S2_S1$Identity <- as.character(Hs_earlyT_S2_S1$`integrated_snn_res.0.75`)
# Hs_earlyT_S2_S1$Identity[Hs_earlyT_S2_S1$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S2_S1$Identity[!Hs_earlyT_S2_S1$Identity %in% c('1','2','6')]<- 'TP2'#'TP1'
Hs_earlyT_S2_S1$Identity[Hs_earlyT_S2_S1$Identity %in% c('1','2','6')]<- 'TP1'#'TP2'
# Hs_earlyT_S2_S1$Identity[Hs_earlyT_S2_S1$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S2_S1$Identity[Hs_earlyT_S2_S1$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S2_S1) <- Hs_earlyT_S2_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S2_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S2_S2 <- Hs_earlyT_S2[,Hs_earlyT_S2$`integrated_snn_res.0.65` %in% c('4')]

DefaultAssay(Hs_earlyT_S2_S2) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S2_S2, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S2_S2, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S2, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S2, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S2, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S2_S2, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S2_S2) <- 'integrated'
Hs_earlyT_S2_S2 <- FindNeighbors(Hs_earlyT_S2_S2, reduction = 'pca', dims = 1:10, verbose = F)



# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S2_S2 <- FindClusters(Hs_earlyT_S2_S2, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')


###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S2_S2, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S2_S2$Identity <- as.character(Hs_earlyT_S2_S2$`integrated_snn_res.0.4`)
# Hs_earlyT_S2_S2$Identity[Hs_earlyT_S2_S2$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S2_S2$Identity[!Hs_earlyT_S2_S2$Identity %in% c('2','3','4')]<- 'TP2'#'TP1'
Hs_earlyT_S2_S2$Identity[Hs_earlyT_S2_S2$Identity %in% c('3','4')]<- 'TP1'#'TP2'
Hs_earlyT_S2_S2$Identity[Hs_earlyT_S2_S2$Identity %in% c('2')]<- 'DN_P'
# Hs_earlyT_S2_S2$Identity[Hs_earlyT_S2_S2$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S2_S2) <- Hs_earlyT_S2_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S2_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S2_S3 <- Hs_earlyT_S2[,Hs_earlyT_S2$`integrated_snn_res.0.65` %in% c('3','9')]

DefaultAssay(Hs_earlyT_S2_S3) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S2_S3, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S2_S3, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S3, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S3, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S2_S3, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S2_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S2_S3, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S2_S3) <- 'integrated'
Hs_earlyT_S2_S3 <- FindNeighbors(Hs_earlyT_S2_S3, reduction = 'pca', dims = 1:10, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S2_S3 <- FindClusters(Hs_earlyT_S2_S3, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')



###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 0.5, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S2_S3, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S2_S3$Identity <- as.character(Hs_earlyT_S2_S3$`integrated_snn_res.0.35`)
# Hs_earlyT_S2_S3$Identity[Hs_earlyT_S2_S3$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S2_S3$Identity[Hs_earlyT_S2_S3$Identity %in% c('2','3')]<- 'TP2'#'TP1'
Hs_earlyT_S2_S3$Identity[Hs_earlyT_S2_S3$Identity %in% c('0','1')]<- 'DN_P'#'TP2'
# Hs_earlyT_S2_S3$Identity[Hs_earlyT_S2_S3$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S2_S3$Identity[Hs_earlyT_S2_S3$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S2_S3) <- Hs_earlyT_S2_S3$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S2_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))




Hs_earlyT_S2$Identity <- as.character(Hs_earlyT_S2$`integrated_snn_res.0.65`)
Hs_earlyT_S2$Identity[Hs_earlyT_S2$Identity %in% c('7')]<- 'TP2'
Hs_earlyT_S2$Identity[Hs_earlyT_S2$Identity %in% c('0','1','2','6','8','10','11')]<- 'DN_P'



Idents(Hs_earlyT_S2) <- Hs_earlyT_S2$Identity
Hs_earlyT_S2 <- Idents_proj(Hs_earlyT_S2_S1, Hs_earlyT_S2)
Hs_earlyT_S2 <- Idents_proj(Hs_earlyT_S2_S2, Hs_earlyT_S2)
Hs_earlyT_S2 <- Idents_proj(Hs_earlyT_S2_S3, Hs_earlyT_S2)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))







Hs_earlyT_S3 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('1')]

DefaultAssay(Hs_earlyT_S3) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S3, ndims = 50, reduction = 'pca')


options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(Hs_earlyT_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S3, features = c('TRDV1'), 
            cols = c('lightgrey', 'red'))



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S3, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S3, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S3) <- 'integrated'
Hs_earlyT_S3 <- FindNeighbors(Hs_earlyT_S3, reduction = 'pca', dims = 1:10, verbose = F)




options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S3 <- FindClusters(Hs_earlyT_S3, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')


options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S3, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S3_S1 <- Hs_earlyT_S3[,Hs_earlyT_S3$`integrated_snn_res.0.35` %in% c('3')]

DefaultAssay(Hs_earlyT_S3_S1) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S3_S1, ndims = 50, reduction = 'pca')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S1, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S1, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S1, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S1, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S1, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))



DefaultAssay(Hs_earlyT_S3_S1) <- 'integrated'
Hs_earlyT_S3_S1 <- FindNeighbors(Hs_earlyT_S3_S1, reduction = 'pca', dims = 1:10, verbose = F)



options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S3_S1 <- FindClusters(Hs_earlyT_S3_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')



###输出不同分辨率的聚类结果
options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 0.3, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S3_S1, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S3_S1$Identity <- as.character(Hs_earlyT_S3_S1$`integrated_snn_res.0.2`)
# Hs_earlyT_S3_S1$Identity[Hs_earlyT_S3_S1$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S3_S1$Identity[!Hs_earlyT_S3_S1$Identity %in% c('0','2')]<- 'DN_P'#'TP1'
Hs_earlyT_S3_S1$Identity[Hs_earlyT_S3_S1$Identity %in% c('0','2')]<- 'TP2'#'TP2'
# Hs_earlyT_S3_S1$Identity[Hs_earlyT_S3_S1$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S3_S1$Identity[Hs_earlyT_S3_S1$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S3_S1) <- Hs_earlyT_S3_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S3_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S3_S2 <- Hs_earlyT_S3[,Hs_earlyT_S3$`integrated_snn_res.0.35` %in% c('4')]

DefaultAssay(Hs_earlyT_S3_S2) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S3_S2, ndims = 50, reduction = 'pca')



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S2, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S2, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S2, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S2, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S2, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S3_S2) <- 'integrated'
Hs_earlyT_S3_S2 <- FindNeighbors(Hs_earlyT_S3_S2, reduction = 'pca', dims = 1:10, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S3_S2 <- FindClusters(Hs_earlyT_S3_S2, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')


options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 0.3, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S3_S2, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S3_S2$Identity <- as.character(Hs_earlyT_S3_S2$`integrated_snn_res.0.2`)
# Hs_earlyT_S3_S2$Identity[Hs_earlyT_S3_S2$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S3_S2$Identity[!Hs_earlyT_S3_S2$Identity %in% c('1')]<- 'DN_P'#'TP1'
Hs_earlyT_S3_S2$Identity[Hs_earlyT_S3_S2$Identity %in% c('1')]<- 'TP2'#'TP2'
# Hs_earlyT_S3_S2$Identity[Hs_earlyT_S3_S2$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S3_S2$Identity[Hs_earlyT_S3_S2$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S3_S2) <- Hs_earlyT_S3_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S3_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S3_S3 <- Hs_earlyT_S3[,Hs_earlyT_S3$`integrated_snn_res.0.35` %in% c('1')]

DefaultAssay(Hs_earlyT_S3_S3) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S3_S3, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S3, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S3, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S3, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S3_S3, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S3_S3) <- 'RNA'
FeaturePlot(Hs_earlyT_S3_S3, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S3_S3) <- 'integrated'
Hs_earlyT_S3_S3 <- FindNeighbors(Hs_earlyT_S3_S3, reduction = 'pca', dims = 1:10, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Hs_earlyT_S3_S3 <- FindClusters(Hs_earlyT_S3_S3, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 0.3, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S3_S3, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S3_S3$Identity <- as.character(Hs_earlyT_S3_S3$`integrated_snn_res.0.3`)
# Hs_earlyT_S3_S3$Identity[Hs_earlyT_S3_S3$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S3_S3$Identity[!Hs_earlyT_S3_S3$Identity %in% c('2')]<- 'DN_P'#'TP1'
Hs_earlyT_S3_S3$Identity[Hs_earlyT_S3_S3$Identity %in% c('2')]<- 'TP2'#'TP2'
# Hs_earlyT_S3_S3$Identity[Hs_earlyT_S3_S3$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S3_S3$Identity[Hs_earlyT_S3_S3$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S3_S3) <- Hs_earlyT_S3_S3$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S3_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S3$Identity <- as.character(Hs_earlyT_S3$`integrated_snn_res.0.35`)
# Hs_earlyT_S3$Identity[Hs_earlyT_S3$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S3$Identity[!Hs_earlyT_S3$Identity %in% c('1','2','3','4','6')]<- 'DN_P'#'TP1'
Hs_earlyT_S3$Identity[Hs_earlyT_S3$Identity %in% c('2')]<- 'TP2'#'TP2'
Hs_earlyT_S3$Identity[Hs_earlyT_S3$Identity %in% c('6')]<- 'γδT_P'#'TP2'
# Hs_earlyT_S3$Identity[Hs_earlyT_S3$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S3$Identity[Hs_earlyT_S3$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S3) <- Hs_earlyT_S3$Identity
Hs_earlyT_S3 <- Idents_proj(Hs_earlyT_S3_S1, Hs_earlyT_S3)
Hs_earlyT_S3 <- Idents_proj(Hs_earlyT_S3_S2, Hs_earlyT_S3)
Hs_earlyT_S3 <- Idents_proj(Hs_earlyT_S3_S3, Hs_earlyT_S3)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S4 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('4')]

DefaultAssay(Hs_earlyT_S4) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S4, ndims = 50, reduction = 'pca')



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S4) <- 'RNA'
FeaturePlot(Hs_earlyT_S4, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S4, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S4, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S4, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S4) <- 'RNA'
FeaturePlot(Hs_earlyT_S4, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))

DefaultAssay(Hs_earlyT_S4) <- 'integrated'
Hs_earlyT_S4 <- FindNeighbors(Hs_earlyT_S4, reduction = 'pca', dims = 1:10, verbose = F)


options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S4 <- FindClusters(Hs_earlyT_S4, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S4, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S4$Identity <- as.character(Hs_earlyT_S4$`integrated_snn_res.0.4`)
# Hs_earlyT_S4$Identity[Hs_earlyT_S4$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S4$Identity[!Hs_earlyT_S4$Identity %in% c('3')]<- 'DN_P'#'TP1'
Hs_earlyT_S4$Identity[Hs_earlyT_S4$Identity %in% c('3')]<- 'TP2'#'TP2'
# Hs_earlyT_S4$Identity[Hs_earlyT_S4$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S4$Identity[Hs_earlyT_S4$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S4) <- Hs_earlyT_S4$Identity

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S4, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))






Hs_earlyT_S5 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('0')]

DefaultAssay(Hs_earlyT_S5) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S5, ndims = 50, reduction = 'pca')



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S5) <- 'RNA'
FeaturePlot(Hs_earlyT_S5, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S5, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S5, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S5, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S5) <- 'RNA'
FeaturePlot(Hs_earlyT_S5, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S5) <- 'integrated'
Hs_earlyT_S5 <- FindNeighbors(Hs_earlyT_S5, reduction = 'pca', dims = 1:10, verbose = F)


options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S5 <- FindClusters(Hs_earlyT_S5, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S5, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S5_S1 <- Hs_earlyT_S5[,Hs_earlyT_S5$`integrated_snn_res.0.3` %in% c('0')]

DefaultAssay(Hs_earlyT_S5_S1) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S5_S1, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S5_S1) <- 'RNA'
FeaturePlot(Hs_earlyT_S5_S1, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


DefaultAssay(Hs_earlyT_S5_S1) <- 'integrated'
Hs_earlyT_S5_S1 <- FindNeighbors(Hs_earlyT_S5_S1, reduction = 'pca', dims = 1:10, verbose = F)



options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S5_S1 <- FindClusters(Hs_earlyT_S5_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')


options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S5_S1, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}


Hs_earlyT_S5_S1$Identity <- as.character(Hs_earlyT_S5_S1$`integrated_snn_res.0.4`)
# Hs_earlyT_S5_S1$Identity[Hs_earlyT_S5_S1$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S5_S1$Identity[!Hs_earlyT_S5_S1$Identity %in% c('1')]<- 'DN_Q'#'TP1'
Hs_earlyT_S5_S1$Identity[Hs_earlyT_S5_S1$Identity %in% c('1')]<- 'TP2'#'TP2'
# Hs_earlyT_S5_S1$Identity[Hs_earlyT_S5_S1$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S5_S1$Identity[Hs_earlyT_S5_S1$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S5_S1) <- Hs_earlyT_S5_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S5_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))




Hs_earlyT_S5_S2 <- Hs_earlyT_S5[,Hs_earlyT_S5$`integrated_snn_res.0.3` %in% c('2')]

DefaultAssay(Hs_earlyT_S5_S2) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S5_S2, ndims = 50, reduction = 'pca')


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S5_S2) <- 'RNA'
FeaturePlot(Hs_earlyT_S5_S2, features = c('TRBV7-2', 'TRDC', 'TRGC1', 'RAG2'), 
            cols = c('lightgrey', 'red'))

DefaultAssay(Hs_earlyT_S5_S2) <- 'integrated'
Hs_earlyT_S5_S2 <- FindNeighbors(Hs_earlyT_S5_S2, reduction = 'pca', dims = 1:20, verbose = F)


options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 20)
Hs_earlyT_S5_S2 <- FindClusters(Hs_earlyT_S5_S2, resolution = seq(0.2, 2, 0.05), verbose = F)
plan('sequential')


options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 2, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S5_S2, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Hs_earlyT_S5_S2$Identity <- as.character(Hs_earlyT_S5_S2$`integrated_snn_res.0.3`)
# Hs_earlyT_S5_S2$Identity[Hs_earlyT_S5_S2$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S5_S2$Identity[!Hs_earlyT_S5_S2$Identity %in% c('4')]<- 'DN_Q'#'TP1'
Hs_earlyT_S5_S2$Identity[Hs_earlyT_S5_S2$Identity %in% c('4')]<- 'DN4'#'TP2'
# Hs_earlyT_S5_S2$Identity[Hs_earlyT_S5_S2$Identity %in% c('3')]<- 'DN_P'
# Hs_earlyT_S5_S2$Identity[Hs_earlyT_S5_S2$Identity %in% c('11')]<- 'DN_P'

Idents(Hs_earlyT_S5_S2) <- Hs_earlyT_S5_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S5_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT_S5$Identity <- as.character(Hs_earlyT_S5$`integrated_snn_res.0.3`)
# Hs_earlyT_S5$Identity[Hs_earlyT_S5$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
Hs_earlyT_S5$Identity[Hs_earlyT_S5$Identity %in% c('1','3','4')]<- 'DN_Q'#'TP1'

Idents(Hs_earlyT_S5) <- Hs_earlyT_S5$Identity
Hs_earlyT_S5 <- Idents_proj(Hs_earlyT_S5_S1, Hs_earlyT_S5)
Hs_earlyT_S5 <- Idents_proj(Hs_earlyT_S5_S2, Hs_earlyT_S5)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S5, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))






Hs_earlyT_S6 <- Hs_earlyT[,Hs_earlyT$`integrated_snn_res.0.2` %in% c('7')]

DefaultAssay(Hs_earlyT_S6) <- 'integrated'

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hs_earlyT_S6, ndims = 50, reduction = 'pca')




options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S6) <- 'RNA'
FeaturePlot(Hs_earlyT_S6, features = c('CD1A', 'CD7', 'CD2','CD34', 'MME', 'CD44', 'PTCRA'), 
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S6, cols = c('lightgrey', 'red'), 
            features = c('HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', 'MS4A1', 'BCL11A') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S6, cols = c('lightgrey', 'red'), 
            features = c('MME', 'CD44', 'HOXA9', 'MEIS1') )

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S6, features =c ('DTX1', 'LYL1', 'BCL11B', 'MPO'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(Hs_earlyT_S6) <- 'RNA'
FeaturePlot(Hs_earlyT_S6, features = c('MPO', 'SPI1', 'LYL1', 'CD7'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Hs_earlyT_S6, features =c ('KLRB1', 'IGLL1', 'PRSS57', 'IGFBP7'), 
            cols = c('lightgrey', 'red'))

DefaultAssay(Hs_earlyT_S6) <- 'integrated'
Hs_earlyT_S6 <- FindNeighbors(Hs_earlyT_S6, reduction = 'pca', dims = 1:10, verbose = F)


options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
Hs_earlyT_S6 <- FindClusters(Hs_earlyT_S6, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')


options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(Hs_earlyT_S6, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}



Hs_earlyT_S6$Identity <- as.character(Hs_earlyT_S6$`integrated_snn_res.0.2`)

Hs_earlyT_S6$Identity[Hs_earlyT_S6$Identity %in% c('0')]<- 'ETP'
Hs_earlyT_S6$Identity[Hs_earlyT_S6$Identity %in% c('1','2')]<- 'TP1'#TP2


Idents(Hs_earlyT_S6) <- Hs_earlyT_S6$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT_S6, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


Hs_earlyT$Identity <- as.character(Hs_earlyT$`integrated_snn_res.0.2`)
Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('5')]<- 'DN_P'
Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('6','8')]<- 'DN_Q'
# Hs_earlyT$Identity[Hs_earlyT$Identity %in% c('8')]<- 'ETP'


Idents(Hs_earlyT) <- Hs_earlyT$Identity
Hs_earlyT <- Idents_proj(Hs_earlyT_S1, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S2, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S3, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S4, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S5, Hs_earlyT)
Hs_earlyT <- Idents_proj(Hs_earlyT_S6, Hs_earlyT)
# Hs_earlyT <- Idents_proj(Hs_earlyT_S7, Hs_earlyT)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))



levels(Hs_earlyT)

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
markers_DN_Q <- FindMarkers(ETP_DN_2, pseudocount.use = 0, logfc.threshold = log(1.5), `ident.2` = 'DN_Q2',
                            `ident.1` = 'DN_Q1')
plan('sequential')
saveRDS(markers_DN_Q, file = 'markers_DN_Q.rds')


DN_Q <- Hs_earlyT[,Idents(Hs_earlyT) %in%c('DN_Q')]



DefaultAssay(DN_Q) <- 'RNA'

DN_Q$RNA@var.features <- rownames(markers_DN_Q)[markers_DN_Q$p_val_adj < 0.01 & 
                                                 abs(markers_DN_Q$avg_logFC) > log(1.5)]
#FindVariableFeatures(DN_Q, selection.method = "vst", nfeatures = 2000, verbose = F)
DN_Q <- ScaleData(DN_Q, verbose = F, features = DN_Q$RNA@var.features)
DN_Q <- RunPCA(DN_Q, npcs = 50, verbose = FALSE, features = DN_Q$RNA@var.features)
###矫正批次效应


DN_Q <- RunHarmony(DN_Q, c('batch', 'source'))



options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(DN_Q, ndims = 50, reduction = 'harmony')



DN_Q <- RunUMAP(DN_Q, reduction = "harmony", dims = 1:20, verbose  = F)


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(DN_Q) <- 'RNA'
FeaturePlot(DN_Q, features = c('MME','HIVEP3','PTCRA','MAL','RAG1','RAG2','BCL11B',
'TRBC1','TRBC2','TRBV7-2','TCF12','IL7R','DNTT','NOTCH3'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))


options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(DN_Q) <- 'RNA'
FeaturePlot(DN_Q, features = c('CD1A','TCF7','CD2','LEF1','AEBP1','LCK','ZAP70',
              'TRGC2','ICAM2',
              'TRDC',
'TRDV1','TRGV4','ID3'), min.cutoff = 'q9',
            cols = c('lightgrey', 'red'))




options(repr.plot.width=8, repr.plot.height=6)
DimPlot(DN_Q, reduction = 'umap',label = T, group.by = 'source',label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


# options(repr.plot.width=8, repr.plot.height=6)
# DimPlot(DN_Q, reduction = 'umap',label = T, group.by = 'old_id',label.size = 20/.pt)+
# scale_color_manual(values = color_used)+
# theme(legend.text = element_text(size = 20))+
# guides(color = guide_legend(override.aes = list(size = 7)))


DefaultAssay(Hs_earlyT_S5) <- 'RNA'
DN_Q <- FindNeighbors(DN_Q, reduction = 'harmony', dims = 1:20, verbose = F)


options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
DN_Q <- FindClusters(DN_Q, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('RNA_snn_res.', i)
    p <- DimPlot(DN_Q, reduction = 'umap',label = T, label.size = 20/.pt,# raster = F,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

DN_Q$Identity <- as.character(DN_Q$`RNA_snn_res.0.55`)
# DN_Q$Identity[DN_Q$Identity %in% c('6', '10', '12', 
#                                                   '0', '3','4','8')]<- 'TP'
DN_Q$Identity[!DN_Q$Identity %in% c('2','6')]<- 'DN_Q1'#'TP1'
DN_Q$Identity[DN_Q$Identity %in% c('2','6')]<- 'DN_Q2'#'TP1'

Idents(DN_Q) <- DN_Q$Identity

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(DN_Q, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))



Hs_earlyT <- Idents_proj(DN_Q, Hs_earlyT)
levels(Hs_earlyT) <- c('ETP', 'TP1', 'TP2', 'DN_P', 'DN_Q1', 'DN_Q2', 'DN4', 
                      'γδT_P')
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Hs_earlyT, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))








saveRDS(Idents(Hs_earlyT), file = 'Hs_earlyT_newid.rds')

saveRDS(Hs_earlyT, file = 'Hs_earlyT.rds')

save(Hs_earlyT_S1, Hs_earlyT_S1_S1, Hs_earlyT_S1_S2, Hs_earlyT_S1_S3,
     Hs_earlyT_S2, Hs_earlyT_S2_S1, Hs_earlyT_S2_S2, Hs_earlyT_S2_S3,
     Hs_earlyT_S3, Hs_earlyT_S3_S1, Hs_earlyT_S3_S2, Hs_earlyT_S3_S3,
     Hs_earlyT_S4,
     Hs_earlyT_S5, Hs_earlyT_S5_S1, Hs_earlyT_S5_S2,
     Hs_earlyT_S6, 
     file = 'Hs_earlyT_subsets.RData')

saveRDS(DN_Q, file = 'DN_Q.rds')













ETP_TP <- Hs_earlyT[,grep('TP',Idents(Hs_earlyT))]

DefaultAssay(ETP_TP) <- 'integrated'
ETP_TP<- RunUMAP(ETP_TP, reduction = 'pca', dims = 1:6, verbose = F, min.dist = 0.3)



options(repr.plot.width=8, repr.plot.height=6)
p <- DimPlot(ETP_TP, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))

p

###Extended Data Fig 2a
DefaultAssay(ETP_TP) <- 'RNA'
gene_order<-unique(c(
'MME', 'CD44', 
'HOPX', 'BAALC', 'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 'IL3RA', #'BCL11A' ,
'HOXA9', 'MEIS1' ,'KLRB1', 'LYL1',  'MPO',# 'SELL',
 'IGLL1', 'PRSS57', 'IGFBP7','CD1A', 'CD7', 'CD2','BCL11B','PTCRA'
            ))

p <- DotPlot(ETP_TP, features = rev(gene_order))

p1 <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_y_discrete(limits = 
                 rev(levels(ETP_TP)), position = 'right', 
                label = function(x){
    x <- gsub('_', '-', x)
})+ scale_size(range=c(0.5, 6.5))+
scale_x_discrete(position = "top") +
theme(axis.title = element_blank(), #legend.direction = c('horizontal'),
      legend.position = 'bottom',legend.box  ='horizontal',
      legend.key.width = unit(0.5, 'cm'),
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.y = element_text(hjust = 1, #vjust = 0.5,
                                 size = 13.75, family = 'sans', face = 'plain'),
     axis.text.x = element_text(hjust =0, angle =90,size = 13.75, family = 'sans', face = 'italic', vjust = 0), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'), breaks = c(-1,0,1))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'), ncol = 2), 
       fill = guide_colorbar(title = "Expression")) 
options(repr.plot.width=5.5, repr.plot.height=3)
p1




saveRDS(ETP_TP, file = 'ETP_TP.rds')



ETP_DN <- Hs_earlyT[,Hs_earlyT$source %in% 'Hu' ]

ETP_DN$combined_idents <- as.character(Idents(ETP_DN))
ETP_DN$combined_idents[grep('ETP|TP', ETP_DN$combined_idents)] <- 'ETP-TP'
ETP_DN$combined_idents <- factor(ETP_DN$combined_idents, levels = c('ETP-TP','DN_P','DN_Q1','DN_Q2','DN4','γδT_P'))

###Extended Data Fig 2a
gene_order<-c('CD34','NFE2', 'BTK', 'SPINK2', 'HHEX', 'BAALC', 'IGLL1', 'IFITM1','TRGC1',
                  'CCNB2','CDC20','CDC6','CENPA',
                  'PCNA','MKI67','CDK1',#'HES1',
              'HES4','MME','HIVEP3','PTCRA','MAL','RAG1','RAG2','BCL11B',
'TRBC1','TRBC2','TRBV7-2','TCF12','IL7R','DNTT','NOTCH3','TSHR',#'CXCR4',
              'TRAC','CD4','CD8B','RORC',
'AQP3','CD8A','CD38','CD1A','TCF7','CD2','LEF1','AEBP1','LCK','ZAP70','GATA3','CCR9','ITGB2',#'ID2',
                  'EGR1',
              'TRGC2','ICAM2',
              'TRDC',
'TRDV1','TRGV4','ID3','CD27')

DefaultAssay(ETP_DN) <- 'RNA'
p <- DotPlot(ETP_DN, features = rev(gene_order), group.by = 'combined_idents')

p1 <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_y_discrete(limits =rev, position = 'right', 
                label = function(x){
    x <- gsub('_', '-', x)
})+ scale_size(range=c(0.5, 6.5))+
scale_x_discrete(position = "top") +
theme(axis.title = element_blank(), 
      legend.position = 'right',legend.box  ='horizontal',
      legend.text = element_text(size = 17.5, family = 'sans'), 
      legend.title = element_text(size = 17.5, family = 'sans'), 
      axis.text.y = element_text(hjust = 1, #vjust = 0.5,
                                 size = 17.5, family = 'sans', face = 'plain'),
     axis.text.x = element_text(hjust =0, angle =90,size = 17.5, family = 'sans', face = 'italic', vjust = 0.5), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression")) 
options(repr.plot.width=17.5, repr.plot.height=4)
p1


DefaultAssay(ETP_DN) <- 'integrated'
ETP_DN<- RunUMAP(ETP_DN, reduction = 'pca', dims = 1:10, verbose = F, min.dist = 0.3)




ETP_DN_tcr_idents <- as.character(Idents(ETP_DN))
names(ETP_DN_tcr_idents) <- names(Idents(ETP_DN))

saveRDS(ETP_DN_tcr_idents, file = 'ETP_DN_tcr_idents.rds')



saveRDS(ETP_DN, file = 'ETP_DN_Hu.rds')


