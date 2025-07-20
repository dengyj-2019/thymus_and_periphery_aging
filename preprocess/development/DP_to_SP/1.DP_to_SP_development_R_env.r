library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(SCopeLoomR)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/T_dev.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/unconventional/unconventional.RData')

library(harmony)
library(ggpubr)


library(openxlsx)

library(ggraph)
library(tidygraph)

ETP_SP_dev <- merge(ETP_DN_2, list(T_dev_3, unconventional_3))

ETP_SP_dev$identity <- Idents(ETP_SP_dev)
ETP_SP_dev$identity <- as.character(ETP_SP_dev$identity)

ETP_SP_dev$identity_2 <- ETP_SP_dev$identity
ETP_SP_dev$identity_2[ETP_SP_dev$identity %in% c('ETP', 'TP1', 'TP2')] <- 'ETP/TP'
ETP_SP_dev$identity_2[ETP_SP_dev$identity %in% c('DN_Q1', 'DN_Q2')] <- 'DN_Q'
ETP_SP_dev$identity_2[ETP_SP_dev$identity %in% c('GNG4+CD8ααT', 'ZNF683+CD8ααT', 'IFN_CD8ααT')] <- 'CD8ααT'

###down size
ETP_SP_sub_dev <- SubsetData(ETP_SP_dev, max.cells.per.ident = 6000, random.seed = 1)

Dev_1 <- ETP_SP_sub_dev[, ETP_SP_sub_dev$identity_2 %in% 
                       c('DP_Q', 'HSP_DP','DP_P','CD8+pre_T', 'CD4+pre_T','CD8ααT','NKT_like_CD8ααT',
                        'pre_Treg1', 'pre_Treg2','CD8+Naive', 'CD4+Naive', 'IFN_response', 
                        'Treg', 'Agonist_1', 'αβ_Entry', 'Agonist_2')]
Dev_1 <- NormalizeData(Dev_1)
Dev_1 <- FindVariableFeatures(Dev_1, nfeatures = 2000)
Dev_1 <- ScaleData(Dev_1, vars.to.regress = c('CC.Difference'))###reduced effect of cell cycle
Dev_1 <- RunPCA(Dev_1, npcs = 50,verbose=F, features = setdiff(Dev_1$RNA@var.features, gogenes))###reduced effect of cell cycle
Dev_1 <- RunHarmony(Dev_1, 'batch')

ElbowPlot(Dev_1, ndims = 50, reduction = 'harmony')

DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 1:2)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 3:4)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 5:6)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 7:8)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 9:10)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 11:12)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 13:14)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 15:16)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 17:18)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 19:20)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 21:22)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 23:24)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 25:26)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 27:28)
DimPlot(Dev_1, reduction = 'harmony', group.by = 'Phase',dims = 29:30)

###passed harmony space to python environment
loom <- build_loom('Dev_1.loom',dgem=Dev_1$RNA@scale.data)
close_loom(loom)
write.csv(Embeddings(Dev_1, reduction = 'harmony')[,setdiff(1:50, c())], 
          quote = F, file = 'Dev_1_harmony.csv')
write.csv(Dev_1[[]], quote = F, file = 'Dev_1_metadata.csv')





###with DP_P exclude to reduce effect of cell cycle and obtain more details of the trajectory
Dev_2 <- ETP_SP_sub_dev[, ETP_SP_sub_dev$identity_2 %in% 
                       c('DP_Q', 'HSP_DP','CD8+pre_T', 'CD4+pre_T','CD8ααT','NKT_like_CD8ααT',
                        'pre_Treg1', 'pre_Treg2','CD8+Naive', 'CD4+Naive', 'IFN_response', 
                        'Treg', 'Agonist_1', 'αβ_Entry', 'Agonist_2')]
Dev_2 <- NormalizeData(Dev_2)
Dev_2 <- FindVariableFeatures(Dev_2, nfeatures = 2000)
Dev_2 <- ScaleData(Dev_2)
Dev_2 <- RunPCA(Dev_2, npcs = 50,verbose=F)
Dev_2 <- RunHarmony(Dev_2, 'batch')

ElbowPlot(Dev_2, ndims = 50, reduction = 'harmony')

DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 1:2)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 3:4)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 5:6)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 7:8)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 9:10)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 11:12)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 13:14)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 15:16)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 17:18)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 19:20)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 21:22)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 23:24)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 25:26)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 27:28)
DimPlot(Dev_2, reduction = 'harmony', group.by = 'Phase',dims = 29:30)

###passed harmony space to python environment
loom <- build_loom('Dev_2.loom',dgem=Dev_2$RNA@scale.data)
close_loom(loom)
write.csv(Embeddings(Dev_2, reduction = 'harmony')[,setdiff(1:50, c())], 
          quote = F, file = 'Dev_2_harmony.csv')
write.csv(Dev_2[[]], quote = F, file = 'Dev_2_metadata.csv')





save(Dev_1, Dev_2, file = 'ETP_SP_Dev.RData')

save(ETP_SP_sub_dev, ETP_SP_dev, file = 'ETP_SP_dev.RData')
