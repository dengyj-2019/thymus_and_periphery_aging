library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
library(harmony)
library(SCopeLoomR)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')

T_NK_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_2.rds')

T_NK_2$identity <- Idents(T_NK_2)

PB_CD4T <- T_NK_2[, Idents(T_NK_2) %in% 
                       c('CD4+CTL', 'CD4+TEFF', 'CD4+TN','CD4+TEM', 'Treg', 
                         'CD4+TCM',
                        'CD4+RTE')]
PB_CD4T <- NormalizeData(PB_CD4T)
PB_CD4T <- FindVariableFeatures(PB_CD4T, nfeatures = 2000)
PB_CD4T <- ScaleData(PB_CD4T)
PB_CD4T <- RunPCA(PB_CD4T, npcs = 50,verbose=F)
PB_CD4T <- RunHarmony(PB_CD4T, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(PB_CD4T, ndims = 50, reduction = 'harmony')

###passed harmony space to python environment
loom <- build_loom('PB_CD4T.loom',dgem=PB_CD4T$RNA@scale.data)
close_loom(loom)
write.csv(Embeddings(PB_CD4T, reduction = 'harmony')[,setdiff(1:50, c())], 
          quote = F, file = 'PB_CD4T_harmony.csv')
write.csv(PB_CD4T[[]], quote = F, file = 'PB_CD4T_metadata.csv')





PB_CD8T <- T_NK_2[, Idents(T_NK_2) %in% 
                       c('CD8+TEM', 'CD8+TCM', 'CD8+TEFF', 'CD8+TN', 'CD8+RTE')]
PB_CD8T <- NormalizeData(PB_CD8T)
PB_CD8T <- FindVariableFeatures(PB_CD8T, nfeatures = 2000)
PB_CD8T <- ScaleData(PB_CD8T)
PB_CD8T <- RunPCA(PB_CD8T, npcs = 50,verbose=F)
PB_CD8T <- RunHarmony(PB_CD8T, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(PB_CD8T, ndims = 50, reduction = 'harmony')

loom <- build_loom('PB_CD8T.loom',dgem=PB_CD8T$RNA@scale.data)
close_loom(loom)
write.csv(Embeddings(PB_CD8T, reduction = 'harmony')[,setdiff(1:50, c())], 
          quote = F, file = 'PB_CD8T_harmony.csv')
write.csv(PB_CD8T[[]], quote = F, file = 'PB_CD8T_metadata.csv')







save(PB_CD4T, PB_CD8T, file = 'PB_Tcell_dev.RData')
