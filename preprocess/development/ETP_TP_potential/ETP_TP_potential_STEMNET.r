library(Seurat)
library(STEMNET)
library(reshape2)
library(ggplot2)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

source('/data1/02.private/dengyj/analysis/thymus/plot/support//plot_helper.r')

library(patchwork)

load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/Bcell/Bcell.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN//ETP_DN.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/Myeloid/Myeloid.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/unconventional//unconventional.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/thymus_step_1.rds')


Hs_earlyT <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/Hs_earlyT.rds')

Ery  <- thymus_step_1[, Idents(thymus_step_1) == '42']
Idents(Ery) <- 'Ery'

ETP_STEMNET <- 
  merge(Ery, list(Hs_earlyT[, Idents(Hs_earlyT) %in% c('ETP','TP1','TP2','DN_P')],
                  Bcell_2[, Idents(Bcell_2) == 'NaiveB'], 
                  Myeloid_2[, Idents(Myeloid_2) %in% c('pDC', 'DC1','DC1P', 'DC2', 'Monocyte')],
                  unconventional_3[, grepl('MatureNK|Lti|ILC', Idents(unconventional_3))]))

ETP_STEMNET$identity <- Idents(ETP_STEMNET)
ETP_STEMNET <- ETP_STEMNET[rowSums(ETP_STEMNET$RNA@counts) > 0, ]

ETP_STEMNET$age_3 <- ''
ETP_STEMNET$age_3[ETP_STEMNET$age <= 12] <- '<=12'
ETP_STEMNET$age_3[ETP_STEMNET$age > 12 & ETP_STEMNET$age < 40] <- '13-39'
ETP_STEMNET$age_3[ETP_STEMNET$age >= 40] <- '>=40'

ETP_STEMNET$Age_group <- ''
ETP_STEMNET$Age_group[ETP_STEMNET$age <=12] <- 'Prepuberal'
ETP_STEMNET$Age_group[ETP_STEMNET$age >= 13 & ETP_STEMNET$age <= 39 ] <- 'Adult'
ETP_STEMNET$Age_group[ETP_STEMNET$age >=40] <- 'Aged'


ETP_STEMNET$idents <- Idents(ETP_STEMNET)
ETP_STEMNET$idents_age <- paste0(ETP_STEMNET$idents,'~', ETP_STEMNET$Age_group)
Idents(ETP_STEMNET) <- ETP_STEMNET$idents_age
ETP_STEMNET_sub <- SubsetData(ETP_STEMNET, max.cells.per.ident = 1000, random.seed = 1)

Idents(ETP_STEMNET) <- ETP_STEMNET$idents
Idents(ETP_STEMNET_sub) <- ETP_STEMNET_sub$idents

DefaultAssay(ETP_STEMNET_sub) <- 'RNA'
ETP_STEMNET_sub <- NormalizeData(ETP_STEMNET_sub)

dev_identity <- as.character(ETP_STEMNET_sub$identity)

dev_identity[grepl('DN', dev_identity)] <- 'T_lineage'
dev_identity[grepl('NaiveB', dev_identity)] <- 'B_lineage'
dev_identity[dev_identity %in% c('DC1P', 'DC1')] <- 'DC1'
dev_identity[grepl('ILC|Lti', dev_identity)] <- 'ILC'
dev_identity[dev_identity %in% c('MatureNK')] <- 'NK'
dev_identity[dev_identity %in% c('ETP','TP1', 'TP2')] <- NA

ETP_STEMNET_sub <- FindVariableFeatures(ETP_STEMNET_sub, nfeatures = 5000)

#cc_features <- ETP_STEMNET_sub$RNA@var.features[moduleColors %in% c('blue', 'turquoise')]
ETP_STEMNET_sub <- ScaleData(ETP_STEMNET_sub,verbose = FALSE, 
                         features = setdiff(ETP_STEMNET_sub$RNA@var.features, c()))
exp_matrix <- t(as.matrix(ETP_STEMNET_sub$RNA@scale.data))


result <- runSTEMNET(exp_matrix, dev_identity)



saveRDS(result, file = 'result_seed1_1000_cells.rds')

saveRDS(ETP_STEMNET, file = 'ETP_STEMNET.rds')
