library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(harmony)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
#library(kBET)
library(openxlsx)
library(ggpubr)
library(ComplexHeatmap)

library(patchwork)

library(Matrix)

library(gam)

library(Matrix)

library(circlize)

library(ComplexHeatmap)

library(Matrix)

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

source('/data1/02.private/dengyj/analysis/mycode//plot_helper.r')

matureSP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/matureSP.rds')
T_NK_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_2.rds')

mycolors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', 
             '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F',
             '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69',
             '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863',
             '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', 
             '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

matureSP$age_3 <- ''
matureSP$age_3[matureSP$age <= 12] <- '<=12'
matureSP$age_3[matureSP$age >= 13 & matureSP$age <= 39] <- '13-39'
matureSP$age_3[matureSP$age >= 40 & matureSP$age <= 99] <- '40-99'

matureSP$tissue <- 'thymus'
T_NK_2$tissue <- 'PB'

cluster_colors <- 
 c('pre_Treg1'='#E5D2DD','CD8+Naive'='#53A85F','CD4+Memory_like'='#808040',
                    'CD8+Memory_like'='#F3B1A0',
                             'CD4+Naive'='#D6E7A3','IFN_response'='#476D87','thymus_Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3','Agonist_2'='#a83117','Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'Agonist_1' = '#a6c050', 'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'pre_Treg2' = '#9566ac' ,'Tex_like' = '#9566ac' , 'None' = 'lightgrey'   ,
       'CD8+SOX4+Naive'='#53A85F','CD8+SOX4-Naive' = "#F2C43D", 
                              'CD4+SOX4-Naive' ="#D76483",
                             'CD4+SOX4+Naive'='#D6E7A3'
                            )

PB_cluster_colors <-  c('CD8+TEM-B'='#cb50a0',
                    'CD8+TEM-K'='#db72b4',#
                    'CD8+RTE'='#276d6d','CD4+TN'='#7bacd9',
                        'CD4+CTL'='#7c75ce',#
                    'CD56_Dim_NK'='#c8af1e',#
                    'CD56_Bright_NK'='#d99333',#
                    'CD27lo_Vδ1T'='#93c88e','CD8+TN'='#DDC4DA',
                    'CD4+TEM'='#868ae3',#'#909ee7',##'
                    'PB_Treg'='#6e7cad','CD4+TCM'='#80a2e6',#'#9cb9e5',
                    'MAIT'='#b1b772',
                    'CD4+RTE'='#66c5b9',
                        'Vδ2Vγ9T'='#A0D7C9',
                    'NKT'='#e3cd82',#
                    'CD27+Vδ1T'='#9cb39b',#
                    'IFN_T'='#ec807c',
                    'T_un'='#747474','Tex'='#A1B4BB', 'CD4+TEFF' = '#55b3db',##
                    'CD8+TCM'='#E59CC4')

cluster_colors <- c(cluster_colors, PB_cluster_colors)





combined_Naive <- merge(matureSP[,grep('Naive', Idents(matureSP))], 
                       T_NK_2[, Idents(T_NK_2) %in% c( 'CD8+TN','CD8+RTE', 'CD4+TN','CD4+RTE')])
combined_Naive$identity <- Idents(combined_Naive)
combined_Naive$age_3 <- ''
combined_Naive$age_3[combined_Naive$age <=12] <- '<=12'
combined_Naive$age_3[combined_Naive$age >= 13 & combined_Naive$age <= 39 ] <- '13-39'
combined_Naive$age_3[combined_Naive$age >=40 & combined_Naive$age <= 99 ] <- '40-99'
combined_Naive$age_3[combined_Naive$age >=100] <- '>=100'

combined_Naive <- combined_Naive[, !combined_Naive$age_3 %in% '>=100']


CD8_combined_Naive <- combined_Naive[, grepl('CD8', Idents(combined_Naive)) & 
                                           !Idents(combined_Naive) %in% c('CD8+SOX4-Naive')]
CD8_combined_Naive$identity <- Idents(CD8_combined_Naive)

CD8_combined_Naive <- NormalizeData(CD8_combined_Naive)
CD8_combined_Naive <- FindVariableFeatures(CD8_combined_Naive)


CD8_combined_Naive <- ScaleData(CD8_combined_Naive,verbose = FALSE)#, vars.to.regress = c('batch'))
CD8_combined_Naive <- RunPCA(CD8_combined_Naive, npcs = 50, verbose = FALSE)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD8_combined_Naive, reduction = 'pca',ndims = 50)


options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_combined_Naive, reduction = 'pca', group.by= 'identity', dims = c(1:2))+
scale_color_manual(values = cluster_colors)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_combined_Naive, reduction = 'pca', group.by= 'batch', dims = c(1:2))

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_combined_Naive, reduction = 'pca', group.by= 'identity', dims = c(2:3))+
scale_color_manual(values = cluster_colors)
options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_combined_Naive, reduction = 'pca', group.by= 'batch', dims = c(2:3))

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(CD8_combined_Naive, reduction = 'pca', dims = c(2:3))+
scale_color_manual(values = cluster_colors)+
theme(panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.position = c(0.54, 0.1)
     )+
guides(color = guide_legend(override.aes = list(size = 7)))



CD4_combined_Naive <- combined_Naive[, grepl('CD4', Idents(combined_Naive)) & 
                                           !Idents(combined_Naive) %in% c('CD4+SOX4-Naive')]
CD4_combined_Naive$identity <- Idents(CD4_combined_Naive)

CD4_combined_Naive <- NormalizeData(CD4_combined_Naive)
CD4_combined_Naive <- FindVariableFeatures(CD4_combined_Naive)


CD4_combined_Naive <- ScaleData(CD4_combined_Naive,verbose = FALSE)#, vars.to.regress = c('batch'))
CD4_combined_Naive <- RunPCA(CD4_combined_Naive, npcs = 50, verbose = FALSE)


options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD4_combined_Naive, reduction = 'pca',ndims = 50)


options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD4_combined_Naive, reduction = 'pca',  dims = c(1:2))+
scale_color_manual(values = cluster_colors)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD4_combined_Naive, reduction = 'pca', group.by= 'batch', dims = c(1:2))

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD4_combined_Naive, reduction = 'pca',  dims = c(2:3))+
scale_color_manual(values = cluster_colors)
options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD4_combined_Naive, reduction = 'pca', group.by= 'batch', dims = c(2:3))

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(CD4_combined_Naive, reduction = 'pca',  dims = c(2:3))+
scale_color_manual(values = cluster_colors)+
theme(panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.position = c(0.54, 0.1)
     )+
guides(color = guide_legend(override.aes = list(size = 7)))


save(CD4_combined_Naive, CD8_combined_Naive, 
       file = 'combined_Naive.RData')




