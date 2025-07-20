library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
library(harmony)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')


load('/data1/02.private/dengyj/analysis/thymus/PB/Bcell/PB_Bcell.RData')
T_NK_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_2.rds')
load('/data1/02.private/dengyj/analysis/thymus/PB/Myeloid/Myeloid.RData')
load('/data1/02.private/dengyj/analysis/thymus/PB/Ery_MK/Ery_MK.RData')

###remove RBC which mostly come from 1 donor
PB_step_2 <- merge(T_NK_2, list(Bcell_3, Myeloid_2, 
                                Ery_MK_2[, ! Idents(Ery_MK_2) %in% c('RBC')]
                               ))

saveRDS(PB_step_2@meta.data[, c('nCount_RNA','nFeature_RNA')], file = 'PB_info_df.rds')

num_of_sample <- table(PB_step_2$age_3)
num_of_sample
save(num_of_sample, file = 'num_of_sample.rds')

PB_step_2$identity <- as.character(Idents(PB_step_2))
saveRDS(PB_step_2$identity, file = 'PB_identity.rds')

PB_step_2 <- NormalizeData(PB_step_2)
PB_step_2 <- FindVariableFeatures(PB_step_2, selection.method = 'vst', nfeatures = 2000)
plan('multisession', workers = 16)
options(future.globals.maxSize = 10 * 1024^3)
PB_step_2 <- ScaleData(PB_step_2, vars.to.regress = c('S.Score', 'G2M.Score'))
plan('sequential')
PB_step_2 <- RunPCA(PB_step_2, npcs = 50) 

PB_step_2 <- RunHarmony(PB_step_2, 'batch') 

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(PB_step_2, ndims = 50, reduction = 'harmony')

PB_step_2 <- RunUMAP(PB_step_2, reduction = "harmony", dims = 1:30)

Tcell_idents <- c('CD4+RTE', 'CD4+TN', 'CD4+TCM', 'CD8+RTE','CD8+TN', 'CD8+TCM',
           'CD4+TEFF', 'Treg', 'CD4+TEM', 'CD8+TEM-K', 
           'CD4+CTL', 'CD8+TEM-B',
           'MAIT', 'Vδ2Vγ9T', 'CD27+Vδ1T',
           'CD27lo_Vδ1T', 'Tex', 'IFN_T','T_un')
NK_idents <- c('CD56_Dim_NK','CD56_Bright_NK', 'LAG3+NK')
Idents(PB_step_2) <- PB_step_2$identity
PB_step_2$anno_1 <- PB_step_2$identity
PB_step_2$anno_1[Idents(PB_step_2) %in% Tcell_idents] <- 'Tcell'
PB_step_2$anno_1[Idents(PB_step_2) %in% NK_idents] <- 'NK'
PB_step_2$anno_1[grep('cDC', PB_step_2$anno_1)] <- 'cDC'

options(repr.plot.width=7, repr.plot.height=7)
Idents(PB_step_2) <- PB_step_2$anno_1
levels(PB_step_2) <- c('HSPC','Naive_B','Memory_B','ABC','Plasma','cDC',
                       'pDC','Mono1','Mono2', 'Mono3','MK','Tcell','NK')
p <- DimPlot(PB_step_2,reduction = "umap",pt.size=0.5)+
scale_color_manual(values = c('NK'='#e0d79e', 'Tcell'='#a3dccb', 'Plasma'= '#847525',
                             'Memory_B'='#809f8c','Naive_B'='#73BDC6','ABC'='#C2B5A2',
                               'Mono2'='#dec5ab','HSPC'='#ffb16b','MK'='#ed9594', 'Mono3' = '#cba266',
'Mono1'='#7c6368','cDC'='#B2C2FF','pDC'='#2762A8'))



new_label <- as.character(1:nlevels(PB_step_2))
names(new_label) <- levels(PB_step_2)

options(repr.plot.width=8, repr.plot.height=9.15)
add_label(p, use_new_label = T, new_label = c(new_label),add_circle_to_label = F,
          location_adj= list('4'=c(0.03, -0.01), #'1'= c(0, 0),
                             '7'=c(0, 0.02)
                             ),
         label_arg = list(color='black', size=6))+
  guides(color=guide_legend(override.aes = list(size=3.3, label = c(new_label)), ncol = 4))+
theme(axis.title = element_blank(),axis.ticks = element_blank(), 
     axis.text = element_blank(), axis.line  = element_blank(), 
     panel.border = element_rect(color = 'black', fill = NA), legend.key.width = unit(2, "lines"),
   legend.key.height = unit(2, "lines"),
      legend.text = element_text(size=20), legend.position = 'bottom', legend.justification='center')



saveRDS(PB_step_2, file = 'PB_step_2.rds')
