library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
library(harmony)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggforce)

library(ggbreak)

library(openxlsx)

library(MySeuratWrappers)

library(patchwork)

library(reshape2)

library(aplot)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

#加载分开聚类的结果

load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/T_dev.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/Myeloid/Myeloid.RData')

load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/Bcell/Bcell.RData')
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/unconventional/unconventional.RData')


T_dev_3 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/T_dev_3.rds')

levels(T_dev_3)

ETP_DN <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN_Hu.rds')


filtered_combined_Niche_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_Niche_2.rds')


filtered_combined_Niche_2[['RNA']] <- filtered_combined_Niche_2[['symbol']]
DefaultAssay(filtered_combined_Niche_2) <- 'RNA'
filtered_combined_Niche_2[['symbol']] <- NULL


filtered_combined_Niche_2 <- filtered_combined_Niche_2[, filtered_combined_Niche_2$data_source=='Hu' & 
                                                      !Idents(filtered_combined_Niche_2) %in% 
                                                       c('Myelin', 'Mesothelial')]



ETP_DN$identity <- as.character(Idents(ETP_DN))
Myeloid_2$identity <- as.character(Idents(Myeloid_2))
T_dev_3$identity <- as.character(Idents(T_dev_3))
filtered_combined_Niche_2$identity <- as.character(Idents(filtered_combined_Niche_2))
Bcell_2$identity <- as.character(Idents(Bcell_2))
unconventional_3$identity <- as.character(Idents(unconventional_3))


identity <- c(ETP_DN$identity, Myeloid_2$identity, T_dev_3$identity, 
              filtered_combined_Niche_2$identity, Bcell_2$identity, unconventional_3$identity)
head(identity)

colors <- c("#060000","#683E25","#6162A2","#DB7E31","#E28F89","#EDC2C9","#D76483",
"#BD3460","#B8498D","#A92B50","#CF3126","#EDE6D4","#D8D15D","#B2C2FF","#496AB5","#4C6256",
"#8F9998","#52A7DD","#AFC7E3","#634490","#D23C86","#F2C43D","#B48A38",
"#706D36","#6FB0FF","#77B248","#8DC296","#0D5B5B","#27576E","#8EA5C4","#4A4391","#5A2E4F",
"#894D4F","#C2B5A2","#D4CB80","#EDF0C5","#6A8C69","#A1B3A7","#73BDC6",
"#3390BB","#2762A8","#E2DEEC","#DDC4DA","#60454A","#A86A31","#F6DC41","#CEDC7D","#555C54",
"#6C3B8A","#83C3C5","#B6D1DA","#75818D","#7C7BB4","#C19FC4","#BDB1B3",
"#543C20","#F6E7AE","#B9CE41","#57A849","#326A3D","#9FCEC8","#A1B4BB","#D8E6F3","#A8A8B2",
"#5A296B","#DECBCF","#E0BD95","#9A7C5A","#7E2E54","#639D48","#3E8946",
"#63B4A7","#A5D3E2","#7896B2","#7E7C91","#D1C2D5","#D5CAD0","#F2E1D1","#A38A36","#E9E778",
"#5D8143","#DAE8DB","#7FC2A7","#518A91","#173655","#26306E","#796185",
"#E6AAC6","#FBF2E9","#B2AEA3","#363F22","#A6C869","#96BB99","#003B2F","#FEFDFE","#8AB9C4",
"#BABDDC","#413549","#7C3B63","#F2CFB3","#E6BC50","#C1C53C","#B3D193","#3B6138","#9DCCBA",
"#399BA8","#455562","#5F698C","#7E9CCE","#ECE2ED")

cluster_colors <- 
  c('ETP' = '#297A6B','TP1'='#C05EA5',
    'TP2'='#e6ba77','DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA',
    'DN4'='#ff927c','γδT_P'='#ee587d','Treg_diff'='#E5D2DD',#'CD8+Naive'='#53A85F',
                    'CD4+Memory'='#808040','CD8+Memory'='#F3B1A0',
                             #'CD4+Naive'='#D6E7A3',
                   # 'CD4+GZMK+T'='#57C3F3',
    'IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3','T_agonist'='#a83117','Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'αβ_Entry_2' = '#a6c050', 'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'Treg_diff_2' = '#9566ac','CD8+SOX4+Naive'='#53A85F','CD8+SOX4-Naive' = "#F2C43D", 
                              'CD4+SOX4-Naive' ="#D76483",
                             'CD4+SOX4+Naive'='#D6E7A3'  )

saveRDS(identity, file = 'thymus_identity.rds')

DefaultAssay(ETP_DN) <- 'RNA'
ETP_DN <- FindVariableFeatures(ETP_DN, verbose = F)

DefaultAssay(T_dev_3) <- 'RNA'
T_dev_3 <- FindVariableFeatures(T_dev_3, verbose = F)

DefaultAssay(Myeloid_2) <- 'RNA'
Myeloid_2 <- FindVariableFeatures(Myeloid_2, verbose = F)


DefaultAssay(filtered_combined_Niche_2) <- 'RNA'
filtered_combined_Niche_2 <- FindVariableFeatures(filtered_combined_Niche_2, verbose = F)

DefaultAssay(Bcell_2) <- 'RNA'
Bcell_2 <- FindVariableFeatures(Bcell_2, verbose = F)

DefaultAssay(unconventional_3) <- 'RNA'
unconventional_3 <- FindVariableFeatures(unconventional_3, verbose = F)



thymus_step_2 <- merge(ETP_DN, list(T_dev_3, Myeloid_2, filtered_combined_Niche_2, Bcell_2, 
                                       #earlyT_and_ILC,
                                      unconventional_3))

thymus_step_2$age_3 <- ''
thymus_step_2$age_3[thymus_step_2$age <= 12] <- '<=12'
thymus_step_2$age_3[thymus_step_2$age > 12 & thymus_step_2$age < 40] <- '13-39'
thymus_step_2$age_3[thymus_step_2$age >= 40] <- '>=40'
thymus_step_2$age_3[is.na(thymus_step_2$age) & thymus_step_2$Age_group %in% c('Prepuberal')] <- '<=12'
thymus_step_2$age_3[is.na(thymus_step_2$age) & thymus_step_2$Age_group %in% c('Adult')] <- '13-39'
thymus_step_2$age_3[is.na(thymus_step_2$age) & thymus_step_2$Age_group %in% c('Aged')] <- '>=40'

thymus_step_2$anno_1 <- ''
thymus_step_2$anno_1[grep('TEC|Ciliated|Ionocyte|Tuft', 
                          thymus_step_2$identity, ignore.case = T)] <- 'Epi'
thymus_step_2$anno_1[grep('Endo', 
                          thymus_step_2$identity, ignore.case = T)] <- 'Endo'
thymus_step_2$anno_1[grep('Fb', thymus_step_2$identity)] <- 'Mes'
thymus_step_2$anno_1[grep('VSMCs', thymus_step_2$identity)] <- 'Mes'
thymus_step_2$anno_1[thymus_step_2$anno_1 == ''] <- 'Hema'

thymus_step_2 <- thymus_step_2[rowSums(thymus_step_2$RNA@counts) > 0,]

vars_features <- unique(c(ETP_DN$RNA@var.features, T_dev_3$RNA@var.features,
                         Myeloid_2$RNA@var.features, filtered_combined_Niche_2$RNA@var.features,
                         Bcell_2$RNA@var.features,unconventional_3$RNA@var.features))

# thymus_step_2 <- NormalizeData(thymus_step_2)
#thymus_step_2 <- FindVariableFeatures(thymus_step_2, selection.method = "vst", nfeatures = 2000)
thymus_step_2$RNA@var.features <- vars_features
thymus_step_2 <- CellCycleScoring(thymus_step_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
thymus_step_2$CC.Difference <- thymus_step_2$S.Score - thymus_step_2$G2M.Score


# options(future.globals.maxSize = 100 * 1024^3)
# plan('multicore', workers = 10)
thymus_step_2 <- ScaleData(thymus_step_2, features = vars_features)#, vars.to.regress = c('S.Score', 'G2M.Score'))


# plan('sequential')

thymus_step_2 <- RunPCA(thymus_step_2, npcs = 50, features = vars_features)

thymus_step_2$batch_2[is.na(thymus_step_2$batch_2)] <- 'sixth'

thymus_step_2 <- RunHarmony(thymus_step_2, 'batch_2')

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(thymus_step_2, ndims = 50, reduction = 'harmony')

thymus_step_2 <- RunUMAP(thymus_step_2, reduction = "harmony", dims = 1:20, verbose = F)

thymus_step_2$anno_2 <- thymus_step_2$identity
thymus_step_2$anno_2[thymus_step_2$anno_1 %in% 'Hema'] <- 'Hema'

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(thymus_step_2,  label = F,pt.size = 0.5, group.by = 'anno_1')+
scale_color_manual(values = colors)+
scale_color_manual(values = c("Endo"="#ec6960","Epi"="#ae80dc",
                              'Mes'="#72bd90","Hema"="#9cbfe7"), 
                  limits = c('Hema', 'Epi', 'Mes', 'Endo'))+
theme(legend.position = c(0.75,0.375),panel.background=element_rect(colour = "black"),
      axis.title=element_blank(),axis.ticks=element_blank(),
      #axis.text=element_blank(),
      axis.line=element_blank(), 
      legend.text = element_text(size=20, family= 'sans'))+
guides(colour = guide_legend(override.aes = list(size=7)))



saveRDS(thymus_step_2, file = 'thymus_step_2.rds')





options(repr.plot.width=16, repr.plot.height=9)
p <- FeaturePlot(thymus_step_2, cols = c('lightgrey', 'red'), pt.size = 0.1,
                 features = c('PTPRC', 'PAX9', 'COL1A1', 'CDH5', 'RAC2', 'EPCAM', 'PDGFRA', 'CLDN5'), 
                 ncol = 4)

p&theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.title=element_text(size = 20, face = 'italic', family = 'sans'),
        legend.text=element_text(size =20, family = 'sans'),
      axis.title=element_blank(),axis.ticks=element_blank(),
      axis.text=element_blank(),axis.line=element_blank(), 
       legend.position = c(0.03, 0.9), 
       legend.direction = 'horizontal')





Hema <- thymus_step_2[, thymus_step_2$anno_1 == 'Hema']

Hema$RNA@var.features <- unique(c(ETP_DN$RNA@var.features, T_dev_3$RNA@var.features,
                         Myeloid_2$RNA@var.features, #filtered_combined_Niche_2$RNA@var.features,
                         Bcell_2$RNA@var.features,unconventional_3$RNA@var.features))
# Hema <- FindVariableFeatures(Hema, verbose = F)




plan('multicore', workers = 10)
options(future.globals.maxSize = 100 * 1024^3)
Hema <- ScaleData(Hema,features = Hema$RNA@var.features)#, vars.to.regress = c('S.Score', 'G2M.Score'))
plan('sequential')


Hema <- RunPCA(Hema, npcs = 50, features = Hema$RNA@var.features)

Hema <- RunHarmony(Hema, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Hema, ndims = 50, reduction = 'harmony')


Bcell<-c('MemoryB','NaiveB','Plasma')
ETP_DN<-c('ETP','DN_Q1', 'DN_Q2','DN4','DN_P','γδT_P','TP1','TP2')
DP<-c('DP_P','DP_Q','HSP_DP','αβ_Entry')
SP<-c('CD4+Naive','CD8+Naive',
      'CD4+Memory_like','CD8+Memory_like', 'Agonist_1',
      'IFN_response','Treg','CD4+pre_T','CD8+pre_T',
      'T_agonist',
     'pre_Treg1','pre_Treg2','Agonist_2')
Innate_T <- c('γδT_3','ZNF683+CD8ααT','GNG4+CD8ααT','γδT_2','IFN_CD8ααT','γδT_1', 'NKT_like_CD8ααT')
Unconventional_cells<- c('MatureNK','MemoryNK','Lti','ImmatureNK')
Myeloid<-c('DC2','aDC','pDC','Macrophage','DC1P','DC1','Monocyte','Neutrophil','Mast')

Hema$hema_class <- ''
Hema$hema_class[Hema$identity %in% Bcell] <- 'Bcell'
Hema$hema_class[Hema$identity %in% ETP_DN] <- 'ETP_DN'
Hema$hema_class[Hema$identity %in% DP] <- 'DP'
Hema$hema_class[Hema$identity %in% SP] <- 'SP'
Hema$hema_class[Hema$identity %in% Unconventional_cells] <- 'Unconventional_cells'
Hema$hema_class[Hema$identity %in% Innate_T] <- 'Innate_T'
Hema$hema_class[Hema$identity %in% Myeloid] <- 'Myeloid'

Hema_colors <- c('Bcell'='#877630','DP'='#d37485','ETP_DN'='#CA9B5D','Innate_T'='#3ca462','Myeloid'='#60a5de','SP'='#6d69ad','Unconventional_cells'='#eb7814')



Hema <- RunUMAP(Hema, reduction = "harmony", dims = 1:30,verbose = F)
#Hema <- RunTSNE(Hema, reduction = "pca", dims = 1:30)





options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Hema,group.by = "hema_class", pt.size=0.5) +
scale_colour_manual(values=Hema_colors, label = function(x) gsub('Unconventional_cells', 'UCs', x),
                    limits = c('ETP_DN', 'DP', 'SP', 'Innate_T', 'Unconventional_cells', 'Bcell', 'Myeloid'))+
d_theme_w(size = 13.75)+
theme(legend.position = c(0.4,0.1),panel.background=element_rect(colour = "black"),
      #legend.spacing = unit(0,'cm'),
      #axis.title=element_blank(),
      axis.ticks=element_blank(),
      legend.text = element_text(size=13.75, family='sans'),axis.line=element_blank())+ 
guides(colour = guide_legend(nrow = 3,override.aes = list(size=7), byrow = T))



saveRDS(Hema, file = 'Hema.rds')




