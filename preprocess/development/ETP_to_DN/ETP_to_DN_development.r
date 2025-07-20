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


library(harmony)
#library(monocle3)




library(ggpubr)


library(openxlsx)

library(ggraph)
library(tidygraph)


library(circlize)

library(ComplexHeatmap, lib.loc = .libPaths()[1])



library(SeuratWrappers)

library(monocle3)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')
source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

library(patchwork)


cluster_colors <- 
c(
    'CD8+Naive'='#53A85F',
                             'CD4+Naive'='#D6E7A3',
'ETP' = '#297A6B','TP1'='#C05EA5',#'ETP'='#a7aad5',
    'TP2'='#e6ba77','DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA',
    'DN4'='#ff927c','γδT_P'='#ee587d','Treg_diff'='#E5D2DD',#'CD8+Naive'='#53A85F',
                    'CD4+Memory'='#808040','CD8+Memory'='#F3B1A0',
                             #'CD4+Naive'='#D6E7A3',
                    'CD4+GZMK+T'='#57C3F3','IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3','T_agonist'='#a83117','Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'αβ_Entry_2' = '#a6c050', 'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'Treg_diff_2' = '#9566ac','CD8+SOX4+Naive'='#53A85F','CD8+SOX4-Naive' = "#F2C43D", 
                              'CD4+SOX4-Naive' ="#D76483",
                             'CD4+SOX4+Naive'='#D6E7A3',
    'CD8ααT' = '#e1812b',
   'MemoryB' = '#e0b139', 
    'NaiveB' = '#698c68', 
    'Plasma'='#ff7765',#'#d87ab3', 
   'MemoryNK'='#7ed670','NKT_like_CD8ααT'='#786085',#'#faae31',
    'MatureNK'= '#b06fa2',#'#d76482',
    'ImmatureNK'='#5065ad',
    'Lti'='#F0E442',#'#71a4c9',
    'rILCP'='#944575',#'#e6bc4f',
    'ZNF683+CD8ααT'='#de243a',
    'GNG4+CD8ααT'='#926c9c',
    'γδT_1'='#009E73','γδT_2'='#c39c62','γδT_3'='#9ba207',#'#c0c53c',
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
'PB_CD4+TN'='#c1d4eb','PB_CD8+TN'='#DDC4DA',
'PB_CD4+RTE'='#52a298'
   )

ETP_DN <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN_Hu.rds')

ETP_DN_dev <- SubsetData(ETP_DN, max.cells.per.ident = 2000)

table(Idents(ETP_DN_dev))

DefaultAssay(ETP_DN_dev) <- 'integrated'

#ETP_DN_dev <- FindVariableFeatures(ETP_DN_dev, selection.method = "vst", nfeatures = 2000, verbose = F)
#ETP_DN_dev <- CellCycleScoring(ETP_DN_dev, s.features = s.genes, g2m.features = g2m.genes)
plan('multicore', workers = 10)
options(future.globals.maxSize = 100 * 1024^3)
ETP_DN_dev <- ScaleData(ETP_DN_dev, vars.to.regress = c('CC.Difference'))
plan('sequential')
ETP_DN_dev <- RunPCA(ETP_DN_dev, npcs = 50, verbose = FALSE)

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(ETP_DN_dev, ndims = 50, reduction = 'pca')

ETP_DN_dev <- RunUMAP(ETP_DN_dev, reduction = "pca", dims = 1:25, verbose  = F)


levels(ETP_DN_dev) <- c('ETP','TP1','TP2','DN_P','DN_Q1','DN_Q2','DN4','γδT_P')


ETP_DN_dev_cds <- as.cell_data_set(ETP_DN_dev)

ETP_DN_dev_cds@clusters$UMAP$clusters <- Idents(ETP_DN_dev)[rownames(colData(ETP_DN_dev_cds))]

#parition的数目根据自己的实际情况处理，如果后面的 learn_graph(ETP_DN_dev_cds, use_partition = F)的use_partition为F，则不用partition来做不同pseudotime推测，没有没有必要设置partition；

ETP_DN_dev_cds@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(ETP_DN_dev_cds)))), levels = 1)

names(ETP_DN_dev_cds@clusters$UMAP$partitions) <- rownames(colData(ETP_DN_dev_cds))

ETP_DN_dev_cds <- estimate_size_factors(ETP_DN_dev_cds)

#创建cds对象的时候，cds要求gene_meta这gedf必须要有一列为gene_short_name

## Add gene names into CDS

ETP_DN_dev_cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(ETP_DN_dev_cds)

#如果出现了order_cells(cds)报错"object 'V1' not found"，否则不执行

#来源：https://github.com/cole-trapnell-lab/monocle3/issues/130](https://github.com/cole-trapnell-lab/monocle3/issues/130

###############################################################
rownames(ETP_DN_dev_cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@reducedDims$UMAP) <- NULL
colnames(ETP_DN_dev_cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
###############################################################
ETP_DN_dev_cds <- learn_graph(ETP_DN_dev_cds, #use_partition = F, 
                    close_loop = T)


options(repr.plot.width = 5.5, repr.plot.height = 5.5)
p <- plot_cells(ETP_DN_dev_cds, label_groups_by_cluster=FALSE,  label_leaves = T, label_roots = T,
                label_branch_points = F,      
                cell_size = 1,#label_groups_by_cluster = F,
            group_label_size=0,             
          graph_label_size = 3)+
d_theme_w(size = 13.75)+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
      legend.background = element_blank(),
      legend.position = c(0.3, 0.77),
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 
      axis.text=element_blank(),
      legend.text = element_text(size = 13.75, margin = margin(l = -2, unit = 'mm')),
      legend.title =element_text(size = 13.75,margin = margin(b = -2, unit = 'mm')))+
guides(color = guide_legend(title = 'Identity', ncol = 2, byrow = F,override.aes = list(size = 7)))+
scale_color_manual(values = cluster_colors,label = thymocyte_lbls_func)+
labs(x = 'UMAP_1',  y = 'UMAP_2')
#p

p

ETP_DN_dev_cds <- order_cells(ETP_DN_dev_cds,root_cells = 'TCATGGAAGGGTGGGA_Adult-Thy1-1118')

options(repr.plot.width = 5.5, repr.plot.height = 5.5)
p <- plot_cells(ETP_DN_dev_cds, label_groups_by_cluster=T,  color_cells_by = "pseudotime",label_leaves = T, 
           label_roots = T,
                label_branch_points = F,      
           group_label_size=8, cell_size = 0.5,
          graph_label_size = 3,label_cell_groups=F)+
d_theme_w(size = 13.75)+
#guides(color = guide_legend(override.aes = list(size = 7)))+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
      legend.position = c(0.275,0.75),
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 

      axis.text=element_blank())+
guides(color = guide_colorbar(title = 'Pseudotime', override.aes = list(size = 7)))+
labs(x = 'UMAP_1',  y = 'UMAP_2')
p




saveRDS(ETP_DN_dev, file = 'ETP_DN_dev.rds')

saveRDS(ETP_DN_dev_cds, file = 'ETP_DN_dev_cds.rds')


