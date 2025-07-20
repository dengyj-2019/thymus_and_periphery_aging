library(ggplot2)
library(Seurat)
library(MySeuratWrappers)


library(pheatmap)

library(ggfun)

library(data.table)

library(ggpubr)

library(openxlsx)

library(patchwork)

library(STEMNET)

library(ComplexHeatmap, lib.loc = .libPaths()[1])
library(aplot)

library(circlize)

library(plyr)
library(reshape2)

library(future)

library(ggforce)

library(rstatix)

library(ggalluvial)

library(ggraph)
library(tidygraph)


library(monocle3)

plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig2_FigS2'

setwd(plot_dir)

cluster_colors <- 
c(
    'CD8+Naive'='#53A85F','ETP' = '#297A6B',
     'TP1'='#C05EA5','TP2'='#e6ba77',#'ETP'='#a7aad5',
    'TP'='#e6ba77',
                             'CD4+Naive'='#D6E7A3','ETP-TP'='#efd703',#'#a7aad5',
'DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA',
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
'other' = 'lightgrey',
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

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

breaks_func <- function(x){
    seq(min(x), max(x), length.out = 3)
}







#Fig.2a
ETP_DN <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN_Hu.rds')
plot_df <- cbind(FetchData(ETP_DN, vars = c('combined_idents')), 
                 Embeddings(ETP_DN, reduction = 'umap'))



options(repr.plot.width=4.3, repr.plot.height=4.3)



p <- ggplot(plot_df, aes(UMAP_2, UMAP_1, color=combined_idents))+
geom_point(size = 1)+
theme(legend.position = c(0.22,0.25),  legend.background=element_blank(),    
      panel.background=element_rect(fill = 'white', color=NA),
      panel.border=element_rect(fill = NA, color='black'),
    axis.title=element_text(size = 13.75),axis.ticks=element_blank(),
      legend.title=element_blank(),
    axis.text=element_blank(),axis.line=element_blank(), 
    legend.text = element_text(size = 13.75))+

guides(colour = guide_legend(override.aes = list(size=7), ncol = 1)) +
scale_color_manual(values = cluster_colors,label = function(x){
    x <- gsub('_', '-', x)
    #x <- gsub('ETP-TP', 'T progenitors',x)
    x <- gsub('DN4', 'ISP-like',x)
})

p
my_saveplot_fun(p, 'ETP_DN_population', width=4.4, height=4.4)

###Extended Data Fig.2a
gene_order<-c('CD34','NFE2', 'BTK', 'SPINK2', 'HHEX', 'BAALC', 'IGLL1', 'IFITM1','TRGC1',
                  'CCNB2','CDC20','CDC6','CENPA',
                  'PCNA','MKI67','CDK1',#'HES1',
              #'HES4','MME','HIVEP3','PTCRA','MAL',
              'RAG1','RAG2','BCL11B',
'TRBC1','TRBC2','TRBV7-2','TCF12','IL7R','DNTT','NOTCH3','TSHR',#'CXCR4',
              'VIPR2','NINL',
              'EPHX2','NLGN4X','SIT1','DDIT4','TRAC','CD4','CD8B',#'RORC',
              #'AQP3',
              'CD8A',
              'CD38','CD1A','TCF7','CD2','LEF1','AEBP1','LCK','ZAP70','GATA3','CCR9','ITGB2',#'ID2',
                  'EGR1','ID3',
              'TRGC2','ICAM2',
              'TRDC',
'TRDV1','TRGV4')

DefaultAssay(ETP_DN) <- 'RNA'
p <- DotPlot(ETP_DN, features = rev(gene_order), group.by = 'combined_idents')

p1 <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_y_discrete(limits =rev, position = 'right', 
                label = function(x){
    x <- gsub('_', '-', x)
    #x <- gsub('ETP-TP', 'T progenitors',x)
    x <- gsub('DN4', 'ISP-like',x)
})+ scale_size(range=c(0, 8))+
scale_x_discrete(position = "bottom") +
theme(axis.title = element_blank(), plot.margin = margin(l=0.5,unit = 'cm'),
      legend.position = 'right',legend.box  ='horizontal',
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.y = element_text(hjust = 1, #vjust = 0.5,
                                 size = 13.75, family = 'sans', face = 'plain'),
     axis.text.x = element_text(hjust =1, angle =45,size = 13.75, family = 'sans', face = 'italic', vjust = 1), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression")) 

options(repr.plot.width=17, repr.plot.height=4)
p1
ggsave('DN_DEGs.svg', width = 17, height = 4)



###Extended Data Fig.2b
ETP_DN_dev_cds <- readRDS('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/ETP_DN/ETP_DN_dev_cds.rds')

ETP_DN_dev <- readRDS('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/ETP_DN/ETP_DN_dev.rds')

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
my_saveplot_fun(p, file = 'ETP_DN_dev_idents', width = 5.5, height = 5.5)

###Extended Data Fig.2c
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
my_saveplot_fun(p, file = 'ETP_DN_dev_pseudotime', width = 5.5, height = 5.5)





###Extended Data Fig.2d
T_dev_3 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/T_dev_3.rds')

ETP_DP <- merge(ETP_DN, T_dev_3[, grepl('DP', Idents(T_dev_3))])

ETP_DP$id <- as.character(Idents(ETP_DP))

ETP_DP$id[grepl('DP', ETP_DP$id)] <- 'DP'
ETP_DP$id[grepl('TP', ETP_DP$id)] <- 'ETP-TP'

DefaultAssay(ETP_DP) <- 'RNA'

saveRDS(ETP_DP, file = 'ETP_DP.rds')

p <- DotPlot(ETP_DP, features = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD40LG', 'CD8A', 'CD8B'), group.by = 'id')

p1 <- ggplot(p$data, aes( id, features.plot,fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_x_discrete(limits =c(   'ETP-TP', 'DN_P','DN_Q1','DN_Q2','DN4','γδT_P','DP'), label = function(x){
    x <- gsub('_', '-', x)
    #x <- gsub('ETP-TP', 'T progenitors',x)
    x <- gsub('DN4', 'ISP-like',x)
})+ scale_size(range=c(3, 8))+
scale_y_discrete(position = "right") +
theme(axis.title = element_blank(), legend.key.height = unit(0.5,'cm'),
      plot.margin = margin(l = 1,unit = 'cm'),
      legend.margin = margin(l = -0.2,unit = 'cm'),
      legend.position = 'right',#legend.box  ='vertical',
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.x = element_text(hjust = 1, angle=45,
                                 size = 13.75, family = 'sans', face = 'plain'),
     axis.text.y = element_text(hjust =1, angle =0,size = 13.75, family = 'sans', face = 'italic', vjust = 0.5), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", ncol = 1,override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression"))

options(repr.plot.width=5.5, repr.plot.height=5.5)
p1
ggsave(file = 'DN_DP_DEGs.svg', width = 5.5, height = 5.5)

#Fig.2b
Hema <- readRDS('/data1/02.private/dengyj/analysis//thymus/clustering_10_20/combined_result/Hema.rds')

Hema_colors <- c('Bcell'='#877630','DP'='#d37485','ETP_DN'='#CA9B5D','Innate_T'='#3ca462','Myeloid'='#60a5de',
                 'SP'='#6d69ad','Unconventional_cells'='#f6df18')#'#eb7814')

Hema$Age_group[is.na(Hema$Age_group)] <- Hema$age_3[is.na(Hema$Age_group)]
Hema$Age_group[Hema$Age_group == '<=12'] <- 'Prepuberal'
Hema$Age_group[Hema$Age_group == '13-39'] <- 'Adult'
Hema$Age_group[Hema$Age_group == '>=40'] <- 'Aged'

IFN_response <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/IFN_response.rds')


ETP_DN_idents <-c('ETP','DN_Q1', 'DN_Q2','DN4','DN_P','γδT_P','TP2','TP1')


Hema$identity2 <- Hema$identity
Hema$identity2[colnames(IFN_response)[Idents(IFN_response) == 'CD4+IFN_response']] <- 'CD4+IFN_response'
Hema$identity2[colnames(IFN_response)[Idents(IFN_response) == 'CD8+IFN_response']] <- 'CD8+IFN_response'
Hema$hema_class2 <- 'other'



Hema$hema_class2[Hema$identity2 %in% ETP_DN_idents] <- 'ETP-DN'

logi <- !Hema$batch_2 %in% c('fifth', 'sixth')
Hema$identity3 <- Hema$identity2
Hema$identity3[grep('ETP|TP1|TP2',Hema$identity3)] <- 'ETP-TP'
plot_df <- 
data.frame(apply(table(Hema$identity3[logi], Hema$batch[logi]), 2, function(x){x/sum(x)}), check.names = F)
plot_df$cluster <- rownames(plot_df)

plot_df <- melt(plot_df, id.vars = 'cluster', value.name = 'pct', variable.name = 'batch')
plot_df$age_3 <- Hema$age_3[match(plot_df$batch, Hema$batch)]
plot_df$age <- Hema$age[match(plot_df$batch, Hema$batch)]
plot_df$pct <- plot_df$pct*100
plot_df <- plot_df[plot_df$cluster %in% Hema$identity3[Hema$hema_class2 %in% c('ETP-DN')], ]

plot_df <- ddply(plot_df, `.variables` = c('cluster', 'age_3'), summarise, pct = mean(pct))
plot_df$cluster <- factor(plot_df$cluster, 
                          levels = c('ETP-TP', 'DN_P', 
                                     'DN_Q1', 'DN_Q2', 'DN4','γδT_P'))


options(repr.plot.width=4, repr.plot.height=4.3)
ggplot(plot_df, aes(age_3, pct, fill  =cluster))+
geom_bar(stat = 'identity', width = 0.5)+
#geom_col(width= 0.5,color = "white",size = 0.1)+

d_theme_w(size = 13.75)+
labs(x = NULL, y = 'of CD45+ %', fill = 'Identity')+
scale_x_discrete(label = age_group_lbls_func, expand = c(0, 0))+
scale_fill_manual(values = cluster_colors, label = function(x){
    x <- gsub('_', '-', x)
    #x <- gsub('ETP-TP', 'T progenitors',x)
})+
scale_y_continuous(expand = c(0,0.1))+
theme(legend.margin = margin(l = 0,unit = 'cm'),
           panel.border = element_blank(),#element_rect(color = 'black', fill = NA), 
      panel.background = element_blank())+#element_rect(color = NA, fill = 'white'),)+
  geom_flow(aes(alluvium = cluster), alpha= .3, color = "white",
            curve_type = "linear", 
            width = .5)
ggsave('ETP_DN_proportion.svg', width=4, height = 4.3)



###Extended Data Fig.2e
get_Age_test_T <- function(
  data
){
  comparison_group <- colnames(data)[grep('_comparison', colnames(data))]
  comparison_group <- unique(as.character(as.matrix(data[, comparison_group, drop=F])))
  data_list <- lapply(comparison_group, function(comparison_group_iter){
    data_iter <- data[, grep(comparison_group_iter, colnames(data), fixed= T), drop  =F]
    colnames(data_iter) <- gsub(paste0(comparison_group_iter, '_'), '', colnames(data_iter))
    data_iter
  })
  data_list <- do.call(rbind, data_list)
  data_list
}

Age_Gene_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/thymus/ETP_DN//Age_Gene_test_list.rds')

DN_Q1_df <- get_Age_test_T(Age_Gene_test_list[['DN_Q1']])

DN_Q2_df <- get_Age_test_T(Age_Gene_test_list[['DN_Q2']])



DN_Q1_selected_genes <- c(#'HES4',
    'TRBC1','TRBV7-2', 'TRBV6-6',
                          'ZAP70','PTCRA','CD3E',
    #'CD7','ID3',
    'LAT','RAG2','IL7R','LEF1','TCF7','IKZF1',
'TOX','RAG1','NFRKB','CD2','LYL1','JUN','FOSL2','NR4A2','EGR1','DUSP1','NR4A1',
'FOS','TOX2','JUNB','FOSB','IGKC','IGHA1','JUN','IGHG1','IGLC2','IGLC3','S100A6',
'CD52','ANXA1','CDKN2A','MME','ICOS','IGLV1-47','HLA-DQA1','KLRB1','CDKN1A','MEIS2',
    'TRAV22'

)

DN_Q1_plot_df <- DN_Q1_df[DN_Q1_df$variable %in% DN_Q1_selected_genes, ]

DN_Q1_plot_df <- plyr::ddply(DN_Q1_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DN_Q1_plot_df$comparison <- gsub('_VS_Rest', '', DN_Q1_plot_df$comparison)

DN_Q1_label_df <- data.frame(Age = unique(DN_Q1_plot_df$comparison))



DN_Q1_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DN_Q1_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



DN_Q1_label_plot <- ggplot(DN_Q1_label_df, aes(1,Age, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
scale_y_discrete(limits = rev)+
labs(fill = 'Age', title = NULL)+
guides(fill = guide_legend(nrow = 3, byrow = T))




DN_Q2_selected_genes <- c('TRBC1','IGKC','IGLC2','IGHG4','IGHG3','IGHG1','IGHA1','TCF7','FOS','TRBV7-2',
'CDKN2A','CD7','ANXA1','IKZF1','JUN','FOSB','DUSP1','DUSP1',
'NR4A2','RAG1','IGLV1-44','IL7R','S100A6','IGKV3-20','FOS',#'JUND',
                          'IGKV3-15','JUNB','CD3E', 'LCK', 'LEF1'
                        
)

DN_Q2_plot_df <- DN_Q2_df[DN_Q2_df$variable %in% DN_Q2_selected_genes, ]

DN_Q2_plot_df <- plyr::ddply(DN_Q2_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DN_Q2_plot_df$comparison <- gsub('_VS_Rest', '', DN_Q2_plot_df$comparison)

DN_Q2_label_df <- data.frame(Age = unique(DN_Q2_plot_df$comparison))



DN_Q2_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DN_Q2_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


DN_Q2_label_plot <- ggplot(DN_Q2_label_df, aes(1,Age, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label =  age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
scale_y_discrete(limits = rev)+
labs(fill = 'Age', title = NULL)+
guides(fill = guide_legend(nrow = 3, byrow = T))



combined_p_val_adj <- c(DN_Q1_plot_df$p_val_adj, DN_Q2_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
DN_Q1_plot_df$p_val_adj[DN_Q1_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
DN_Q2_plot_df$p_val_adj[DN_Q2_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- round(seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5), digits = 0)
#signif(seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5), digits = 1)

combined_exp <- c(DN_Q1_plot_df$scale_exp, DN_Q2_plot_df$scale_exp)

combined_exp_break <- round(seq(min(combined_exp), max(combined_exp), length.out = 5), 1)

DN_Q1_p <- ggplot(DN_Q1_plot_df, aes(variable,comparison, # size = -log10(p_val_adj), 
                                     fill = scale_exp))+
geom_tile(color = 'lightgrey')+
# scale_size_continuous(range = c(0.5,7), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f6be53'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_x_discrete(limits=DN_Q1_genes_order, position = 'bottom')+
scale_y_discrete(limits = rev(c('<=12', '13-39', '>=40')), position = 'right',
                 label = age_group_lbls_func)+
guides(fill = guide_colorbar(title = "Expression", order = 0,barwidth = unit(1, 'cm'))#,
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 1)
      )+
labs(x='Age', y=NULL, title = 'DN-Q1')+
theme(plot.title = element_text(size = 17.5, hjust=0.5, face = 'bold'),
      axis.text.y = element_text(size = 13.75, family = 'sans', angle = 0, hjust = 1,  vjust = 1), 
      #axis.text.x = element_blank(),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_text(size = 13.75, family = 'sans', face = 'italic', angle = 60, vjust = 1, hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans', vjust = 1),


      axis.line =element_blank())




options(repr.plot.width=14*1/3, repr.plot.height=11)
DN_Q2_p <- ggplot(DN_Q2_plot_df, aes(variable, comparison, #size = -log10(p_val_adj), 
                                     fill = scale_exp))+
geom_tile(color = 'lightgrey')+
# scale_size_continuous(range = c(0.5,7), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f6be53'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_x_discrete(limits=DN_Q2_genes_order, position = 'bottom')+
scale_y_discrete(limits = rev(c('<=12', '13-39', '>=40')), position = 'right',
                 label = age_group_lbls_func)+
guides(fill = guide_colorbar(title = "Expression", order = 0, barwidth = unit(1, 'cm'))#,
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 1)
      )+
labs(x='Age', y=NULL, title = 'DN-Q2')+
theme(plot.title = element_text(size = 17.5, hjust=0.5, face = 'bold'),
      axis.text.y = element_text(size = 13.75, family = 'sans', angle = 0, hjust = 1, vjust = 1), 
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_text(size = 13.75, family = 'sans', face = 'italic', angle = 60, vjust = 1, hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans', vjust = 1),


      axis.line =element_blank())
options(repr.plot.width=17, repr.plot.height=5)
p <- wrap_plots(DN_Q1_label_plot,  DN_Q1_p,DN_Q2_label_plot, 
                DN_Q2_p)+
plot_layout(ncol = 2, guide = 'collect', widths  = c(0.5,19))&
theme(legend.position = 'right', plot.margin = margin(l=0.3,unit = 'cm'), 
      plot.title = element_text(size = 13.75, hjust = 0.5,family='sans'))
p
ggsave('DN_Q1_Q2_DEGs.svg', width = 17, height = 5)









###Fig.2c
ETP_TP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_TP.rds')
options(repr.plot.width=4.3, repr.plot.height=4.3)
plot_df <- cbind(Embeddings(ETP_TP, reduction = 'umap'), data.frame(identity=Idents(ETP_TP)))


p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = identity))+
geom_point(size=1)+
  theme(legend.position = c(0.78,0.89),
        panel.background=element_rect(colour = NA, fill = 'white'),
        panel.border=element_rect(colour = 'black', fill = NA),
        axis.title=element_text(size = 13.75),axis.ticks=element_blank(),
        legend.title=element_blank(),
        legend.background=element_blank(),
        axis.text=element_blank(),axis.line=element_blank(), 
        legend.text = element_text(size = 13.75))+#ylim(-6,6)+xlim(-6,7.5)+
  guides(colour = guide_legend(override.aes = list(size=7), ncol = 1)) +
  scale_color_manual(values=cluster_colors)#+
#stat_ellipse(level = 0.82, linetype='dashed', size = 2, show.legend =F)
p
my_saveplot_fun(p, 'ETP_TP_umap', width = 4.4, height = 4.4)






###Fig.2d
DefaultAssay(ETP_TP) <- 'RNA'
gene_order<-unique(c(
'MME', 'CD44', 
'HOPX',  'CCR9', 'CCR7', 'TYROBP', 'MEF2C', 
    #'IL3RA', #'BCL11A' ,
'HOXA9', 
    'KLRB1',   'FCER1G','CD33','IL1B',  #'LYL1', 
    'MPO',# 'SELL',

 'PRSS57', 'IGFBP7',#'CD1A',
    'IGLL1', 'CD7', 'CD2','BCL11B','PTCRA','RAG1', 'TCF7'
            ))

p <- DotPlot(ETP_TP, features = rev(gene_order))

p1 <- ggplot(p$data, aes(id, features.plot,fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_x_discrete(limits = 
                 levels(ETP_TP), position = 'bottom', 
                label = function(x){
    x <- gsub('_', '-', x)
})+ scale_size(range=c(0.5, 6.5))+
scale_y_discrete(position = "right", limits=rev) +
d_theme_w(size = 13.75, italic_text = '',color = 'black')+
theme(axis.title = element_blank(), plot.title  = element_text(hjust = 0.5), 
      legend.position = 'right',#legend.box  ='horizontal',
      legend.key.height = unit(0.75, 'cm'),
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.x = element_text(hjust = 0.5, angle=0,#vjust = 0.5,
                                 size = 13.75, family = 'sans', face = 'plain'),
     axis.text.y = element_text(hjust =0, angle =0,size = 13.75, family = 'sans', face = 'italic', vjust = 0.5), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'), breaks = c(-1,0,1))+
guides(size = guide_legend(title = "% Cells", title.position = 'top',
                           override.aes = list(fill = 'black'), ncol= 1), 
       fill = guide_colorbar(title = "Expression",title.position = 'top')) #+
#labs(title = 'Celltype specific DEGs of ETP-TP')
options(repr.plot.width=4, repr.plot.height=4.3)
p1
ggsave('ETP_TP_DEGs.svg', width=4, height = 4.3)












###Fig.2e
Age_Gene_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/thymus/ETP_TP/Age_Gene_test_list.rds')

get_Age_test_T <- function(
  data
){
  comparison_group <- colnames(data)[grep('_comparison', colnames(data))]
  comparison_group <- unique(as.character(as.matrix(data[, comparison_group, drop=F])))
  data_list <- lapply(comparison_group, function(comparison_group_iter){
    data_iter <- data[, grep(comparison_group_iter, colnames(data), fixed= T), drop  =F]
    colnames(data_iter) <- gsub(paste0(comparison_group_iter, '_'), '', colnames(data_iter))
    data_iter
  })
  data_list <- do.call(rbind, data_list)
  data_list
}

Age_Gene_test_list_T <- lapply(Age_Gene_test_list, get_Age_test_T)

ETP_df <- Age_Gene_test_list_T[['ETP']]


TP1_df <- Age_Gene_test_list_T[['TP1']]


ETP_selected_genes <- c('LCK','CD2','SELL','CD34','PTCRA','TRGC2','TRBC1','CD3D','CD3G','CD7','RUNX3','FOSL2',
'BCL11A','DUSP2','BHLHE40','MME','PROM1','HLA-DRB5','HLA-DRB1','SPI1',
'IGHD','IGHM','JUNB','JUND','IL3RA','IGHV1-69','LMO4','IGKC','NR4A2',#'JCHAIN',
                        'DUSP1',
'TSC22D3',
                        'FOS','MEF2A','JUND','FOSB','KLRB1','IGHA1','IGHG1')

ETP_plot_df <- ETP_df[ETP_df$variable %in% ETP_selected_genes, ]

ETP_plot_df <- plyr::ddply(ETP_plot_df, c('variable'), transform, scale_exp = scale(means_1))
ETP_plot_df$comparison <- gsub('_VS_Rest', '', ETP_plot_df$comparison)

ETP_label_df <- data.frame(Age = unique(ETP_plot_df$comparison))



ETP_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(ETP_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



ETP_label_plot <- ggplot(ETP_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
#scale_x_discrete(limits = rev)+
labs(fill = 'Age', title= 'ETP')+
guides(fill = guide_legend(nrow = 3, byrow = T))




TP1_selected_genes <- c('CD2','CD28','LEF1','TCF7','TRGC1','TRBC1','CD3D','CD3G','CD27','CD7','LMO4',
'PBX1','MEIS1','MME','HOPX','MEF2C','HOXA9','FOXO1','FOS','IGHD','IGHM',
'TCF4','CEBPB','JUNB','JUND','IL3RA','KLRB1','S100A6','IGKC','NR4A2','LMO2',
'EGR1','JCHAIN','PRDM1','MS4A1','IRF8','JUNB','FOSB','SPIB','S100A9','S100A8',
    'KLRD1','IGHA1','IGHG1'

                          
)

TP1_plot_df <- TP1_df[TP1_df$variable %in% TP1_selected_genes, ]

TP1_plot_df <- plyr::ddply(TP1_plot_df, c('variable'), transform, scale_exp = scale(means_1))
TP1_plot_df$comparison <- gsub('_VS_Rest', '', TP1_plot_df$comparison)

TP1_label_df <- data.frame(Age = unique(TP1_plot_df$comparison))



TP1_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(TP1_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


TP1_label_plot <- ggplot(TP1_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label =  age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
#scale_x_discrete(limits = rev)+
labs(fill = 'Age', title= 'TP1')+
guides(fill = guide_legend(nrow = 3, byrow = T))



combined_p_val_adj <- c(ETP_plot_df$p_val_adj, TP1_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
ETP_plot_df$p_val_adj[ETP_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
TP1_plot_df$p_val_adj[TP1_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- round(seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5), digits = 0)
#signif(seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5), digits = 1)

combined_exp <- c(ETP_plot_df$scale_exp, TP1_plot_df$scale_exp)

combined_exp_break <- round(seq(min(combined_exp), max(combined_exp), length.out = 5), 1)

ETP_p <- ggplot(ETP_plot_df, aes(comparison, variable,# size = -log10(p_val_adj),
                                 fill = scale_exp))+
geom_tile()+
# scale_size_continuous(range = c(0.5,7), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f6be53'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(ETP_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(fill = guide_colorbar(title = "Expression", order = 0,barwidth = unit(0.7, 'cm'))#,
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 1)
      )+
labs(y='Age', y=NULL)+
theme(axis.text.x = element_text(size = 13.75, family = 'sans', angle = 45, hjust = 1,  vjust = 1), 
      #axis.text.x = element_blank(),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans', face = 'italic', angle = 0, vjust = 0.5, hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans', vjust = 1),


      axis.line =element_blank())




options(repr.plot.width=14*1/3, repr.plot.height=11)
TP1_p <- ggplot(TP1_plot_df, aes(comparison, variable, #size = -log10(p_val_adj),
                                 fill = scale_exp))+
geom_tile()+
# scale_size_continuous(range = c(0.5,7), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f6be53'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(TP1_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(fill = guide_colorbar(title = "Expression", order = 0, barwidth = unit(0.7, 'cm'))#,
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 1)
      )+
labs(y='Age', y=NULL)+
theme(axis.text.x = element_text(size = 13.75, family = 'sans', angle = 45, hjust = 1, vjust = 1), 
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans', face = 'italic', angle = 0, vjust = 0.5, hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans', vjust = 1),


      axis.line =element_blank())
options(repr.plot.width=6.2, repr.plot.height=7.8)
p <- wrap_plots(ETP_label_plot,TP1_label_plot,   ETP_p,
                TP1_p)+
plot_layout(ncol = 2, guide = 'collect', heights = c(1,19))&
theme(legend.position = 'right', plot.margin = margin(l=0.3,unit = 'cm'), 
      plot.title = element_text(size = 13.75, hjust = 0.5,family='sans'))
p
ggsave('ETP_TP1_DEGs.svg', width=6.2, height = 7.8)







###Extended Data Fig.2f
ETP_STEMNET_mature <- ETP_STEMNET[, !grepl('ETP|TP', Idents(ETP_STEMNET))]

ETP_STEMNET_mature$idents <- as.character(Idents(ETP_STEMNET_mature))

ETP_STEMNET_mature$idents[ETP_STEMNET_mature$idents %in% c('DN_P')] <- 'T_lineage'
ETP_STEMNET_mature$idents[ETP_STEMNET_mature$idents %in% c('NaiveB')] <- 'B_lineage'
ETP_STEMNET_mature$idents[ETP_STEMNET_mature$idents %in% c('DC1P', 'DC1')] <- 'DC1'
ETP_STEMNET_mature$idents[ETP_STEMNET_mature$idents %in% c('Lti', 'rILCP')] <- 'ILC'
ETP_STEMNET_mature$idents[ETP_STEMNET_mature$idents %in% c('MatureNK')] <- 'NK'
Idents(ETP_STEMNET_mature) <- ETP_STEMNET_mature$idents
levels(ETP_STEMNET_mature) <- c('T_lineage', 'B_lineage', 'NK', 'ILC', 'pDC', 'Monocyte', 'DC1', 'DC2', 'Ery')

p <- VlnPlot(ETP_STEMNET_mature,stacked=T,pt.size=0,direction = 'horizontal',
             features = c( 'CD3G','CD3D', 
'MS4A1', 'PAX5', 'NCR1', 'KLRF1', 'KLRB1','IL7R','TNFSF11','KIT', 'LILRA4', 'IL3RA', 'CD14', 'VCAN' ,
                     'CLEC9A', 'XCR1', 'CD1C','CLEC10A', 
                     'GYPA', 'KLF1'),
       cols = plot_colors[levels(ETP_STEMNET_mature)], line.size = 0 )+
theme(axis.text.y = element_text(size = 13.75, family = 'sans'),
      axis.ticks.y = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_blank(),
      strip.background.x = element_blank(), 
     strip.text.x = element_text(size = 13.75, hjust = 0, face = 'italic', family = 'sans', angle =90),
      panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     plot.margin=unit(c(0, 0, 0, 1),'lines'))+
scale_x_discrete(position = 'top',label = function(x){ 
    x <- gsub('_', '-', x)
    x
}, limits = rev)+
scale_y_reverse()
options(repr.plot.width=17, repr.plot.height=4)
#p

p+facet_wrap(.~features, strip.position = 'bottom', nrow = 1)+
theme(strip.text.x = element_text(size = 13.75, hjust = 0.5, face = 'italic', family = 'sans', angle =45))
ggsave('ETP_STEMNET_mature_DEGs.svg', width =17, height = 4)





result <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Development/STEMNET_for_ETP/result_seed1_1000_cells.rds')

ETP_STEMNET <- readRDS('/data1/02.private/dengyj/analysis/thymus/Development/STEMNET_for_ETP/ETP_STEMNET.rds')
options(repr.plot.width=8, repr.plot.height=8)
plot_info <- STEMNET::plot(result)
plot_df <- plot_info$PlotData
plot_df <- cbind(plot_df, ETP_STEMNET[[]][rownames(plot_df), ], result@posteriors[rownames(plot_df), ])
plot_df$identity <- as.character(plot_df$identity)
plot_colors <- c('T_lineage' = '#007dbe', 'B_lineage' = '#8752a8', 'NK' = '#71b4e4',
                 'ILC' = '#dcd423',  'pDC' = '#dd7fba',
 'Monocyte' = '#d9b8db',
                 'DC1' = '#65c84f', 'DC2' = '#6e9a8e','Ery' = '#f8000d')



####Fig.2f

options(repr.plot.width=3.8, repr.plot.height=3.8)
p <- ggplot(plot_df, aes(x, y, shape = class, color = direction)) + geom_point(size = 0.5)+
  scale_shape_manual(values = c('FALSE' = 2, 'TRUE' = 19), labels = c('Mature Cells', 'ETP-TP'))+
  scale_color_manual(values = plot_colors, label = function(x){
      x <- gsub('_lineage', '', x)
  }) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        panel.background = element_rect(color= NA, fill ='white'),
        panel.border = element_rect(color= 'black', fill =NA),
        legend.key=element_rect(fill= 'white'),
        legend.text = element_text(size = 13.75, family = 'sans'),
        legend.position = 'bottom',legend.box='vertical',legend.box.just='left',
       legend.title = element_text(size = 13.75,family = 'sans')) + 
  labs(shape = 'Class', colour = 'Lineage') + 
guides(color = F, 
       shape = F)+
stat_ellipse(aes(group = direction), data = 
             plot_df[!plot_df$class & !plot_df$direction %in% 
                     c('T_lineage','B_lineage', 'Ery', 'NK'), ], 
             color = 'black',linetype = 2, type = 'norm', size = 0.5, level = 0.7) + 
stat_ellipse(aes(group = direction), data = 
             plot_df[!plot_df$class & plot_df$direction %in% 
                     c('T_lineage'), ], 
            color = 'black',linetype = 2, type = 'norm', size = 0.5, level = 0.9999) + 
stat_ellipse(aes(group = direction), data = 
             plot_df[!plot_df$class & plot_df$direction %in% 
                     c('B_lineage'), ], 
            color = 'black',linetype = 2, type = 'norm', size = 0.5, level = 0.99) + 
stat_ellipse(aes(group = direction), data = 
             plot_df[!plot_df$class & plot_df$direction %in% 
                     c('Ery'), ], 
            color = 'black',linetype = 2, type = 'euclid', size = 0.5, level = 0.12)+
stat_ellipse(aes(group = direction), data = 
             plot_df[!plot_df$class & plot_df$direction %in% 
                     c('NK'), ], 
            color = 'black',linetype = 2, type = 'norm', size = 0.5, level = 0.9999) 
p
my_saveplot_fun(p, 'ETP_STEMNET', width = 3.8, height = 3.8)

###Fig.2f legend
lgd1 = Legend(labels = gsub('_', '-', names(plot_colors)), title = "Lineage", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(2, "mm"),
              type = "points", pch = 19,size=unit(0.4, "cm"), 
              ncol = 1,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = plot_colors))

lgd2 = Legend(labels = c('Mature Cells', 'ETP-TP'),size=unit(0.4, "cm"),
              title_gap = unit(2, "mm"),
              title = "Class", type = "points", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              ncol= 1,title_position = 'topleft',
              legend_gp = gpar(cex = 5),
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              pch = c(2,19))

pd = packLegend(lgd2,lgd1, direction = "vertical")


options(repr.plot.width=1.3, repr.plot.height=3.2)
ggplotify::as.ggplot(grid::grid.grabExpr(
draw(pd,x = unit(0, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))))
ggsave('ETP_stemnet_lgd.svg', width = 1.3, height = 3.2)

#Fig.2g
options(repr.plot.width=6, repr.plot.height=6)
tmp_plot_df1 <- plot_df[plot_df$Age_group != 'Prepuberal'| ! plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- plot_df[plot_df$Age_group == 'Prepuberal' & plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity!='TP1',], 
                     tmp_plot_df2[tmp_plot_df2$identity=='TP1',])
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity=='TP2',], 
                     tmp_plot_df2[tmp_plot_df2$identity!='TP2',])

p1 <- ggplot(tmp_plot_df1, aes(x,y))+
  geom_point(color = 'lightgrey', size = 1 ,shape =16)+
  geom_point(data = tmp_plot_df2,
             mapping = aes(color = identity), size = 1,shape =15)+
  scale_color_manual(values = 
                       c('ETP' = '#297A6B','TP1'='#C05EA5',
                         'TP2'='#e6ba77', 'other' = 'lightgrey'), 
                     breaks = c('ETP', 'TP1', 'TP2'))+
  guides(color = F)+
theme(panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
     axis.ticks = element_blank(), 
     axis.title = element_blank(),
     axis.text = element_blank())+
annotation_custom(grid::textGrob(label = 'Group I',hjust = 0.5,
                                   x=grid::unit(0.8,"npc") ,                                
                                   y=grid::unit(0.9,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))

options(repr.plot.width=6, repr.plot.height=6)
tmp_plot_df1 <- plot_df[plot_df$Age_group != 'Adult'| ! plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- plot_df[plot_df$Age_group == 'Adult' & plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity!='TP1',], 
                     tmp_plot_df2[tmp_plot_df2$identity=='TP1',])
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity=='TP2',], 
                     tmp_plot_df2[tmp_plot_df2$identity!='TP2',])
p2 <- ggplot(tmp_plot_df1, aes(x,y))+
  geom_point(color = 'lightgrey', size = 1 ,shape =16)+
  geom_point(data = tmp_plot_df2,
             mapping = aes(color = identity), size = 1,shape =16)+
  scale_color_manual(values = 
                       c('ETP' = '#297A6B','TP1'='#C05EA5',
                         'TP2'='#e6ba77', 'other' = 'lightgrey'), 
                     breaks = c('ETP', 'TP1', 'TP2'))+
  guides(color = F)+
theme(panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
     axis.ticks = element_blank(), 
     axis.title = element_blank(),
     axis.text = element_blank())+
annotation_custom(grid::textGrob(label = 'Group II',hjust = 0.5,
                                   x=grid::unit(0.8,"npc") ,                                
                                   y=grid::unit(0.9,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))

options(repr.plot.width=6, repr.plot.height=6)
tmp_plot_df1 <- plot_df[plot_df$Age_group != 'Aged'| ! plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- plot_df[plot_df$Age_group == 'Aged' & plot_df$identity %in% c('ETP', 'TP1', 'TP2'), ]
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity!='TP1',], 
                     tmp_plot_df2[tmp_plot_df2$identity=='TP1',])
tmp_plot_df2 <- rbind(tmp_plot_df2[tmp_plot_df2$identity=='TP2',], 
                     tmp_plot_df2[tmp_plot_df2$identity!='TP2',])
p3 <- ggplot(tmp_plot_df1, aes(x,y))+
  geom_point(color = 'lightgrey', size = 1 ,shape =16)+
  geom_point(data = tmp_plot_df2,
             mapping = aes(color = identity), size = 1,shape =17)+
  scale_color_manual(values = 
                       c('ETP' = '#297A6B','TP1'='#C05EA5',
                         'TP2'='#e6ba77', 'other' = 'lightgrey'), 
                     breaks = c('ETP', 'TP1', 'TP2'))+
  guides(color = F)+
theme(panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
     axis.ticks = element_blank(), 
     axis.title = element_blank(),
     axis.text = element_blank())+
annotation_custom(grid::textGrob(label = 'Group III',hjust = 0.5,
                                   x=grid::unit(0.8,"npc") ,                                
                                   y=grid::unit(0.9,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 20)))

# my_saveplot_fun(p, 'ETP_STEMNET_age', width = 12, height = 4)

options(repr.plot.width=3.8*3, repr.plot.height=3.8)
p <- wrap_plots(list(p1,p2,p3))+plot_layout(ncol = 3)
p
my_saveplot_fun(p, 'ETP_STEMNET_age', width = 3.8*3, height = 3.8)

####Fig.2g legend
lgd1 = Legend(labels = c('ETP', 'TP1', 'TP2'), title = "Identity", 
              title_gp = gpar(fontsize = 17.5, fontfamily='sans'),
              grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), gap = unit(6, "mm"),
              title_gap = unit(6, "mm"),
              nrow = 1,title_position = 'leftcenter',
              labels_gp = gpar(fontsize = 17.5, fontfamily = 'sans'), 
              legend_gp = gpar(fill = c(cluster_colors[ c('ETP', 'TP1', 'TP2')])))

lgd2 = Legend(labels = c(paste0('Group ', as.roman(1:3))),size = unit(6, "mm"), gap = unit(6, "mm"),
              title_gap = unit(6, "mm"),
              title = "Age", type = "points", 
              title_gp = gpar(fontsize = 17.5, fontfamily='sans'),
              nrow = 1,title_position = 'leftcenter',
              legend_gp = gpar(cex = 5),
              labels_gp = gpar(fontsize = 17.5, fontfamily = 'sans'), 
              pch = c(15,16,17))

pd = packLegend(lgd1, lgd2, #row_gap = unit(1, "cm"), 
                column_gap = unit(0.5, "cm"), direction = "horizontal")





options(repr.plot.width=9, repr.plot.height=0.4)
ggplotify::as.ggplot(grid::grid.grabExpr(
draw(pd,x = unit(0, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))))
ggsave('ETP_STEMNET_legend.svg', width =9, height = 0.4)



###Fig.2h left
plot_df$Age_group <- factor(plot_df$Age_group, levels = c('Prepuberal', 'Adult' ,'Aged'))

tmp_plot_df <- plot_df[plot_df$identity %in% c('ETP' ,'TP1'), c('identity','Age_group',
                                                                'T_lineage', 'ILC')]
tmp_plot_df$Age_group <- age_group_lbls_func(tmp_plot_df$Age_group)
saveRDS(tmp_plot_df, file = 'ETP_TP1_aging_potential.rds')
tmp_plot_df <- plot_df[plot_df$identity %in% c('ETP'),]
stat_test <- wilcox_test(
 T_lineage ~ Age_group, data = tmp_plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.4,fun = 'mean')
stat_test$y.position <- stat_test$y.position + 0.3

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 3)


p1 <- ggplot(tmp_plot_df, aes(Age_group, T_lineage))+
geom_boxplot(mapping = aes(fill = Age_group), outlier.colour = NA)+
scale_x_discrete(label = age_group_lbls_func)+
scale_y_continuous(limits = c(0,1.05), oob = scales::squish)+
geom_point(size = 0.3, alpha = 0.3,stroke = 0, position = 'jitter')+
# geom_point( mapping = aes(x=Age_group, y=T_lineage), 
#            position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors)+
labs(x = NULL,  fill = 'Age')+
guides(fill = F)+
#scale_y_continuous(expand = c(0, 0.005))+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt, tip.length = 0)+
labs(title = 'T-lineage',y = NULL)
options(repr.plot.width=3, repr.plot.height=3)



stat_test <- wilcox_test(
 ILC ~ Age_group, data = tmp_plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.3,fun = 'mean')
stat_test$y.position <- stat_test$y.position + 0.1

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 3)

p2 <- ggplot(tmp_plot_df, aes(Age_group, ILC))+
geom_boxplot(mapping = aes(fill = Age_group), outlier.colour = NA)+
scale_x_discrete(label = age_group_lbls_func)+
geom_point(size = 0.3, alpha = 0.3,stroke = 0, position = 'jitter')+
# geom_point( mapping = aes(x=Age_group, y=ILC), 
#            position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors)+
labs(x = NULL,  fill = 'Age')+
guides(fill = F)+
scale_y_continuous(limits = c(0,0.22), oob = scales::squish)+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt, tip.length = 0)+
labs(title = 'ILC',y = NULL)
options(repr.plot.width=3, repr.plot.height=4)
wrap_plots(list(p1,p2))+plot_layout()

ggsave('ETP_potential_score_aging.svg', width =3, height = 4)



###Fig.2h right
tmp_plot_df <- plot_df[plot_df$identity %in% c('TP1'),]
stat_test <- wilcox_test(
 T_lineage ~ Age_group, data = tmp_plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.25,fun = 'mean')
stat_test$y.position <- stat_test$y.position + 0.1
logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 3)

p1 <- ggplot(tmp_plot_df, aes(Age_group, T_lineage))+
geom_boxplot(mapping = aes(fill = Age_group), outlier.colour = NA)+
scale_x_discrete(label = age_group_lbls_func)+
geom_point(size = 0.3, alpha = 0.3,stroke = 0, position = 'jitter')+
# geom_point( mapping = aes(x=Age_group, y=T_lineage), 
#            position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors)+
labs(x = NULL,  fill = 'Age')+
guides(fill = F)+
#scale_y_continuous(expand = c(0, 0.005))+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt, tip.length = 0)+
labs(title = 'T-lineage',y = NULL)
# options(repr.plot.width=3, repr.plot.height=4.5)





stat_test <- wilcox_test(
 ILC ~ Age_group, data = tmp_plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.25,fun = 'mean')
stat_test$y.position <- stat_test$y.position + 0.07
logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 3)

p2 <- ggplot(tmp_plot_df, aes(Age_group, ILC))+
geom_boxplot(mapping = aes(fill = Age_group), outlier.colour = NA)+
scale_x_discrete(label = age_group_lbls_func)+
geom_point(size = 0.3, alpha = 0.3,stroke = 0, position = 'jitter')+
# geom_point( mapping = aes(x=Age_group, y=ILC), 
#            position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors)+
labs(x = NULL,  fill = 'Age')+
guides(fill = F)+
scale_y_continuous(limits = c(0,0.13), oob = scales::squish)+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt, tip.length = 0)+
labs(title = 'ILC',y = NULL)
options(repr.plot.width=3, repr.plot.height=4)


wrap_plots(list(p1,p2))+plot_layout()
ggsave('TP1_potential_score_aging.svg', width =3, height = 4)





Hs_earlyT_multiomics <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Hs_earlyT_combined/Hs_earlyT_multiomics.rds')

Hs_earlyT <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Hs_earlyT_combined/Hs_earlyT_query.rds')

load('/data1/02.private/dengyj/analysis/thymus/Integrated/Hs_earlyT_combined/Thy1_3_df.RData')

Thy_df <- rbind(Thy1_df, Thy2a_df, Thy2b_df, Thy2c_df, Thy3_df)

tmp_cluster_colors <- c('Thy1' = '#DC143C', 
                        'Thy2b' = '#808000', 
                        'Thy2a' = '#20B2AA',
                        'Thy2c' = '#FFA500',
                       'Thy3' = '#800080')

options(repr.plot.width=3, repr.plot.height=3)

p <- DimPlot(Hs_earlyT, reduction = 'umap',label = F, label.size = 13.75/.pt)+
scale_color_manual(values = tmp_cluster_colors)+
d_theme_w(size = 13.75)+
theme(plot.margin = ggplot2::margin(rep(0), unit = 'cm'),
      axis.ticks = element_blank(),
      axis.text=element_blank(),axis.title = element_blank(),
      legend.text = element_text(size = 13.75, margin = ggplot2::margin(l=-0.3,unit = 'cm')),
      legend.position = c(0.015,0.19))+
guides(color = guide_legend(override.aes = list(size = 3)))#+

p

 

###Fig.2i

options(repr.plot.width=3, repr.plot.height=3)
new_p <- p + geom_mark_hull(linetype = 'dashed',color = 'black',
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

#ggsave('Thy1_3_PID.svg', width = 3.7, height = 3.7) 
my_saveplot_fun(new_p, 'Thy1_3_PID', width = 3, height = 3) 
new_p


###Fig.2j
plot_df <- cbind(data.frame(identity = Idents(Hs_earlyT)), Embeddings(Hs_earlyT, reduction = 'umap'))
plot_df$ID <- as.character(Idents(ETP_TP))[match(rownames(plot_df), colnames(ETP_TP))]
plot_df$ID[is.na(plot_df$ID)] <- 'other'
options(repr.plot.width=3, repr.plot.height=3)

p <- ggplot(mapping = aes(UMAP_1, UMAP_2))+
geom_point(data = plot_df, size = 0.5,stroke = 0,
           mapping = aes(color = ID))+
scale_color_manual(values = cluster_colors, breaks = c('ETP', 'TP1', 'TP2'))+
geom_mark_hull(linetype = 'dashed',color = 'black',
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
                    con.border= 'none')+
d_theme_w(size = 13.75)+
theme(legend.position = c(0.13, 0.17), 
      axis.text=element_blank(),legend.background = element_blank(),
      axis.title = element_blank(), axis.ticks = element_blank())+
guides(color = guide_legend(override.aes = 
                           list(size = 5)))+
labs(color = NULL)
p
my_saveplot_fun(p, 'Thy1_3_ETP_TP_loc', width = 3, height = 3) 







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

p1 <- ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, size = 0.5,mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 13.75)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = gsub('-TotalSeqB', ' Protein', tmp_SM),hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75)))+
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

p2 <- ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, size = 0.5,mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 13.75)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = gsub('-TotalSeqB', ' Protein', tmp_SM),hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75)))+
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
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.95)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.05)

tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

p3 <- ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, size = 0.5,mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 13.75)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = gsub('-TotalSeqB', ' Protein', tmp_SM),hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))


options(repr.plot.width=5, repr.plot.height=4)

options(repr.plot.width=7.5, repr.plot.height=6)
tmp_SM <- 'CD44-TotalSeqB'
tmp_plot_df <- plot_df
max_val <- quantile(tmp_plot_df[, tmp_SM], 0.99)
min_val <- quantile(tmp_plot_df[, tmp_SM], 0.9)
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] > max_val] <- max_val
tmp_plot_df[, tmp_SM][tmp_plot_df[, tmp_SM] < min_val] <- min_val

p4 <- ggplot(tmp_plot_df, aes_string('UMAP_1', 'UMAP_2'))+
geom_point(shape = 21, stroke = 0, size = 0.5,mapping = aes_string(fill = paste0('`',tmp_SM, '`' )))+
scale_fill_gradientn(colours = c("lightgrey", "blue"), 
                        breaks = c(min(tmp_plot_df[,tmp_SM]), 
                                   mean(c(min(tmp_plot_df[,tmp_SM]), max(tmp_plot_df[,tmp_SM]))), 
                                   max(tmp_plot_df[,tmp_SM])),  
    labels = c("Low", "Median", "High"))+
d_theme_w(size  = 13.75)+
theme(#legend.position = c(0.12,0.15)
     )+
labs(title = NULL, fill = NULL)+
annotation_custom(grid::textGrob(label = gsub('-TotalSeqB', ' Protein', tmp_SM),hjust = 0.5,
                                   x=grid::unit(0.25,"npc") ,                                
                                   y=grid::unit(0.93,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75)))+
scale_y_continuous(breaks = c(-4,0,4))+
 geom_mark_hull( linetype = 'dashed',size = 0.7,
                    mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'ident'),data = Thy_df,
                    expand = unit(0.3,'mm'),concavity = 5,
                    radius = unit(1, "mm"), con.type = 'straight',
                    con.border= 'none')+
scale_color_manual(values = tmp_cluster_colors)+
guides(fill = guide_colorbar(order = 1))


###Fig.2k
options(repr.plot.width=7.2, repr.plot.height=6)
p <- wrap_plots(list(p1,p2,p3,p4))+plot_layout(ncol =2, guides = 'collect')&
theme(axis.text = element_blank(), axis.ticks = element_blank(),
     axis.title = element_blank())&
labs(color = 'Identity')
p
my_saveplot_fun(p, 'Thy1_3_defined_surface_markers', width = 7.2, height = 6) 







###Fig.2l
ETP_TP <- AddMetaData(ETP_TP, Idents(Hs_earlyT), col.name = 'PID')

logi <- !is.na(ETP_TP$PID)
plot_df <- data.frame(table(ETP_TP$PID[logi], Idents(ETP_TP)[logi]))
plot_df <- plyr::ddply(plot_df, 'Var2', 'transform', percentage = Freq/sum(Freq))
plot_df$Var1 <- factor(plot_df$Var1, levels = rev(levels(plot_df$Var1)))
options(repr.plot.width=6, repr.plot.height=1.6)
ggplot(plot_df, aes(percentage*100, Var2, fill = Var1)) + geom_bar(stat = 'identity') + 
#scale_x_discrete(limits = c('ETP_TP_1', 'ETP_TP_2', 'ETP_TP_3'), label = c('ETP_TP', 'TP1', 'TP2'))+
labs(y= 'Transcriptional ID', x = '% Proportion', fill = 'Phenotypical ID')+
d_theme_w(size = 13.75)+
theme(plot.margin = margin(t = 0.85,unit = 'cm'),
      axis.text.y = element_text(size =13.75, angle = 0, hjust = 1),legend.position = 'right',
      axis.text.x = element_text(size =13.75),
      axis.title = element_text(size =13.75), 
     legend.title = element_text(size =13.75), 
      legend.text = element_text(size = 13.75, margin = ggplot2::margin(l=-0.2,unit = 'cm')),
     panel.background = element_rect(fill = NA, color = 'black'))+
scale_y_discrete(limits = rev(c('ETP', 'TP1', 'TP2')))+
scale_fill_manual(values = tmp_cluster_colors, limits =rev)+guides(fill = guide_legend(ncol=2))
ggsave('ETP_TP_mapping.svg', width = 6, height = 1.6) 

















