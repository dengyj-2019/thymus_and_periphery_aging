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

library(ggpubr)

plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig1_FigS1'

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



thymus_step_2 <- 
readRDS('/data1/02.private/dengyj/analysis//thymus/clustering_10_20/combined_result/thymus_step_2.rds')

###Extended Data Fig.1a
plot_df <- cbind(FetchData(thymus_step_2, vars = c('anno_1')), 
                 Embeddings(thymus_step_2, reduction = 'umap'))


options(repr.plot.width=5.6, repr.plot.height=4)
p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color=anno_1))+
geom_point(size=1, stroke = 0)+
scale_color_manual(values = colors)+
scale_color_manual(values = c("Endo"="#ec6960","Epi"="#ae80dc",
                              'Mes'="#72bd90","Hema"="#9cbfe7"), 
                  limits = c('Hema', 'Epi', 'Mes', 'Endo'), 
                  label = function(x){
                      x <- gsub('Fibro', 'Mes', x)
                  })+
theme(legend.position = 'right',legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.ticks=element_blank(),
      axis.text=element_blank(),axis.line=element_blank(), axis.title=element_text(size = 17.5),
      legend.text = element_text(size=17.5, family= 'sans'))+
guides(colour = guide_legend(override.aes = list(size=10)))+
labs(color = NULL)
p
my_saveplot_fun(p, 'super_population', width=5.6, height=4)


###Extended Data Fig.1b
p <- VlnPlot(thymus_step_2,stacked=T,pt.size=0, direction = 'horizontal',group.by = 'anno_1',
             features = c('PTPRC', 'RAC2', 'LCK', 'MAL','LAT',
                          'PAX9', 'SIX1','EPCAM','KRT17', 'KRT5','COL1A1', 'PDGFRA','PDGFRB',
                          'FBN1', 'FGF7', 
                          'CDH5',  'CLDN5', 'CD36', 'ESAM', 'ADGRL4'),
       cols = c("Endo"="#ec6960","Epi"="#ae80dc",
                              'Mes'="#72bd90","Hema"="#9cbfe7"), line.size = 0 )+
theme(axis.text.y = element_text(size = 17.5, family = 'sans'),
      axis.ticks.y = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_blank(),
      strip.background.x = element_blank(), 
     strip.text.x = element_text(size = 17.5, hjust = 0, face = 'italic', family = 'sans', angle =90),
      panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     plot.margin=unit(c(0, 0, 0, 1),'lines'))+
scale_x_discrete(limits = rev(c('Hema', 'Epi', 'Mes', 'Endo')),
                 position = 'top')+
scale_y_reverse()
options(repr.plot.width=11, repr.plot.height=4)
p 

ggsave('super_population_DEGs.svg', width=11, height=4)



###Fig.1c
Hema <- readRDS('/data1/02.private/dengyj/analysis//thymus/clustering_10_20/combined_result/Hema.rds')

Hema_colors <- c('Bcell'='#877630','DP'='#d37485','ETP_DN'='#CA9B5D','Innate_T'='#3ca462','Myeloid'='#60a5de',
                 'SP'='#6d69ad','Unconventional_cells'='#f9e21d')#'#eb7814')

Hema$Age_group[is.na(Hema$Age_group)] <- Hema$age_3[is.na(Hema$Age_group)]
Hema$Age_group[Hema$Age_group == '<=12'] <- 'Prepuberal'
Hema$Age_group[Hema$Age_group == '13-39'] <- 'Adult'
Hema$Age_group[Hema$Age_group == '>=40'] <- 'Aged'
plot_df <- cbind(FetchData(Hema, vars = c('hema_class', 'Age_group', 'batch_2')), 
                 Embeddings(Hema, reduction = 'umap'))


options(repr.plot.width=6, repr.plot.height=6)
p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color=hema_class))+
geom_point(size=0.1)+
scale_colour_manual(values=Hema_colors, label = function(x){
    x <- gsub('Unconventional_cells', 'NK & ILC', x)
    x <- gsub('ETP_DN', 'DN', x)                
    x <- gsub('_', '-', x)    
},
                    limits = c('ETP_DN', 'DP', 'SP', 'Innate_T', 'Bcell', 'Myeloid', 'Unconventional_cells'))+
d_theme_w(size = 13.75)+
theme(legend.position = c(0.69,0.1),legend.spacing.y = unit(0.1,'cm'),
      legend.spacing.x = unit(0,'cm'),legend.key = element_blank(),
      panel.background=element_rect(colour = "black"),
      #axis.text = element_blank(),
      legend.background = element_blank(),
      #legend.spacing = unit(0,'cm'),
      #axis.title=element_blank(),
      axis.ticks=element_blank(),
      legend.text = element_text(size=13.75, family='sans'),axis.line=element_blank())+ 
guides(colour = guide_legend(nrow = 3,override.aes = list(size=4), byrow = F, ))+
labs(color = NULL)
p
my_saveplot_fun(p, 'hema_population', width=4.4, height=4.4)

###Fig.1d

p <- VlnPlot(Hema,stacked=T,pt.size=0, group.by = 'hema_class',direction = 'horizontal',
             features = c('HES4','MME', 'HIVEP3','NOTCH1',  'IGLL1','SPINK2',
'AQP3', 'ELOVL4', 'RORC','CD4', 'CD8A','CD8B','CD3D', 'CCR7', 'TOX2', 'CD27' ,
'LY9', 'TRDV1', 'TRGC2','TRGC1', 'ZNF683','PDCD1',
                     'NCR1', 'TYROBP', 'KLRD1','CTSW','KLRF1',  
                     'CD19', 'PAX5','MS4A1', 'CD79A','FCRLA', 'CD24', 
                     'SPI1', 'LYZ', 'IRF8', 
                     'MS4A6A','CD68', 'FCER1G'),
       cols = Hema_colors, line.size = 0 )+
theme(axis.text.y = element_text(size = 17.5, family = 'sans'),
      axis.ticks.y = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_blank(),
      strip.background.x = element_blank(), 
     strip.text.x = element_text(size = 17.5, hjust = 0, face = 'italic', family = 'sans', angle =90),
      panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     plot.margin=unit(c(0, 0, 0, 1),'lines'))+
scale_x_discrete(limits = rev(c('ETP_DN', 'DP', 'SP', 'Innate_T', 'Bcell', 'Myeloid',
                                'Unconventional_cells')),
                 position = 'top',
                label = function(x){
    x <- gsub('Unconventional_cells', 'NK & ILC', x)
    x <- gsub('ETP_DN', 'DN', x)                
    #x <- gsub('Innate_T', 'innate T', x)
    x <- gsub('_', '-', x)
})+
scale_y_reverse()
options(repr.plot.width=11.5, repr.plot.height=6)
p                
ggsave('hema_DEGs.svg', width = 11.5, height = 6)                 
# my_saveplot_fun(p, 'hema_DEGs', width = 11.5, height = 6)                 



###Fig.1e
IFN_response <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/IFN_response.rds')

Bcell_idents <-c('MemoryB','NaiveB','Plasma')
ETP_DN_idents <-c('ETP','DN_Q1', 'DN_Q2','DN4','DN_P','γδT_P','TP2','TP1')
DP_idents<-c('DP_P','DP_Q','HSP_DP', 'αβ_Entry')
CD4_SP_idents<-c('pre_Treg1','CD4+Naive',#'CD4+SOX4+Naive',
                 'CD4+Memory_like','Agonist_2',
     'Treg','CD4+pre_T','pre_Treg2',  'CD4+IFN_response')
CD8_SP_idents<-c('CD8+Naive',#'CD8+SOX4+Naive', 
'CD8+Memory_like','CD8+pre_T','Agonist_1', 'CD8+IFN_response')
Innate_T_idents <- c('γδT_3','ZNF683+CD8ααT','GNG4+CD8ααT','γδT_2','IFN_CD8ααT','γδT_1', 'NKT_like_CD8ααT')
Unconventional_cells_idents<- c('MatureNK','MemoryNK','Lti','ImmatureNK')
Myeloid_idents<-c('DC2','aDC','pDC','Macrophage','DC1P','DC1','Monocyte','Neutrophil','Mast')



Hema$identity2 <- Hema$identity
Hema$identity2[colnames(IFN_response)[Idents(IFN_response) == 'CD4+IFN_response']] <- 'CD4+IFN_response'
Hema$identity2[colnames(IFN_response)[Idents(IFN_response) == 'CD8+IFN_response']] <- 'CD8+IFN_response'
Hema$hema_class2 <- 'other'


Hema$hema_class2[Hema$identity2 %in% Bcell_idents] <- 'Bcell'
Hema$hema_class2[Hema$identity2 %in% ETP_DN_idents] <- 'ETP-DN'
Hema$hema_class2[Hema$identity2 %in% DP_idents] <- 'DP'
Hema$hema_class2[Hema$identity2 %in% CD4_SP_idents] <- 'CD4+SP'
Hema$hema_class2[Hema$identity2 %in% CD8_SP_idents] <- 'CD8+SP'
Hema$hema_class2[Hema$identity2 %in% Unconventional_cells_idents] <- 'Unconventional_cells'
Hema$hema_class2[Hema$identity2 %in% Innate_T_idents] <- 'Innate_T'
Hema$hema_class2[Hema$identity2 %in% Myeloid_idents] <- 'Myeloid'

###remove enriched thymocytes
logi <- !Hema$batch_2 %in% c('fifth', 'sixth')
plot_df <- 
data.frame(apply(table(Hema$hema_class2[logi], Hema$batch[logi]), 2, function(x){x/sum(x)}), check.names = F)
plot_df$cluster <- rownames(plot_df)

plot_df <- melt(plot_df, id.vars = 'cluster', value.name = 'pct', variable.name = 'batch')
plot_df$age_3 <- Hema$age_3[match(plot_df$batch, Hema$batch)]
plot_df$age <- Hema$age[match(plot_df$batch, Hema$batch)]
plot_df$pct <- plot_df$pct*100
write.xlsx(plot_df, file = 'scRNAseq_hema_cells_pct_aging.xlsx')

saveRDS(plot_df, file = 'scRNAseq_hema_cells_pct_aging.rds')

clu <- c('ETP-DN', 'DP', 'CD4+SP', 'CD8+SP', 'Bcell', 'Myeloid')
plot_list <- lapply(clu, function(x){
    tmp_df <- plot_df[plot_df$cluster == x, ]
    pairs <- get_pairs(unique(tmp_df$age_3))
    aov_res <- aov(pct ~ age_3, data = tmp_df)
    stat_df <- rstatix::pairwise_t_test(pct~age_3, data = tmp_df, pool.sd = F)
    stat_df$p.signif <- 'ns'
    stat_df$p.signif[stat_df$p < 0.05] <- '*'
    stat_df$p.signif[stat_df$p < 0.01] <- "**"
    stat_df$p.signif[stat_df$p < 0.001] <- "***"
    stat_df$p.signif[stat_df$p < 0.0001] <- "****"
    logi <- stat_df$p.signif == 'ns'
    stat_df$p.signif[logi] <- stat_df$p[logi]
    max_val <- max(tmp_df$pct)
    len <- nrow(stat_df)
    by_val <- 0.2* max_val
    stat_df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = len)
    p <- ggplot(tmp_df, aes(age_3, pct))+
    
  stat_summary(fun.data = 'get_sem', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge( .9))+
    stat_summary(color = 'black',
        fun.y = mean, 
        geom = "errorbar", 
        aes(ymax = ..y.., ymin = ..y..), 
        width = 0.5)+
    scale_color_manual(values= age_colors, label = age_group_lbls_func)+
     scale_x_discrete(label = age_group_lbls_func)+
    stat_pvalue_manual(stat_df, label = "p.signif", label.size = 15/.pt)+
#     stat_compare_means(comparisons = pairs, method="wilcox.test",  size = 13.75/.pt, vjust = 0.15, tip.length = 0,
#                        step.increase = 0.07,
#                        label ='p.signif',
#                        label.x.npc = 0.3)+
  
    d_theme_w(size = 15, italic_text = '')+
    theme(panel.border = element_rect(fill = NA, color = NA),axis.line = element_line(),
        plot.title = element_text(hjust = 0.5),plot.margin = margin(t = 0, unit = 'cm'),
          axis.text = element_text(hjust = 0.5,size = 15))+      
    guides(color = F,shape = F)+
    labs(y='of CD45+ (%)', x=  NULL, color = 'Age group', title = x)+
    geom_point(mapping = aes(color = age_3, shape = age_3))+
    scale_shape_manual(values = c('<=12' = 16, '13-39' = 15, '>=40' = 17))+
    scale_y_continuous(expand = expansion(mult = .1))
    return(p)
})
options(repr.plot.width=3.2*6, repr.plot.height=3.5)
wrap_plots(plot_list)+plot_layout(nrow = 1)






