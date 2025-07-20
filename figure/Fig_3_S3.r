library(ggplot2)
library(Seurat)

library(ggraph)
library(tidygraph)

library(data.table)

library(ggpubr)


library(patchwork)

library(ComplexHeatmap, lib.loc = .libPaths()[1])

library(circlize)

library(reshape2)

library(MySeuratWrappers)

library(Matrix)

library(openxlsx)

plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot//Fig3_FigS3'

setwd(plot_dir)

cluster_colors <- 
c(
    'CD8+Naive'='#53A85F',
                             'CD4+Naive'='#D6E7A3',
'ETP_1' = '#297A6B','ETP_2'='#C05EA5','ETP'='#a7aad5',
    'ETP_3'='#e6ba77','DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA',
    'DN4'='#ff927c','γδT_P'='#ee587d','pre_Treg1'='#E5D2DD','Treg_diff'='#E5D2DD',
    #'CD8+Naive'='#53A85F',
                    'CD4+Memory_like'='#808040','CD8+Memory_like'='#F3B1A0',
                             #'CD4+Naive'='#D6E7A3',
                    'CD4+GZMK+T'='#57C3F3','IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3',
    'Agonist_2'='#a83117','T_agonist'='#a83117',
    'Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'Agonist_1' = '#509cbe','αβ_Entry_2' = '#509cbe', #'#a6c050',
                    'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'pre_Treg2' = '#9566ac','Treg_diff_2'='#9566ac',
    'Tex_like'='#9566ac',
    'CD8+SOX4+Naive'='#53A85F','CD8+SOX4-Naive' = "#F2C43D", 
                              'CD4+SOX4-Naive' ="#D76483",
                             'CD4+SOX4+Naive'='#D6E7A3' ,
    'CD8ααT' = '#e1812b',
   'MemoryB' = '#e0b139', 
    'NaiveB' = '#698c68', 'NKT'='#786085',
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
'PB_CD4+TN'='#7bacd9','PB_CD8+TN'='#DDC4DA',
'PB_CD4+RTE'='#66c5b9'
   )

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode//plot_helper.r')




source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

T_dev_3 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/T_dev_3.rds')

levels(T_dev_3) <- 
c('DP_P','DP_Q','HSP_DP', 'αβ_Entry','Agonist_1', 'CD4+pre_T', 'CD8+pre_T', 'Agonist_2', 
                    'CD4+Naive', 
'CD8+Naive',  'IFN_response', 'pre_Treg1', 'pre_Treg2', 
  'Treg','Tex_like', 'CD4+Memory_like',
                    'CD8+Memory_like')

T_dev_3$identity <- Idents(T_dev_3)
plot_df <- as.data.frame(cbind(FetchData(T_dev_3, 'identity'), Embeddings(T_dev_3, 'umap')))

###Fig.3a
options(repr.plot.width=5, repr.plot.height=7)
p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = identity))+
geom_point(size = 0.5)+
theme(legend.position = c('bottom'),panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      axis.title=element_text(size = 13.75, family = 'sans',color='black'),axis.ticks=element_blank(),
      axis.text=element_blank(),axis.line=element_blank(),legend.key = element_rect(fill= 'white'), 
      legend.title=element_text(size = 13.75, margin = margin(t = -0.5, unit = 'cm')),
     legend.text = element_text(size = 13.75, family = 'sans', 
                                color='black', 
                               margin = margin(l = -0.2, unit = 'cm')))+
#scale_y_reverse()+
guides(colour = guide_legend(override.aes = list(size=5,alpha=1), ncol=3, title.position = 'top')) +
scale_color_manual(values=cluster_colors,label = function(x){
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('CD4.*Memory_like', 'memory-like CD4 T', x)
    x <- gsub('CD8.*Memory_like', 'memory-like CD8 T', x)
    #x <- gsub('Agonist', 'agonist', x)
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x
})+labs(color = 'Identity')#+scale_y_reverse()+scale_x_reverse()
p
my_saveplot_fun(p, 'T_dev_idents', width=5,height=7)



###Fig.3b
p <- DotPlot(T_dev_3, features = rev(c('CD4', 'CD40LG','CD8A', 'CD8B','MKI67','TCF19','CDK1', 'ELOVL4','DNTT', 'RAG1',
                                'RORC', 'CD1E','DNAJB1', 'CHORDC1','HSPH1', 'HSPA1A',  'CD28','TRAC',
                                       'TOX','TOX2', 
                                'BACH2',
                                'CCR9', 'CCR4','ZBTB7B','RUNX3',  
                                'TCF7',  'CCR7','SELL','IL7R', 'KLF2','S1PR1', 
                                'IKZF2','HIVEP3', 'NFATC1','NR4A3','PDCD1', 
                                'CD69', 'NR4A1','EGR2','MIR155HG','ICOS', 'IL2RA',
                                'FOXP3','CTLA4','TNFRSF4','CD27','LAG3',
                                'FAS','PRDM1', 'CD58',
                                'ANXA1', 'AHNAK', 'PRF1','GZMK', 'NKG7','CRTAM','IFIT1', 'OAS1')),

       cols = c('lightgrey', 'red'))

new_p <- ggplot(p$data, aes(features.plot,id,  fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21)+
scale_size_continuous(range = c(0,5))+
scale_fill_gradientn(colours = c('lightgrey', 'red')) + 
theme(axis.text.y = element_text(size = 13.75, family = 'sans', hjust = 1), 
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'horizontal',
      legend.position = 'bottom',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_text(size = 13.75, family = 'sans', face = 'italic',angle =60, hjust =1, vjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),
        legend.margin = margin(l = 2, r = 1, unit = 'cm'),
      axis.line =element_blank(),
     plot.margin  = unit(c(0,0,0,2), 'lines'))+
guides(
       fill = guide_colorbar(title = "Expression", order = 1, 
                            label.position = "bottom", 
            title.position = "left", title.vjust=0.9),  
size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'), order = 2, 
                   label.position = "bottom", title.vjust=0.9 ))+
scale_y_discrete(limits = rev(levels(T_dev_3)), position = 'right',label = function(x){
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('CD4.*Memory_like', 'memory-like CD4 T', x)
    x <- gsub('CD8.*Memory_like', 'memory-like CD8 T', x)
    #x <- gsub('Agonist', 'agonist', x)
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x
})
options(repr.plot.width=13, repr.plot.height=7)
new_p
#my_saveplot_fun(new_p, 'T_dev_3_CD8AAT_DEGs', width = 11, height = 5)
ggsave('T_dev_3_DEGs.svg', width = 13, height = 7)







load('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/ETP_SP_Dev.rds')

###Fig.3c
Dev_2_FDG_df <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/Dev_2_FDG_df.csv', 
                         row.names= 1)



options(repr.plot.width=7.5, repr.plot.height=7.5)
connection <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/Dev_2_paga_con_df.csv', 
                       row.names = 1, check.names = F)
pos <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/ETP_SP/Dev_2_paga_pos_df.csv',
                row.names = 1)
colnames(pos) <- c('x', 'y')
connection$Identity <- rownames(connection)
plot_graph <- melt(connection, id.vars = 'Identity', variable.name = 'Identity2', value.name = 'Connectivity')
plot_graph$threshold <- 'no'
plot_graph$threshold[plot_graph$Connectivity > 0.33] <- 'yes' #0.13, 0.11
graph <- as_tbl_graph(plot_graph, directed = F) %>%
  mutate(size = as.data.frame(table(Dev_2$identity_2), row.names = 1)[name, ])



options(repr.plot.width=6.85, repr.plot.height=5)
lmls <- c('DP_Q','HSP_DP', 'αβ_Entry','αβ_Entry_2','CD4+pre_T', 'CD8+pre_T', 'T_agonist', 
                    'CD4+Naive', 'CD8+Naive', 'CD8ααT', 
                     'IFN_response',#'CD4+Memory',
          'Treg_diff','Treg_diff_2', 'Treg','NKT_like_CD8ααT'
                    )
p <- ggraph(graph, layout = pos)+ 
  geom_edge_link(aes(width = Connectivity, colour = threshold, linetype = threshold), show.legend = F)+ 
  scale_edge_width_continuous(range = c(0.25, 2)) + 
scale_edge_color_manual(values = c('yes' = '#000000', 'no' = 'lightgrey')) + 
  scale_edge_linetype_manual(values = c('yes' = 'solid','no' = 'dashed')) + 
  geom_node_point(aes(#size = size, 
                      fill = name), shape=21, color = 'black', size = 8) + 
  scale_size(range=c(5, 10))  + 
scale_fill_manual(values = cluster_colors) + 
scale_color_manual(values = cluster_colors, limits = lmls,label = function(x){
    x <- gsub('Memory', 'Memory-like', x)
    x <- gsub('αβ_Entry_2', 'Agonist-1', x)
    x <- gsub('T_agonist', 'Agonist-2', x)
    x <- gsub('Treg_diff_2', 'pre-Treg2', x)
    x <- gsub('Treg_diff', 'pre-Treg1', x)
    x <- gsub('NKT_like_CD8ααT', 'NKT', x)   
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('CD4.*Memory_like', 'memory-like CD4 T', x)
    x <- gsub('CD8.*Memory_like', 'memory-like CD8 T', x)
    #x <- gsub('Agonist', 'agonist', x)                      
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)                  
    x
}) + d_theme_w(size =13.75)+
  theme(panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title=element_text(size=13.75, color ='black'),
       ) + 
# coord_flip()+
labs(x = 'FA_y', y = 'FA_x')+
guides(size = F, fill = F, color = guide_legend(override.aes = list(size = 5)))+
  labs(colour = 'Identity') +
geom_point(data = Dev_2_FDG_df,mapping = aes(x,y, color = identity_2),size = 0.5)


p$layers <- c(p$layers[[length(p$layers)]], p$layers[1:(length(p$layers)-1)])

p
my_saveplot_fun(p,'T_Dev_2_PAGA', width = 6.85,height = 5)

###Fig.3d
options(repr.plot.width=5, repr.plot.height=5)
lmls <- c('DP_Q','HSP_DP', 'αβ_Entry','αβ_Entry_2','CD4+pre_T', 'CD8+pre_T', 'T_agonist', 
                    'CD4+Naive', 'CD8+Naive', 'CD8ααT', 
                     'IFN_response',#'CD4+Memory',
          'Treg_diff','Treg_diff_2', 'Treg','NKT_like_CD8ααT'
                    )

p <- ggraph(graph, layout = pos)+ 
  geom_edge_link(aes(width = Connectivity, colour = threshold, linetype = threshold), show.legend = F)+ 
  scale_edge_width_continuous(range = c(0.25, 2)) + 
scale_edge_color_manual(values = c('yes' = '#000000', 'no' = 'lightgrey')) + 
  scale_edge_linetype_manual(values = c('yes' = 'solid','no' = 'dashed')) + 
  geom_node_point(aes(#size = size, 
                      fill = name), shape=21, color = 'black', size = 8) + 
  scale_size(range=c(5, 10))  + scale_fill_manual(values = cluster_colors) + 
# scale_color_manual(values = cluster_colors) + 
scale_color_gradientn(colours = viridis::plasma(20), limits = c(0, 0.5), oob=scales::squish)+
  theme(panel.background = element_rect(color = NA, fill = 'white'), legend.position = c(0.8, 0.2),
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title=element_text(size=13.75, color ='black'),
        legend.title=element_text(size=13.75, color ='black'),
        legend.text=element_text(size=13.75, color ='black'),
       ) + 
labs(x = 'FA_x', y = 'FA_y')+
guides(size = F, fill = F)+
  labs(colour = 'Pseudotime') +
geom_point(data = Dev_2_FDG_df,mapping = aes(x,y, color = dpt_pseudotime),size = 0.5)

p$layers <- c(p$layers[[length(p$layers)]], p$layers[1:(length(p$layers)-1)])

p
my_saveplot_fun(p,'T_Dev_3_PAGA_time', width = 5,height = 5)

###Extended Data Fig.3a
load('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/unconventional/unconventional.rds')

CD8aaT_NKT <- unconventional_3[, grep('CD8ααT|NKT', Idents(unconventional_3))]
CD8aaT_NKT$idents <- as.character(Idents(CD8aaT_NKT))
CD8aaT_NKT$idents[grepl('NKT', CD8aaT_NKT$idents)] <- 'NKT'
CD8aaT_NKT$idents[grepl('CD8ααT', CD8aaT_NKT$idents)] <- 'CD8ααT'
Idents(CD8aaT_NKT) <- CD8aaT_NKT$idents

CD8_CT_NCT <- merge(CD8aaT_NKT, T_dev_3[, Idents(T_dev_3) %in% c('CD8+pre_T', 'CD8+Naive')])

levels(CD8_CT_NCT) <- c('CD8+pre_T', 'CD8+Naive', 'CD8ααT', 'NKT')



options(repr.plot.width=3.5, repr.plot.height=12)
p <- VlnPlot(CD8_CT_NCT,stacked=T,pt.size=0, direction = 'vertical',cols =  cluster_colors,features = c(
    'CD3G','CD3D', 'CD8A', 'CD8B', 'CCR9','CD1E','CD1A','TOX2','SOX4','CCR7','SELL','LRRN3','ACTN1','IL7R','CD27',
     'MME', 'PDCD1','TRGC1','TRGC2','ZNF683','XCL1', 'KLRB1', 'NKG7','KLRD1','GZMH','FCGR3A','GZMB'), 
             line.size = 0 )+
theme(axis.text.x = element_text(size = 13.75, family = 'sans'),strip.background.y = element_blank(),
      axis.ticks.x = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_blank(),#element_text(size = 13.75),
     strip.text.y = element_text(size = 13.75, hjust = 0, face = 'italic', family = 'sans', angle =0),
      panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     plot.margin=unit(c(0, 0, 0, 2.5),'lines'))+
scale_x_discrete(label = function(x){ 
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x
})

p     
#my_saveplot_fun(p, 'CD8_CT_NCT_DEGs', width = 3, height = 12)
ggsave( 'CD8_CT_NCT_DEGs.svg', width = 3.5, height = 12)





###Extended Data Fig.3b
plot_df = plyr::ddply(data.frame(table(Idents(T_dev_3), T_dev_3$age_3)), 
                      'Var2', transform,percentage=Freq/sum(Freq) *100)

colnames(plot_df) <- c('Identity', 'Age', 'Freq', 'Percentage')
plot_df$Identity <- factor(plot_df$Identity)
p <- ggplot(plot_df,aes(x=Age,y=Percentage,
                   fill=Identity,
                   group=Identity))+
geom_area(colour="black",size=0.1)+


scale_fill_manual(values = cluster_colors, 
                          label = function(x){
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('CD4.*Memory_like', 'memory-like CD4 T', x)
    x <- gsub('CD8.*Memory_like', 'memory-like CD8 T', x)
    #x <- gsub('Agonist', 'agonist', x)
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x
})+
theme(plot.margin = margin(l = 0, t = 0.2, unit = 'cm'),
      panel.border = element_blank(),#element_rect(color = 'black', fill = NA), 
      panel.background = element_blank(),#element_rect(color = NA, fill = 'white'), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),
      axis.ticks = element_blank(),
    axis.title.x=element_text(size=13.75, family='sans'),
    axis.title.y=element_text(size=13.75, angle =90, family='sans'),
    axis.text.x=element_text(size=13.75, family='sans', angle =45, hjust=1, vjust = 1),
    axis.text.y=element_text(size=13.75, family='sans'),
    legend.text = element_text(size = 13.75, family='sans'),legend.spacing.y = unit(1e-3, 'npc'),
    legend.title = element_text(size = 13.75, family='sans')
)+guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 0.1)))+
scale_x_discrete(limits = c('<=12', '13-39','>=40'), expand=c(0,0),
                 label = age_group_lbls_func)+
scale_y_continuous( expand = c(0, 0))+
labs(x = NULL, y = 'Percentage %')
options(repr.plot.width=5, repr.plot.height=5)
p 
#my_saveplot_fun(p, 'T_dev_3_proportion', width = 5, height = 5)
ggsave('T_dev_3_proportion.svg', width = 5, height = 5)



###Extended Data Fig.3c
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



Age_Gene_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/thymus/DP_SP/Age_Gene_test_list.rds')


DP_P_df <- get_Age_test_T(Age_Gene_test_list[['DP_P']])




DP_Q_df <- get_Age_test_T(Age_Gene_test_list[['DP_Q']])





DP_P_selected_genes <- c('IL7R','CXCR4','CCR9','SYK',#'TPT1',
                         'TOX2','GRAP2',#'RAG1',
                         'RAG2','TRBC1','FKBP5',
'TSC22D3','FOS','JUN','ID3','GAS2','KLF6','RMI2','DUSP1','PIK3IP1','MRGBP','SHLD1','CXCR3')

DP_Q_selected_genes <- c('TSC22D3','FKBP5','PIK3IP1','TOX2','SHLD1','ID3','JUN','S100A10','JUNB',
'KLF6','IL7R','CXCR4','CCR9','TPT1','GRAP2','RAG2','TRBC1','FOS','GAS2')


DP_P_plot_df <- DP_P_df[DP_P_df$variable %in% DP_P_selected_genes, ]

DP_P_plot_df <- plyr::ddply(DP_P_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DP_P_plot_df$comparison <- gsub('_VS_Rest', '', DP_P_plot_df$comparison)

DP_P_label_df <- data.frame(Age = unique(DP_P_plot_df$comparison))



DP_P_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DP_P_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



DP_P_label_plot <- ggplot(DP_P_label_df, aes(Age,1, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'DP-P')

options(repr.plot.width=14*1/3, repr.plot.height=11)



DP_Q_plot_df <- DP_Q_df[DP_Q_df$variable %in% DP_Q_selected_genes, ]

DP_Q_plot_df <- plyr::ddply(DP_Q_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DP_Q_plot_df$comparison <- gsub('_VS_Rest', '', DP_Q_plot_df$comparison)

DP_Q_label_df <- data.frame(Age = unique(DP_Q_plot_df$comparison))



DP_Q_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DP_Q_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


DP_Q_label_plot <- ggplot(DP_Q_label_df, aes(Age,1, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'DP-Q')



combined_p_val_adj <- c(DP_P_plot_df$p_val_adj, DP_Q_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
DP_P_plot_df$p_val_adj[DP_P_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
DP_Q_plot_df$p_val_adj[DP_Q_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- 
as.numeric(sprintf("%0.0f",seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5)/10))*10

DP_P_plot_df$p_val_adj <- scales::squish(DP_P_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))
DP_Q_plot_df$p_val_adj <- scales::squish(DP_Q_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))


combined_exp <- c(DP_P_plot_df$scale_exp, DP_Q_plot_df$scale_exp)

combined_exp_break <- as.numeric(sprintf("%0.1f",  seq(min(combined_exp), max(combined_exp), length.out = 5)))

DP_P_p <- ggplot(DP_P_plot_df, aes(comparison, variable, #size = -log10(p_val_adj), 
                                   fill = scale_exp))+
geom_tile()+
# scale_size_continuous(range = c(0.5,8), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f16e62'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(DP_P_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(
       fill = guide_colorbar(title = "Expression", order = 1)#,
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2)
)+
labs(x='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =0, hjust = 0), 
      plot.margin = margin(t = 0, unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans',face = 'italic'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),


      axis.line =element_blank())




DP_Q_p <- ggplot(DP_Q_plot_df, aes( comparison,variable, #size = -log10(p_val_adj), 
                                   fill = scale_exp))+
geom_tile()+

scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f16e62'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(DP_Q_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(
       fill = guide_colorbar(title = "Expression", order = 1)#,

)+
labs(x='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =0, hjust = 0), 
      plot.margin = margin(t = 0, unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans', face = 'italic'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=6, repr.plot.height=5)
p <- wrap_plots(DP_P_label_plot, DP_P_p,  DP_Q_label_plot,DP_Q_p, byrow = F)+
plot_layout(ncol = 2, guide = 'collect', heights = c(1,19))&
theme(plot.margin = margin(l = 0.36, unit = 'cm'))
p
ggsave("DP_P_Q_aging_DEGs.svg", width = 6, height = 5)




###Extended Data Fig.3d
Age_GS_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Enrichment/ssGSEA/thymus/Age_GS_test_list.rds')

DP_P_df <- get_Age_test_T(Age_GS_test_list[['DP_P']])
DP_Q_df <- get_Age_test_T(Age_GS_test_list[['DP_Q']])




DP_P_GS <- c('Apoptotic cleavage of cellular proteins',
'GO_POSITIVE_REGULATION_OF_CELL_AGING',
'GO_NEGATIVE_REGULATION_OF_CELL_DEVELOPMENT',
'Eukaryotic Translation Elongation',
'GO_EXPORT_ACROSS_PLASMA_MEMBRANE',
#'GO_K63_LINKED_POLYUBIQUITIN_MODIFICATION_DEPENDENT_PROTEIN_BINDING',
'GO_REGULATION_OF_TELOMERE_MAINTENANCE',
#'NOTCH2 intracellular domain regulates transcription',
'GO_V_D_J_RECOMBINATION',
'GO_T_CELL_HOMEOSTASIS',
'GO_POSITIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION',
#'GO_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H_OXYGEN_AS_ACCEPTOR',
'GO_REACTIVE_OXYGEN_SPECIES_BIOSYNTHETIC_PROCESS',
'KEGG_GLYCOLYSIS_GLUCONEOGENESIS',
'GO_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE', 

'GO_CELL_CYCLE_DNA_REPLICATION',
'GO_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION',
'GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
'GO_T_CELL_APOPTOTIC_PROCESS',
'GO_GLYCOLYTIC_PROCESS',
'GO_REGULATION_OF_T_CELL_DIFFERENTIATION_IN_THYMUS',
'GO_CELLULAR_SENESCENCE',
'GO_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE',
'GO_MHC_CLASS_II_BIOSYNTHETIC_PROCESS')


DP_Q_GS <-  c('Eukaryotic Translation Elongation',
'GO_REGULATION_OF_TELOMERE_MAINTENANCE',
'Apoptotic cleavage of cellular proteins',
'GO_POSITIVE_REGULATION_OF_CELL_AGING',
'KEGG_ASTHMA',
'GO_T_CELL_RECEPTOR_COMPLEX',
#'Phosphorylation of CD3 and TCR zeta chains',
'GO_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION',
'TCR signaling',
'GO_MHC_PROTEIN_BINDING',
'GO_T_CELL_HOMEOSTASIS',
'GO_CELLULAR_SENESCENCE',

#'GO_PROGRAMMED_CELL_DEATH_IN_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES',
'CD28 co-stimulation',
'GO_REGULATION_OF_T_CELL_CHEMOTAXIS',              
'GO_T_CELL_PROLIFERATION',
'GO_NADH_DEHYDROGENASE_ACTIVITY', 
'GO_REGULATION_OF_CELL_CYCLE_CHECKPOINT',

'GO_T_CELL_PROLIFERATION',
'GO_POSTREPLICATION_REPAIR',
'GO_GLYCOLYTIC_PROCESS',
'GO_REGULATION_OF_T_CELL_CHEMOTAXIS',
'GO_GLYCEROLIPID_METABOLIC_PROCESS')

DP_P_plot_df <- DP_P_df[DP_P_df$variable %in% DP_P_GS, ]

DP_P_plot_df <- plyr::ddply(DP_P_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DP_P_plot_df$comparison <- gsub('_VS_Rest', '', DP_P_plot_df$comparison)

DP_P_label_df <- data.frame(Age = unique(DP_P_plot_df$comparison))



DP_P_GS_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DP_P_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



DP_P_label_plot <- ggplot(DP_P_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, hjust = 0.5),
            plot.margin = margin(t = 0.3, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'DP-P')

options(repr.plot.width=14*1/3, repr.plot.height=11)



DP_Q_plot_df <- DP_Q_df[DP_Q_df$variable %in% DP_Q_GS, ]

DP_Q_plot_df <- plyr::ddply(DP_Q_plot_df, c('variable'), transform, scale_exp = scale(means_1))
DP_Q_plot_df$comparison <- gsub('_VS_Rest', '', DP_Q_plot_df$comparison)

DP_Q_label_df <- data.frame(Age = unique(DP_Q_plot_df$comparison))



DP_Q_GS_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(DP_Q_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


DP_Q_label_plot <- ggplot(DP_Q_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, hjust = 0.5),
            plot.margin = margin(t = 0.3, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'DP-Q')



combined_p_val_adj <- c(DP_P_plot_df$p_val_adj, DP_Q_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
DP_P_plot_df$p_val_adj[DP_P_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
DP_Q_plot_df$p_val_adj[DP_Q_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- 
as.numeric(sprintf("%0.0f",seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5)/10))*10

DP_P_plot_df$p_val_adj <- scales::squish(DP_P_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))
DP_Q_plot_df$p_val_adj <- scales::squish(DP_Q_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))


combined_exp <- c(DP_P_plot_df$scale_exp, DP_Q_plot_df$scale_exp)

combined_exp_break <- as.numeric(sprintf("%0.1f",  seq(min(combined_exp), max(combined_exp), length.out = 5)))

DP_P_p <- ggplot(DP_P_plot_df, aes(comparison, variable, size = -log10(p_val_adj), color = scale_exp))+
geom_point()+
scale_size_continuous(range = c(0.5,10), 
                     breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_color_gradientn(colours = c("#481055","#32538f", "#00918c", '#00be63', '#fae800'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(DP_P_GS_order), label = function(x){
    x <- gsub('GO_|KEGG_', '', x)
    x <- gsub('_', ' ', x)
    x <-  str_to_sentence_func(x)
},position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), label = age_group_lbls_func)+
guides(
       color = guide_colorbar(title = "Enrichment\nscore", order = 1),
size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2))+
labs(X='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =45, hjust = 1), 
      plot.margin=unit(c(0, 0, 0, 0),'lines'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans',hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=14*1/3, repr.plot.height=11)
DP_Q_p <- ggplot(DP_Q_plot_df, aes(comparison, variable, size = -log10(p_val_adj), color = scale_exp))+
geom_point()+
scale_size_continuous(range = c(0.5,10), 
                     breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_color_gradientn(colours = c("#481055","#32538f", "#00918c", '#00be63', '#fae800'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(DP_Q_GS_order), label = function(x){
    x <- gsub('GO_|KEGG_', '', x)
    x <- gsub('_', ' ', x)
    x <-  str_to_sentence_func(x)
},position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), label = age_group_lbls_func)+
guides(
       color = guide_colorbar(title = "Enrichment\nscore", order = 1),
size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2))+
labs(X='Age', y=NULL)+
theme(
    axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =45, hjust = 1), 
    plot.margin=unit(c(0, 0, 0, 0),'lines'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans',hjust = 1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=8, repr.plot.height=12)
p <- wrap_plots(DP_P_label_plot,  DP_P_p, DP_Q_label_plot,DP_Q_p)+
plot_layout(ncol = 1, guide = 'collect', heights = c(1,19, 1,19))
p
my_saveplot_fun(p, 'DP_P_Q_age_genesets', width = 8, height = 12)












###Fig.3e
thymus_Naive_age_GS_combined_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Enrichment/ssGSEA/thymus/Naive/thymus_Naive_age_GS_combined_list.rds')

CD8_Naive_df <- get_Age_test_T(thymus_Naive_age_GS_combined_list[['CD8+Naive']])
CD4_Naive_df <- get_Age_test_T(thymus_Naive_age_GS_combined_list[['CD4+Naive']])


CD4_Naive_GS <- c('Mitotic Anaphase',
'GO_NUCLEOSOMAL_DNA_BINDING',
'TCR signaling',
'Signaling by NOTCH',
'Golgi-to-ER retrograde transport',
'Signaling by WNT',
'DNA Repair',
'ISG15 antiviral mechanism',
#'TCF dependent signaling in response to WNT',

'GO_CYTOSKELETON_DEPENDENT_CYTOKINESIS',

'GO_NEGATIVE_REGULATION_OF_MICROTUBULE_POLYMERIZATION',
'PIP3 activates AKT signaling',
#'GO_CELL_CYCLE_G2_M_PHASE_TRANSITION',

'Interleukin-27 signaling',
'GO_RESPONSE_TO_OXYGEN_LEVELS',
'Programmed Cell Death',
'MAPK family signaling cascades',

'GO_REGULATION_OF_T_HELPER_CELL_DIFFERENTIATION',                  
'GO_INTERLEUKIN_18_PRODUCTION',
'Activation of the AP-1 family of transcription factors',
'GO_INTERLEUKIN_6_PRODUCTION',
'Interleukin-10 signaling',
'GO_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION',
'GO_INTERLEUKIN_8_SECRETION',                  
'GO_INTERLEUKIN_17_PRODUCTION',
#'GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY',
'GO_INTERLEUKIN_1_SECRETION',
'GO_T_HELPER_17_TYPE_IMMUNE_RESPONSE',
'GO_MHC_PROTEIN_COMPLEX',                  
'GO_INTERFERON_GAMMA_PRODUCTION', 
'GO_T_CELL_HOMEOSTASIS',
'GO_T_CELL_PROLIFERATION',
'GO_T_CELL_DIFFERENTIATION_IN_THYMUS'                 )


CD8_Naive_GS <- c('GO_DNA_PACKAGING_COMPLEX',
'GO_LIPOPEPTIDE_BINDING',
#'GO_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
#'GO_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY',
'GO_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM',
'GO_ENDOCYTIC_VESICLE_LUMEN',
'GO_NUCLEOSOME_BINDING',
'GO_CYTOSKELETAL_ADAPTOR_ACTIVITY',                  
'GO_THYMOCYTE_APOPTOTIC_PROCESS',
'Golgi-to-ER retrograde transport',
'GO_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN',
'GO_CELLULAR_RESPONSE_TO_CALCIUM_ION',
'GO_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING',
'Activation of the AP-1 family of transcription factors',
'GO_NEGATIVE_REGULATION_OF_CHEMOTAXIS',
'GO_MHC_PROTEIN_COMPLEX',
'GO_NEGATIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION',
'GO_NEGATIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS',
'rRNA processing in the mitochondrion',
'Eukaryotic Translation Elongation',
'GO_INTERFERON_GAMMA_PRODUCTION', 
'TRAF6 mediated NF-kB activation',
#'GO_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE',
'GO_NEGATIVE_REGULATION_OF_CELL_KILLING')

CD4_Naive_plot_df <- CD4_Naive_df[CD4_Naive_df$variable %in% CD4_Naive_GS, ]

CD4_Naive_plot_df <- plyr::ddply(CD4_Naive_plot_df, c('variable'), transform, scale_exp = scale(means_1))
CD4_Naive_plot_df$comparison <- gsub('_VS_Rest', '', CD4_Naive_plot_df$comparison)

CD4_Naive_label_df <- data.frame(Age = unique(CD4_Naive_plot_df$comparison))



CD4_Naive_GS_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(CD4_Naive_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



CD4_Naive_label_plot <- ggplot(CD4_Naive_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'mature CD4 T')

options(repr.plot.width=14*1/3, repr.plot.height=11)



CD8_Naive_plot_df <- CD8_Naive_df[CD8_Naive_df$variable %in% CD8_Naive_GS, ]

CD8_Naive_plot_df <- plyr::ddply(CD8_Naive_plot_df, c('variable'), transform, scale_exp = scale(means_1))
CD8_Naive_plot_df$comparison <- gsub('_VS_Rest', '', CD8_Naive_plot_df$comparison)

CD8_Naive_label_df <- data.frame(Age = unique(CD8_Naive_plot_df$comparison))



CD8_Naive_GS_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(CD8_Naive_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


CD8_Naive_label_plot <- ggplot(CD8_Naive_label_df, aes(Age, 1,fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'mature CD8 T')



combined_p_val_adj <- c(CD4_Naive_plot_df$p_val_adj, CD8_Naive_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
CD4_Naive_plot_df$p_val_adj[CD4_Naive_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
CD8_Naive_plot_df$p_val_adj[CD8_Naive_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- 
as.numeric(sprintf("%0.0f",seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5)/10))*10

CD4_Naive_plot_df$p_val_adj <- scales::squish(CD4_Naive_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))
CD8_Naive_plot_df$p_val_adj <- scales::squish(CD8_Naive_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))


combined_exp <- c(CD4_Naive_plot_df$scale_exp, CD8_Naive_plot_df$scale_exp)

combined_exp_break <- as.numeric(sprintf("%0.1f",  seq(min(combined_exp), max(combined_exp), length.out = 5)))

CD4_Naive_p <- ggplot(CD4_Naive_plot_df, aes(comparison, variable, size = -log10(p_val_adj), color = scale_exp))+
geom_point()+
scale_size_continuous(range = c(0.5,6), 
                     breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_color_gradientn(colours = c("#481055","#32538f", "#00918c", '#00be63', '#fae800'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(CD4_Naive_GS_order), label = function(x){
    x <- gsub('GO_|KEGG_', '', x)
    x <- gsub('_', ' ', x)
    x <-  str_to_sentence_func(x)
    x <- gsub('CD4 positive alpha beta T-cell', 'CD4 αβ T-cell', x)
    x <- gsub('Activation of the AP-1 family of transcription factors', 
             'Activation of the AP-1 family', x)
},position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), age_group_lbls_func)+
guides(
       color = guide_colorbar(title = "Enrichment\nscore", order = 1),
size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2))+
labs(X='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =45, hjust = 1), 
      #plot.margin=margin(r = 1,unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=14*1/3, repr.plot.height=11)
CD8_Naive_p <- ggplot(CD8_Naive_plot_df, aes(comparison, variable, size = -log10(p_val_adj), color = scale_exp))+
geom_point()+
scale_size_continuous(range = c(0.5,6), 
                     breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_color_gradientn(colours = c("#481055","#32538f", "#00918c", '#00be63', '#fae800'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(CD8_Naive_GS_order), label = function(x){
    x <- gsub('GO_|KEGG_', '', x)
    x <- gsub('_', ' ', x)
    x <-  str_to_sentence_func(x)
    x <- gsub('CD4 positive alpha beta T-cell', 'CD4 αβ T-cell', x)
    x <- gsub('Activation of the AP-1 family of transcription factors', 
             'Activation of the AP-1 family', x)
},position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), age_group_lbls_func)+
guides(
       color = guide_colorbar(title = "Enrichment\nscore", order = 1),
size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2))+
labs(X='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =45, hjust = 1), 
      #plot.margin=margin(r = 1,unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=11, repr.plot.height=7)
p <- wrap_plots(CD4_Naive_label_plot, CD8_Naive_label_plot, CD4_Naive_p, CD8_Naive_p)+
plot_layout(ncol = 2, guide = 'collect', widths = c(1,1,38, 38), heights = c(1,19))&
theme(plot.margin=margin(l = 0,r = 0,unit = 'cm'),legend.background = element_blank(),
      legend.margin = margin(l = 0, unit = 'cm'))
p
#my_saveplot_fun(p, 'Naive_age_genesets', width = 18, height = 8)
ggsave("CD4_8_SP_aging_GS.svg", width = 11, height = 7)

###Extended Data Fig.3e
thymus_Naive_age_gene_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/thymus/Naive/thymus_Naive_age_gene_list.rds')

CD4_Naive_df <- thymus_Naive_age_gene_list[['CD4+Naive']]
CD8_Naive_df <- thymus_Naive_age_gene_list[['CD8+Naive']]

result <- lapply(names(thymus_Naive_age_gene_list), function(x){
    df <- thymus_Naive_age_gene_list[[x]]
    df_iter_result <- lapply(unique(df$comparison), function(y){
        df_iter <- df[df$comparison %in% y & df$p_val_adj < 0.05 & abs(df$log2FC) > 0.5, ]
        df_iter <- df_iter[order(df_iter$p_val_adj,-df_iter$log2FC), ]
        genes <- df_iter$variable[1:20]
        genes <- genes[!is.na(genes)]
        genes
    })
    unique(unlist(df_iter_result))
})
names(result) <- names(thymus_Naive_age_gene_list)


CD4_Naive_selected_genes <- result[['CD4+Naive']]

CD8_Naive_selected_genes <- result[['CD8+Naive']]

CD4_Naive_plot_df <- CD4_Naive_df[CD4_Naive_df$variable %in% CD4_Naive_selected_genes, ]

CD4_Naive_plot_df <- plyr::ddply(CD4_Naive_plot_df, c('variable'), transform, scale_exp = scale(means_1))
CD4_Naive_plot_df$comparison <- gsub('_VS_Rest', '', CD4_Naive_plot_df$comparison)

CD4_Naive_label_df <- data.frame(Age = unique(CD4_Naive_plot_df$comparison))



CD4_Naive_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(CD4_Naive_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))



CD4_Naive_label_plot <- ggplot(CD4_Naive_label_df, aes(Age,1, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'mature CD4 T')

options(repr.plot.width=14*1/3, repr.plot.height=11)



CD8_Naive_plot_df <- CD8_Naive_df[CD8_Naive_df$variable %in% CD8_Naive_selected_genes, ]

CD8_Naive_plot_df <- plyr::ddply(CD8_Naive_plot_df, c('variable'), transform, scale_exp = scale(means_1))
CD8_Naive_plot_df$comparison <- gsub('_VS_Rest', '', CD8_Naive_plot_df$comparison)

CD8_Naive_label_df <- data.frame(Age = unique(CD8_Naive_plot_df$comparison))



CD8_Naive_genes_order <- 
rownames(order_heatmap2(as.matrix(data.frame(
    dcast(CD8_Naive_plot_df, variable~comparison,value.var = 'scale_exp'), check.names=F,row.names = 'variable'))))


CD8_Naive_label_plot <- ggplot(CD8_Naive_label_df, aes(Age,1, fill = Age))+
geom_tile()+
scale_fill_manual(values = age_colors,
                     limits = c('<=12', '13-39', '>=40'), 
                     label = age_group_lbls_func)+
theme(legend.text = element_text(size = 13.75, family='sans'),
      plot.title = element_text(size = 13.75, family='sans',hjust = 0.5),
      plot.margin = margin(t = 0, unit = 'cm'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())+
labs(title = 'mature CD8 T')



combined_p_val_adj <- c(CD4_Naive_plot_df$p_val_adj, CD8_Naive_plot_df$p_val_adj)
combined_p_val_adj[combined_p_val_adj == 0] <- min(combined_p_val_adj[combined_p_val_adj!=0])
CD4_Naive_plot_df$p_val_adj[CD4_Naive_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
CD8_Naive_plot_df$p_val_adj[CD8_Naive_plot_df$p_val_adj==0] <- min(combined_p_val_adj)
combined_p_val_adj <- -log10(combined_p_val_adj)

combined_p_val_adj_break <- 
as.numeric(sprintf("%0.0f",seq(min(combined_p_val_adj), max(combined_p_val_adj), length.out = 5)/10))*10

CD4_Naive_plot_df$p_val_adj <- scales::squish(CD4_Naive_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))
CD8_Naive_plot_df$p_val_adj <- scales::squish(CD8_Naive_plot_df$p_val_adj, range(10^-(combined_p_val_adj_break)))


combined_exp <- c(CD4_Naive_plot_df$scale_exp, CD8_Naive_plot_df$scale_exp)

combined_exp_break <- as.numeric(sprintf("%0.1f",  seq(min(combined_exp), max(combined_exp), length.out = 5)))

CD4_Naive_p <- ggplot(CD4_Naive_plot_df, aes(comparison, variable,# size = -log10(p_val_adj), 
                                             fill = scale_exp))+
geom_tile()+
# scale_size_continuous(range = c(0.5,8), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f16e62'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(CD4_Naive_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(
       fill = guide_colorbar(title = "Expression", order = 1)
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2)
)+
labs(x='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =0, hjust = 0), 
      plot.margin = margin(t = 0, unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans',face = 'italic'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=14*1/3, repr.plot.height=11)
CD8_Naive_p <- ggplot(CD8_Naive_plot_df, aes( comparison,variable,# size = -log10(p_val_adj), 
                                             fill = scale_exp))+
geom_tile()+
# scale_size_continuous(range = c(0.5,8), 
#                      breaks = combined_p_val_adj_break, limits = range(combined_p_val_adj_break))+
scale_fill_gradientn(colours = c('#290124', '#185ea6', '#f16e62'), 
                      oob = scales::squish,
                     breaks = combined_exp_break, limits = range(combined_exp_break))+
scale_y_discrete(limits=rev(CD8_Naive_genes_order), position = 'right')+
scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                 label = age_group_lbls_func)+
guides(
       fill = guide_colorbar(title = "Expression", order = 1)
#size = guide_legend(title = "-log10 P.adj", override.aes = list(fill = 'black'), order = 2)
)+
labs(x='Age', y=NULL)+
theme(axis.text.x = element_blank(),#element_text(size = 13.75, family = 'sans', angle =0, hjust = 0), 
      plot.margin = margin(t = 0, unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 13.75, family = 'sans', face = 'italic'),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())



options(repr.plot.width=6, repr.plot.height=7)
p <- wrap_plots(CD4_Naive_label_plot, CD4_Naive_p,  CD8_Naive_label_plot,CD8_Naive_p, byrow = F)+
plot_layout(ncol = 2, guide = 'collect', heights = c(1,19))&
theme(plot.margin = margin(l = 0.36, unit = 'cm'))
p
ggsave("CD4_CD8_SP_aging_DEGs.svg", width = 6, height = 7)




matureSP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/matureSP.rds')

CD4_matureSP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/CD4_matureSP.rds')

CD8_matureSP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/CD8_matureSP.rds')

lbls_func <- function(x){
    x <- gsub('CD8\\+', 'CD8 ', x)
    x <- gsub('CD4\\+', 'CD4 ', x)
    x <- gsub('Naive', 'SP', x)
}



###Fig.3f
plot_df <- cbind(Embeddings(matureSP, reduction = 'umap'),
                 FetchData(matureSP, vars = 'age_3'))
plot_df$age_group <-age_group_lbls_func(plot_df$age_3)
nbins <- 12
options(repr.plot.width = 3.7*3, repr.plot.height = 3.7)
p <- ggplot(plot_df, aes(UMAP_1, UMAP_2))+
stat_density_2d_filled(bins=nbins)+
d_theme_w(size = 13.75)+
theme(strip.background = element_blank(),axis.text = element_blank(),
      legend.margin = margin(l = -0.1, unit = 'cm'),
      axis.ticks = element_blank())+
geom_point(size = 0.3, stroke = 0, color = 'white')+
facet_wrap(.~age_group)+
labs(fill = 'Density')
lvls <- levels(ggplot_build(p)$data[[1]]$level)
lvls_used <- c(lvls[1], lvls[ceiling(length(lvls)/2)], lvls[length(lvls)])
p <- p+scale_fill_manual(values  = colorRampPalette(c(viridis::viridis(nbins)[1:(nbins-1)]))(nbins), 
                    breaks = lvls_used, label = c('Low', 'Mid', 'High'))
p
my_saveplot_fun(p, 'matureSP_age', width=3.7*3,height=3.7)

###Fig.3g
options(repr.plot.width=4, repr.plot.height=4)
p <- DimPlot(matureSP, reduction = 'umap', label = T)+
guides(color = F)+
scale_color_manual(values = cluster_colors, label = lbls_func)+
d_theme_w(size = 13.75)+
theme(axis.text = element_blank(), axis.ticks = element_blank())

p$layers[[2]]$data$ident <- lbls_func(p$layers[[2]]$data$ident)

p$layers[[2]]$geom_params$size <- 13.75/.pt

p
my_saveplot_fun(p, 'matureSP_idents', width=4,height=4)





###Fig.3h
levels(matureSP) <- 
c('CD4+SOX4+Naive','CD8+SOX4+Naive','CD4+SOX4-Naive','CD8+SOX4-Naive')
exp <- my_group_mean(expm1(matureSP$RNA@data[unique(c('SOX4','SATB1','TOX','TOX2','LCP2','LAT',#'CD28',
                                'STMN1','CCR9','ID3','CD38',
'GZMM','LRRN3','HSPB1','BCL11A','BCL11B','THEMIS','BACH2','IKZF2','LEF1','ID2',
                                'STAT1',
                                
'KLF2', 'S1PR1','CD52','JUNB','JUN','FOS','FOSB','NOSIP','SELL','GPR183','ANXA1','AHNAK','IL6R',
'IL10RA','S100A4','S100A6','FCER1G','CXCR3','KLRK1')), ]), Idents(matureSP))
exp <- my_scale(exp, 'row')
rows_lmls <- rownames(order_heatmap(exp))


exp <- data.frame(exp, check.names = F)
exp$genes <- rownames(exp)
exp <- melt(exp, id.vars = 'genes', variable.name = 'cluster', value.name = 'exp')



options(repr.plot.width=11, repr.plot.height=4)
p <- ggplot(exp, aes( genes, cluster,fill = exp))+
geom_tile()+ 
d_theme_w(size = 13.75)+
theme(plot.margin = margin(l = 0.7, unit = 'cm'),
    axis.text.y =element_text(angle= 0, hjust = 1),legend.key.width = unit(1.5,'cm'),
     axis.text.x =element_text(angle = 45, hjust = 1, face = 'italic'), 
     axis.title = element_blank(), axis.line = element_blank(), 
     panel.background = element_rect(fill = 'white', colour = NA), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     legend.position = 'top', legend.justification='center')+ 
#scale_size_continuous(range=c(2, 7.5)) +coord_flip()+
scale_y_discrete( limits =rev,label = lbls_func, position = 'right')+
scale_x_discrete( limits =c('SOX4','TOX2',
'CD38','STAT1','ID3','BCL11B','BACH2','LAT','SATB1','TOX','CCR9','BCL11A','LEF1',
'HSPB1','THEMIS','GZMM','LCP2','LRRN3','IKZF2','ID2','STMN1','ANXA1','IL6R','KLF2','IL10RA',
                            'AHNAK','S100A6','CD52','S100A4', 'S1PR1',
                            'JUNB','JUN','FOSB','FOS',
                            'NOSIP','SELL','GPR183','CXCR3','FCER1G',
'KLRK1'
))+
#scale_y_discrete(limits = c('T_CD8+SOX4+TN','T_CD8+SOX4-TN','P_CD8+SOX4+TN','P_CD8+SOX4-TN'))+
scale_fill_gradientn(colours = c("lightgray",'#00bade','#f8fc00', '#ff2100'), limits = c(-1.5, 1.5))+
guides( 
       fill = guide_colorbar(title = "Expression"))

ggsave("matureSP_DEGs.svg", width = 11, height = 4)

p





###Fig.3i left
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/may20/CD4_SOX4_pos_Naive_versus_CD4_SOX4_neg_Naive.Gsea.1747710346972//gsea_report_for_CD4_SOX4_pos_Naive_1747710346972.tsv.xlsx')
data1$cluster <- 'CD4+SOX4+Naive'


data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/may20/CD4_SOX4_pos_Naive_versus_CD4_SOX4_neg_Naive.Gsea.1747710346972/gsea_report_for_CD4_SOX4_neg_Naive_1747710346972.tsv.xlsx')
data2$cluster <- 'CD4+SOX4-Naive'


plot_df <- rbind(data1, data2)


pathways <- c(#'GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_UP',

'SIGNALING BY WNT',
#'GO_RAC_PROTEIN_SIGNAL_TRANSDUCTION',
#'GO_NUCLEOSOME_ORGANIZATION',
'KEGG_REGULATION_OF_ACTIN_CYTOSKELETON',
#'GO_STRUCTURAL_CONSTITUENT_OF_CYTOSKELETON',
#'GO_CHROMATIN_REMODELING',
'TCR SIGNALING',
'G2/M TRANSITION',
'GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP',
'GO_DNA_RECOMBINATION',
'GO_T_CELL_DIFFERENTIATION_IN_THYMUS',

#'GO_THYMOCYTE_APOPTOTIC_PROCESS',
#'GO_APOPTOTIC_DNA_FRAGMENTATION',
#'GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_DN',
'EUKARYOTIC TRANSLATION INITIATION',
'GO_LIPOPROTEIN_METABOLIC_PROCESS',
'GO_MHC_PROTEIN_COMPLEX',
#'GO_REGULATION_OF_INCLUSION_BODY_ASSEMBLY',
'INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING',
'INTERLEUKIN-10 SIGNALING',
'GO_T_HELPER_17_TYPE_IMMUNE_RESPONSE',
    'GO_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM',
'GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_DN')

#plot_df$NES <- -1*plot_df$NES
plot_df <- plot_df[abs(plot_df$NES) > 1 & plot_df$NOM.p.val < 0.05 & plot_df$FDR.q.val < 0.25, ]
plot_df <- plot_df[plot_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
plot_df <- plot_df[order(plot_df$NES, decreasing = F), ]
plot_df$GS.br..follow.link.to.MSigDB <- factor(plot_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = plot_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(plot_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color ='black')+
scale_y_discrete(label = function(x){
    x <- gsub('GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_DN', 
             'NAIVE_CD4_TCELL_ADULT_BLOOD_VS_CD4_THYMOCYTE_UP', x)
    x <- gsub('GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_DN', 
             'GSE11057_MEMORY_VS_NAIVE_CD4_TCELL_UP', x)
    x <- gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β', x)
    x <- gsub('KAECH', '', x)
    x <- gsub('GSE\\d*', '', x)
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)

    x <- str_to_sentence_func(x)
    x <- gsub('TUMOR NECROSIS FACTOR', 'TNF', x, ignore.case = T)
    x <- gsub('interferon GAMMA', 'INF-γ', x, ignore.case = T) 
    x <- gsub('tcell', 'Tcell', x)  
    x <- gsub('G2/m', 'G2/M', x)
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), axis.ticks=element_blank(),#legend.box = 'vertical',
#plot.margin = margin(l = 0.2, unit = 'cm'),
      legend.margin = margin(b= -0.2,unit = 'cm'),
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c(cluster_colors), limits= rev,label = lbls_func)+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=6, repr.plot.height=4)
p
ggsave("CD4_SP_GSEA.svg", width = 6, height = 4)

###Fig.3i right
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/may20/CD8_SOX4_pos_Naive_versus_CD8_SOX4_neg_Naive.Gsea.1747710355157/gsea_report_for_CD8_SOX4_pos_Naive_1747710355157.tsv.xlsx')
data1$cluster <- 'CD8+SOX4+Naive'


data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/Enrichment/may20/CD8_SOX4_pos_Naive_versus_CD8_SOX4_neg_Naive.Gsea.1747710355157/gsea_report_for_CD8_SOX4_neg_Naive_1747710355157.tsv.xlsx')
data2$cluster <- 'CD8+SOX4-Naive'


plot_df <- rbind(data1, data2)


pathways <- c(
'GO_T_CELL_RECEPTOR_COMPLEX',
'SIGNALING BY WNT',
'GO_CHROMATIN_REMODELING',
#'GO_PROTEIN_ACYLATION',
'MITOTIC ANAPHASE',
    #'GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_UP',
    #'GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_DN',
    
#'GO_THYMOCYTE_APOPTOTIC_PROCESS',
'GO_STRUCTURAL_CONSTITUENT_OF_CYTOSKELETON',
#'DNA REPAIR',
'GO_T_CELL_DIFFERENTIATION_IN_THYMUS',
'KAECH_NAIVE_VS_MEMORY_CD8_TCELL_UP',
'PEPTIDE CHAIN ELONGATION',
'GSE32423_MEMORY_VS_NAIVE_CD8_TCELL_UP',
'GO_MHC_PROTEIN_COMPLEX',
'INTERFERON GAMMA SIGNALING',
'GO_PROTEIN_LIPID_COMPLEX',
'INTERLEUKIN-10 SIGNALING',
'GO_T_CELL_MIGRATION',
#'GO_NEGATIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE',
#'INTERFERON ALPHA/BETA SIGNALING',
#'GO_INTERFERON_GAMMA_PRODUCTION',
'GO_INTERLEUKIN_6_PRODUCTION',
'GO_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE')

#plot_df$NES <- -1*plot_df$NES
plot_df <- plot_df[abs(plot_df$NES) > 1 & plot_df$NOM.p.val < 0.05 & plot_df$FDR.q.val < 0.25, ]
plot_df <- plot_df[plot_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
plot_df <- plot_df[order(plot_df$NES, decreasing = F), ]
plot_df$GS.br..follow.link.to.MSigDB <- factor(plot_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = plot_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(plot_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color ='black')+
scale_y_discrete(label = function(x){
    x <- gsub('GSE1460_CD4_THYMOCYTE_VS_NAIVE_CD4_TCELL_ADULT_BLOOD_DN', 
             'NAIVE_CD4_TCELL_ADULT_BLOOD_VS_CD4_THYMOCYTE_UP', x)    
    x <- gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β', x)
    x <- gsub('KAECH', '', x)
    x <- gsub('GSE\\d*', '', x)
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)

    x <- str_to_sentence_func(x)
    x <- gsub('TUMOR NECROSIS FACTOR', 'TNF', x, ignore.case = T)
    x <- gsub('interferon GAMMA', 'INF-γ', x, ignore.case = T)
    x <- gsub('tcell', 'Tcell', x)
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), axis.ticks=element_blank(),#legend.box = 'vertical',
#plot.margin = margin(l = 0.2, unit = 'cm'),
      legend.margin = margin(b= -0.2,unit = 'cm'),
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c(cluster_colors), limits= rev,label = lbls_func)+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=6, repr.plot.height=4)
p
ggsave("CD8_SP_GSEA.svg", width = 6, height = 4)







###Fig.3j left
plot_df <- 
data.frame(apply(table(Idents(CD4_matureSP), CD4_matureSP$age_3), 2, function(x){x/sum(x)}), check.names = F)
plot_df$cluster <- rownames(plot_df)

plot_df <- melt(plot_df, id.vars = 'cluster', value.name = 'pct', variable.name = 'age_3')


plot_df$pct <- plot_df$pct*100


options(repr.plot.width=3, repr.plot.height=4)
ggplot(plot_df, aes(age_3, pct, fill  =cluster))+
geom_bar(stat = 'identity', width = 0.5)+
#geom_col(width= 0.5,color = "white",size = 0.1)+

d_theme_w(size = 13.75)+
guides(fill = guide_legend(nrow = 2))+
labs(x = NULL, y = 'Proportion %', fill = 'Identity')+
scale_x_discrete(label = age_group_lbls_func, expand = c(0, 0))+
scale_fill_manual(values = cluster_colors, limits =rev,label = lbls_func)+
scale_y_continuous(expand = c(0,0.1))+
theme(plot.margin = margin(r = 0.4,unit = 'cm'),legend.position = 'top',
           panel.border = element_blank(),#element_rect(color = 'black', fill = NA), 
      panel.background = element_blank())+#element_rect(color = NA, fill = 'white'),)+
  ggalluvial::geom_flow(aes(alluvium = cluster), alpha= .3, color = "white",
            curve_type = "linear", 
            width = .5)
ggsave("CD4_SP_pct.svg", width = 3, height = 4)

###Fig.3j right
plot_df <- 
data.frame(apply(table(Idents(CD8_matureSP), CD8_matureSP$age_3), 2, function(x){x/sum(x)}), check.names = F)
plot_df$cluster <- rownames(plot_df)

plot_df <- melt(plot_df, id.vars = 'cluster', value.name = 'pct', variable.name = 'age_3')


plot_df$pct <- plot_df$pct*100


options(repr.plot.width=3, repr.plot.height=4)
ggplot(plot_df, aes(age_3, pct, fill  =cluster))+
geom_bar(stat = 'identity', width = 0.5)+
#geom_col(width= 0.5,color = "white",size = 0.1)+

d_theme_w(size = 13.75)+
guides(fill = guide_legend(nrow = 2))+
labs(x = NULL, y = 'Proportion %', fill = 'Identity')+
scale_x_discrete(label = age_group_lbls_func, expand = c(0, 0))+
scale_fill_manual(values = cluster_colors, limits =rev,label = lbls_func)+
scale_y_continuous(expand = c(0,0.1))+
theme(plot.margin = margin(r = 0.4,unit = 'cm'),legend.position = 'top',
           panel.border = element_blank(),#element_rect(color = 'black', fill = NA), 
      panel.background = element_blank())+#element_rect(color = NA, fill = 'white'),)+
  ggalluvial::geom_flow(aes(alluvium = cluster), alpha= .3, color = "white",
            curve_type = "linear", 
            width = .5)
ggsave("CD8_SP_pct.svg", width = 3, height = 4)






