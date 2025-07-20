library(ggplot2)
library(Seurat)


library(data.table)

library(ggpubr)

library(openxlsx)

library(patchwork)

library(ComplexHeatmap, lib.loc = .libPaths()[1])
library(aplot)

library(circlize)

library(plyr)
library(reshape2)

library(clusterProfiler)

library(AUCell)
library(SCopeLoomR)
library(SCENIC)

library(Matrix)




plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig5_FigS6'

setwd(plot_dir)

breaks_func <- function(x){
    seq(min(x), max(x), length.out = 3)
}

lbls_func <- function(x){
    x <- gsub('CD27\\+Vδ1T', 'Vδ1T-1', x)
    x <- gsub('CD27lo_Vδ1T', 'Vδ1T-2', x)
    x <- gsub('CD4\\+', 'CD4 ', x)
    x <- gsub('CD8\\+', 'CD8 ', x)
    x <- gsub('_', '-', x)
    return(x)
}

cluster_colors <- 
c(
    'CD8+Naive'='#53A85F', 'CD4+Naive'='#D6E7A3', 
    'thymus_CD8+Naive'='#53A85F','thymus_CD4+Naive'='#D6E7A3','PB_CD8+RTE'='#276d6d',
'PB_CD4+TN'='#7bacd9','PB_CD8+TN'='#DDC4DA',
'PB_CD4+RTE'='#66c5b9','CD8+SOX4+Naive'='#53A85F','CD8+SOX4-Naive' = "#F2C43D", 
                              'CD4+SOX4-Naive' ="#D76483",
                             'CD4+SOX4+Naive'='#D6E7A3')

PB_cluster_colors <-  c('CD8+TEM-B'='#cb50a0','CD8+TEFF'='#cb50a0',
                    'CD8+TEM-K'='#db72b4','CD8+TEM'='#db72b4',
                    'CD8+RTE'='#276d6d','CD4+TN'='#7bacd9',
                        'CD4+CTL'='#7c75ce',#
                    'CD56_Dim_NK'='#c8af1e',#
                    'CD56_Bright_NK'='#d99333',#
                    'CD27lo_Vδ1T'='#93c88e','CD8+TN'='#DDC4DA',
                    'CD4+TEM'='#868ae3',#'#909ee7',##'
                    'Treg'='#6e7cad','CD4+TCM'='#80a2e6',#'#9cb9e5',
                    'MAIT'='#b1b772',
                    'CD4+RTE'='#66c5b9',
                        'Vδ2Vγ9T'='#A0D7C9',
                    'NKT'='#e3cd82',#
                    'CD27+Vδ1T'='#9cb39b',#
                    'IFN_T'='#ec807c',
                    'T_un'='#747474','Tex'='#A1B4BB', 'CD4+TEFF' = '#55b3db',##
                    'CD8+TCM'='#E59CC4', 

 'CD8+RTE'='#276d6d','CD4+TN'='#7bacd9',
                       'CD8+TN'='#DDC4DA', 'CD4+RTE'='#66c5b9')

cluster_colors <- c(cluster_colors, PB_cluster_colors)

module_colors <- c('1' = '#be6b63', '2' = '#86a866', 
                            '3' = '#bda25b', '4' = '#6a82aa', '5' ='#6e578c','6' ='#e5cca9', '7' ='#7cc3cd')

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

PB_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/combined_result/PB_step_2.rds')

PB_step_2$anno_1 <- factor(PB_step_2$anno_1, levels = c('HSPC','MK','cDC',
                      'Mono1','Mono2', 'Mono3', 'pDC','Naive_B','Memory_B','ABC','Plasma','Tcell','NK'))
Idents(PB_step_2) <- PB_step_2$anno_1

######Extended Data Fig.6a
plot_df <- cbind(FetchData(PB_step_2,vars = c('anno_1')), Embeddings(PB_step_2, reduction = 'umap'))

p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = anno_1))+
geom_point(size=0.5)+
scale_color_manual(values = c('NK'='#e0d79e', 'Tcell'='#a3dccb', 'Plasma'= '#847525',
                             'Memory_B'='#809f8c','Naive_B'='#73BDC6','ABC'='#C2B5A2',
                               'Mono2'='#dec5ab','HSPC'='#ffb16b','MK'='#ed9594', 'Mono3' = '#cba266',
                              'Mono1'='#7c6368','cDC'='#B2C2FF','pDC'='#2762A8'), 
                  label = function(x){
                      x <- gsub('_', '-', x)
                  })

new_label <- as.character(1:nlevels(PB_step_2))
names(new_label) <- levels(PB_step_2)


p <- add_label(p, use_new_label = T, new_label = c(new_label),add_circle_to_label = F,
          location_adj= list('10'=c(0.01, -0.02)
                             ),
         label_arg = list(color='black', size=13.75/.pt))+
  guides(color=guide_legend(override.aes = list(size=2, label = c(new_label)), ncol = 4))+
theme(axis.title = element_text(size=13.75),axis.ticks = element_blank(), 
     axis.text = element_blank(), axis.line  = element_blank(), 
       panel.background = element_rect(color = NA, fill = 'white'),
      legend.spacing.x = unit(0.1, 'cm'),
     panel.border = element_rect(color = 'black', fill = NA), #legend.key.width = unit(2, "lines"),
   #legend.key.height = unit(2, "lines"),
      legend.title=element_text(size=13.75),
      
      legend.text = element_text(size=13.75), legend.position = 'bottom', legend.justification='center')+
labs(color = 'Identity')
options(repr.plot.width=5, repr.plot.height=6)
p
my_saveplot_fun(p, 'PB_step_2_idents', width = 5, height = 6)



######Extended Data Fig.6b
Idents(PB_step_2) <- PB_step_2$anno_1

features <- rev(unique(c('CD34', 'KIT', 'PF4', 'TUBB1','SPI1', 'CD33', 'ITGAX',
                     'CD1C', 'CLEC10A',  
                     'CLEC7A','CD14', 'FCGR3A','FCGR3B',
                     'IL3RA', 'LILRA4','CLEC4C',
                     'CD19', 'MS4A1', 'IGHM','TCL1A','CD27', 
                     'TBX21', 'TNFRSF17', 'DERL3', 
                      'CD3E', 'CD3D', 'NCR1', 'NCAM1')))
features <- intersect(features, rownames(PB_step_2))
levels(PB_step_2) <- c('HSPC','MK','cDC',
                      'Mono1','Mono2', 'Mono3', 'pDC','Naive_B','Memory_B','ABC','Plasma','Tcell','NK')
exp <- AverageExpression(PB_step_2, assays = 'RNA', slot = 'data', features = features, verbose = F)$RNA


scale_exp <- t(my_scale(exp, 'row'))
features_lvls <- colnames(scale_exp)
identis_lvls <- rownames(scale_exp)
label_df <- data.frame(Identity = rownames(scale_exp))

scale_exp <- data.frame(scale_exp)
scale_exp$Identity <- rownames(scale_exp)


scale_exp <- melt(scale_exp, id.vars = 'Identity', value.name = 'Expression', variable.name = 'features')
scale_exp$features <- factor(scale_exp$features, levels = features_lvls)
scale_exp$Identity <- factor(scale_exp$Identity, levels = identis_lvls)
label_df$Identity <- factor(label_df$Identity, levels = identis_lvls)
p <- ggplot(scale_exp, aes(features, Identity, fill = Expression))+
geom_tile(color='grey60')+
scale_fill_gradientn(colours = c('#00add9', '#040600', '#fdfc00'), limits = c(-1.5,1.5), oob = scales::squish)+
theme( plot.margin = margin(0,0,0,0, 'cm'),legend.box = 'vertical',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),axis.ticks =element_blank(),
    panel.background = element_rect(color = NA, fill = 'white'), 
     #panel.border = element_rect(fill = NA, colour = 'black'),
    axis.title.x=element_text(size=13.75, family='sans'),
    axis.title.y=element_text(size=13.75, angle =90, family='sans'),
    axis.text.x=element_text(size=13.75, face='italic', angle=45, hjust=1,vjust = 1),
    #axis.text.y=element_text(size=13.75, family='sans'),
    axis.text.y=element_blank(),
    legend.text = element_text(size = 13.75, family='sans', angle = 90,vjust = 0.5),
    legend.title = element_text(size = 13.75, family='sans')
)+
labs(y=NULL, x =NULL, fill = 'Expression')+
guides(fill = guide_colorbar(title.vjust = 0.9))+
scale_x_discrete(position = "bottom") 

label_plot <- ggplot(label_df, aes(1, Identity, fill = Identity))+
geom_tile()+
scale_fill_manual(values = c('NK'='#e0d79e', 'Tcell'='#a3dccb', 'Plasma'= '#847525',
                             'Memory_B'='#809f8c','Naive_B'='#73BDC6','ABC'='#C2B5A2','Mono3' = '#cba266',
                             'Mono2'='#dec5ab','HSPC'='#ffb16b','MK'='#ed9594', 
                             'Mono1'='#7c6368','cDC'='#B2C2FF','pDC'='#2762A8'),
                 label = function(x){
                      x <- gsub('_', '-', x)
                  })+
guides(fill = guide_legend(title.vjust = 0.9))+
theme(legend.text = element_text(size = 13.75, family='sans'),
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
     text=element_blank(), axis.ticks=element_blank())
new_p <- wrap_plots(label_plot, p)+plot_layout(nrow = 1, widths = c(0.5,19), guides = 'collect')&
theme(legend.position = 'top')
options(repr.plot.width=9, repr.plot.height=6)
new_p
#my_saveplot_fun(new_p, 'PB_step_2_DEGs',width=9, height=6)
ggsave('PB_step_2_DEGs.svg',width=9, height=6)



######Extended Data Fig.6c
plot_df = plyr::ddply(data.frame(table(Idents(PB_step_2), PB_step_2$age_3)), 
                      'Var2', transform,percentage=Freq/sum(Freq) *100)

colnames(plot_df) <- c('Cluster', 'Age', 'Freq', 'Percentage')
p <- ggplot(plot_df,aes(x=Age,y=Percentage,fill=Cluster,group=Cluster))+
geom_area(colour="black",size=0.1)+
scale_fill_manual(values = c('NK'='#e0d79e', 'Tcell'='#a3dccb', 'Plasma'= '#847525',
                             'Memory_B'='#809f8c','Naive_B'='#73BDC6','ABC'='#C2B5A2',
                               'Mono2'='#dec5ab','HSPC'='#ffb16b','MK'='#ed9594', 'Mono3' = '#cba266',
'Mono1'='#7c6368','cDC'='#B2C2FF',#'DC1'="#b8b637",'DC3'='#326A3D','DC5'="#BD3460",
                               'pDC'='#2762A8'),
                 label = function(x){
                      x <- gsub('_', '-', x)
                  })+
#ggtitle("Proportion of thymocyte during aging ")+
theme(panel.background = element_blank(),panel.border = element_blank(),
      plot.margin = margin(t = 0.3,l = 0.5,unit = 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),
    axis.title.x=element_text(size=13.75),axis.title.y=element_text(size=13.75),
        axis.text.x=element_text(size=13.75, family='sans', angle =30, hjust=1, vjust = 1),,
      axis.text.y=element_text(size=13.75),
    legend.title = element_text(size = 13.75),
    legend.text = element_text(size = 13.75)
)+
scale_x_discrete(labels = age_group_lbls_func, expand = c(0, 0),
                 limits = c('<=12', '13-39', '40-99',  '>=100'))+
scale_y_continuous(expand = c(0, 0))+
labs(y = 'Percentage %', x=  NULL, fill = 'Identity')
#

write.csv(reshape2::dcast(plot_df, Cluster~Age, value.var = 'Percentage'), file = 'PB_All_cells_Proportion.csv')
options(repr.plot.width=4, repr.plot.height=6)

p
#my_saveplot_fun(p, 'PB_step_2_proportion', width = 4, height = 6)
ggsave( 'PB_step_2_proportion.svg', width = 4, height = 6)



###combined analysis of our data and CITE-seq data
combined_PB <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE164378/combined_PB.rds')

###project a represents CITE-seq data objtained from GEO with accession GSE164378
sub_T1 <- combined_PB[, combined_PB$project=='a']
DefaultAssay(sub_T1) <- 'ADT'
sub_T1 <- NormalizeData(sub_T1,normalization.method = 'CLR')
DefaultAssay(sub_T1) <- 'RNA'

###project b represents our data
sub_T2 <- combined_PB[, combined_PB$project == 'b']


levels(sub_T2) <- c('CD4+RTE', 'CD4+TN', 'CD4+TCM', 'CD4+TEFF',  'CD4+TEM', 'CD4+CTL',
                    'Treg','CD8+RTE','CD8+TN', 'CD8+TCM',
            'CD8+TEM-K', 
           'CD8+TEM-B',
           'MAIT', 'Vδ2Vγ9T', 'CD27+Vδ1T',
           'CD27lo_Vδ1T', 'Tex', 'IFN_T','T_un')

sub_T2$Identity <- Idents(sub_T2)
new_label <- as.character(1:nlevels(sub_T2))
names(new_label) <- levels(sub_T2)
idents_plot_df <- cbind(FetchData(sub_T2, vars = 'Identity'), 
                 Embeddings(sub_T2, reduction = 'umap'))

###Fig.5a


options(repr.plot.width=4.25, repr.plot.height=6.5)
g <- ggplot(idents_plot_df, aes(UMAP_1, UMAP_2, color = Identity))+
geom_point(size= 0.4)+
guides(color = guide_legend(override.aes = list(size = 7)))+
theme(legend.background = element_blank(),plot.margin = margin(r = 0, unit = 'cm'),
      legend.direction = 'vertical',
      legend.key=element_rect(fill='white'),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.ticks=element_blank(),legend.title=element_text(size=13.75, family= 'sans'),
      axis.text=element_blank(),axis.line=element_blank(), axis.title=element_text(size = 13.75),
      legend.text = element_text(size=13.75, family= 'sans'))+
scale_color_manual(values = cluster_colors, label = lbls_func)

g <- add_label(g, use_new_label = T, new_label = c(new_label),
          location_adj= list('4'=c(0.01, 0.035), '10' = c(-0.035, 0), '5' = c(0.02,0),
                            '6' = c(0.02,0.02)),
         label_arg = list(color='black', size=13.75/.pt))+
  guides(color=guide_legend(override.aes = list(size=2, label = c(new_label), alpha=1), ncol=3, 
                           title.hjust = 0.5,
                            title.position = 'top'))+
theme(axis.ticks = element_blank(), 
     axis.text = element_blank(), axis.line  = element_blank(),  legend.key.width = unit(1, "lines"),
   legend.key.height = unit(1.35, "lines"),
      legend.text = element_text(size=13.75), legend.position = 'bottom', 
      legend.justification='center')
g

my_saveplot_fun(g, 'combined_PB_idents',width=4.25, height=6.5)


###Fig.5b
DefaultAssay(sub_T1) <- 'ADT'
plot_df <- cbind(FetchData(sub_T1, vars = c('CD45RA', 'CD45RO')), 
                 Embeddings(sub_T1, reduction = 'umap'))

plot_list <- lapply(c('CD45RA', 'CD45RO'
                     ), function(feature){
    val <- plot_df[, feature]
    if(feature == 'CD45RO'){
        max_val <- quantile(val, 0.5)
        min_val <- quantile(val, 0.2)
    }else{
        max_val <- quantile(val, 0.8)
        min_val <- quantile(val, 0.5)        
    }
    
    
    val[val > max_val] <- max_val
    val[val < min_val] <- min_val
    plot_df[, paste0(feature)] <- val
    p <- ggplot(plot_df, aes_string('UMAP_1', 'UMAP_2', color =paste0(feature)))+
    geom_point(size = 0.25)+
    scale_color_gradientn(colours = c('#d1cdd5', '#2100ff'),  breaks=breaks_func, 
                      label = function(x){sprintf(paste0('%0.', 1, 'f'), x)})+
    theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.title=element_blank(),
        legend.text=element_text(size =13.75, angle=90, vjust=0, hjust=0,family = 'sans',color = 'black'),
      axis.ticks=element_blank(),axis.title=element_blank(),#element_text(size = 15),
      axis.text=element_blank(),axis.line=element_blank(), legend.key.width = unit(5, 'mm') ,
           legend.position = c(0.75,0.795),
        legend.title=element_blank(),legend.background =element_blank(),
       legend.direction = 'horizontal')+
    labs(title = feature)+
    annotation_custom(
        grid::textGrob(label = feature,hjust = 0.5,vjust=0.5,
                       x=grid::unit(0.77,"npc") , 
                       y=grid::unit(0.925,"npc"), 
                       gp=gpar(col = 'black', fontsize = 13.75)))
    return(p)
})



options(repr.plot.width=8, repr.plot.height=4)
p <- wrap_plots(plot_list)+plot_layout(nrow = 1)
p
my_saveplot_fun(p, 'CD45RA_RO', width = 8, height = 4)

###Extended Data Fig.6d
p <- DotPlot(sub_T2,
        features = rev(c(#'CD3E', 
                         'CD4', 'CD40LG',#'ZBTB7B',
                         'CD8A', 'CD8B','RUNX3','SOX4' ,'TOX',
                    'PECAM1', 'CCR7', 'LEF1', 'TCF7','SELL','PASK','ITGB1','ANXA1', 
                      'GATA3', 
                     'CCR6','RORC', 'GZMK', 'CXCR3', 'CXCR6', 
                     'CCL5','NKG7','CST7',#'IFNG', 
                         'TBX21','EOMES','GZMB','GZMA', 'GZMH', 'FCGR3A','PRF1','FOXP3', 'IL2RA',
            'SLC4A10','ZBTB16',
                     'TRDC', 'TRGC1','TRDV2', 'TRGV9','TRGC2','TRDV1','ZNF683',
                     'CD27','CD28',
                     'HAVCR2', 'ISG15','RPL36A')))


g <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black', key_glyph=draw_key_point)+ 
theme(axis.text.x =element_text(angle = 45, hjust = 1, size = 13.75, face = 'italic', vjust=1),
     axis.text.y =element_text( size = 13.75), 
     axis.title = element_blank(), axis.line = element_blank(), 
     panel.background = element_rect(fill = 'white', colour = NA), 
     panel.border = element_rect(fill = NA, colour = 'black'),legend.margin = margin(l = 1.5, r = 1, unit = 'cm'),
     legend.position = 'top', legend.justification='center', legend.text = element_text( size = 15),
     legend.title = element_text( size = 13.75))+ 
scale_size(range=c(0, 6)) +
  
scale_y_discrete(limits = rev(c('CD4+RTE', 'CD4+TN', 'CD4+TCM', 
           'CD4+TEFF',  'CD4+TEM', 
           'CD4+CTL','Treg','CD8+RTE','CD8+TN', 'CD8+TCM','CD8+TEM-K',  'CD8+TEM-B',
           'MAIT', 'Vδ2Vγ9T', 'CD27+Vδ1T',
           'CD27lo_Vδ1T', 'Tex', 'IFN_T','T_un')), label = lbls_func)+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'),order = 2, 
                           label.position = 'bottom', title.vjust = 0.9), 
       fill = guide_colorbar(title = "Expression", title.vjust = 0.9, order = 1))
g
options(repr.plot.width=17, repr.plot.height=5.5)
#my_saveplot_fun(g, 'Tcell_DEG', width = 9, height = 6.5)
ggsave('Tcell_DEG.svg', width = 17, height = 5.5)

###Fig.5c

plot_df = plyr::ddply(data.frame(table(Idents(sub_T2), sub_T2$age_3)), 
                      'Var2', transform,percentage=Freq/sum(Freq) *100)

colnames(plot_df) <- c('Cluster', 'Age', 'Freq', 'Percentage')
p <- ggplot(plot_df,aes(x=Age,y=Percentage,fill=Cluster,group=Cluster))+
geom_area(colour="black",size=0.1)+
scale_fill_manual(values = cluster_colors, label = lbls_func)+
theme(panel.background = element_blank(), panel.border = element_blank(),
      plot.margin = margin(l = 0.5, t = 0.5,unit = 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),
     plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 18,face = "italic"),
    axis.title.x=element_text(size=13.75),axis.title.y=element_text(size=13.75),
    axis.text.x=element_text(size=13.75, angle=30, hjust=1,vjust = 1,
                            ),axis.text.y=element_text(size=13.75),
    legend.title = element_text(size = 13.75),
    legend.text = element_text(size = 13.75)
)+
scale_x_discrete(limits = c('<=12', '13-39', '40-99', '>=100'), expand = c(0, 0),
                 labels = age_group_lbls_func)+
scale_y_continuous(expand = c(0, 0))+
guides(fill = guide_legend(ncol =1))+
labs(x = NULL, y  ='Percentage %', fill = 'Identity')

write.csv(reshape2::dcast(plot_df, Cluster~Age, value.var = 'Percentage'), file = 'PB_T_NK_Proportion.csv')
options(repr.plot.width=4.6, repr.plot.height=6.5)
p

#my_saveplot_fun(p, 'PB_T_NK_2.proportion', width = 5, height = 7)
ggsave('PB_T_NK_2.proportion.svg', width =4.6, height = 6.5)









gene_exp <- gene_exp_raw <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/PB//gene_exp_raw.rds')

gene_exp <- my_scale(gene_exp, 'row')

gene_exp_var <- apply(gene_exp_raw, 1, var)

load('/data1/02.private/dengyj/analysis/thymus/Gene/PB/gene_exp_rl_cl.RData')

plot_features <- c(
'CCR7',
    'ID3','GLS','LFNG',

    'ETS1', 'BCL11B', 'JUND',
'KLF2','LEF1',
    'SATB1','SELL',
'TCF7',
    'TGFB1','S100A8','FYN',
'MT-CO1','MT-CO2','MT-ND1','MT-ND2','S100A9','CCL5','CST7','GZMA',
'GZMK','HLA-DMA','HLA-DQA2','HLA-DQB1','HLA-DRB1',
'NKG7','S100A4','CXCR3','CD69','DUSP1','FOS','JUN','JUNB',#'S100A6',
'HIST1H1E','HMGA1','SELPLG','CD6','RELA', 'FOXO1',
'PTEN','DGKZ','GRK6','YBX1', 'CEBPB', 'RAC1')


gene_exp_rl_filtered <- lapply(unique(gene_exp_rl), function(iter){
    gene_exp_rl_iter <- gene_exp_rl[gene_exp_rl==iter]
    gene_exp_var_filtered <- gene_exp_var[names(gene_exp_rl_iter)]
    gene_exp_var_filtered <- gene_exp_var_filtered[order(gene_exp_var_filtered, decreasing = T)]
    result <- names(gene_exp_var_filtered)[1:60]
    result <- union(result, intersect(c(plot_features), 
                                          names(gene_exp_var_filtered)))    
    result <- result[!is.na(result)]
    result
})

pathway_list <- list('1' = c('GO_CHROMATIN_BINDING',
'Transcriptional regulation by RUNX1',
'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
'KEGG_WNT_SIGNALING_PATHWAY',
'PIP3 activates AKT signaling',
'GO_HISTONE_BINDING',
'KEGG_MTOR_SIGNALING_PATHWAY',
'GO_T_CELL_DIFFERENTIATION',
'KEGG_TGF_BETA_SIGNALING_PATHWAY',
'GO_T_CELL_ACTIVATION'),
                     '2' =c(#'GO_MHC_CLASS_II_PROTEIN_COMPLEX',
'GO_MHC_PROTEIN_COMPLEX',
'GO_FATTY_ACID_BINDING',
'Translocation of ZAP-70 to Immunological synapse',
'Phosphorylation of CD3 and TCR zeta chains',
'PD-1 signaling'
                     ),
                            
                              '3' = c('GO_MHC_CLASS_II_PROTEIN_COMPLEX',
#'GO_RESPONSE_TO_INTERFERON_GAMMA',
'Translocation of ZAP-70 to Immunological synapse',
'Phosphorylation of CD3 and TCR zeta chains',
'GO_ER_TO_GOLGI_TRANSPORT_VESICLE_MEMBRANE',
'PD-1 signaling',
'GO_S100_PROTEIN_BINDING',
'Interferon gamma signaling',
#'Interferon Signaling',
'GO_T_CELL_MIGRATION'),
                             '4' = c(#'GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES',
'Activation of the AP-1 family of transcription factors',
'KEGG_MAPK_SIGNALING_PATHWAY',
#'Interleukin-4 and Interleukin-13 signaling',
'Senescence-Associated Secretory Phenotype (SASP)'
                             ) , '5' =c('KEGG_CHEMOKINE_SIGNALING_PATHWAY',
'GO_POSITIVE_REGULATION_OF_CELL_ADHESION',
'GO_HISTONE_ACETYLTRANSFERASE_BINDING',
'GO_UBIQUITIN_LIKE_PROTEIN_LIGASE_BINDING',
'GO_RESPONSE_TO_INTERLEUKIN_6'), 
'6' = c(#'tRNA processing in the mitochondrion',
#'rRNA processing in the mitochondrion',
'KEGG_OXIDATIVE_PHOSPHORYLATION',
'GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
'GO_MITOCHONDRIAL_MEMBRANE_PART',
'GO_MITOCHONDRION_ORGANIZATION'))

options(repr.plot.width=12, repr.plot.height=20)
gene_exp_p <- 
Heatmap(gene_exp, cluster_columns = F, show_row_names = T,heatmap_legend_param=list(title='Expression'),
        column_split = gene_exp_cl, 
        row_split = gene_exp_rl,row_dend_reorder=T,
        row_names_gp = gpar(fontsize = 4*4/5, fontfamily = 'sans', fontface = 'italic'), 
       top_annotation = HeatmapAnnotation(
               df=data.frame(Age = sapply(strsplit(colnames(gene_exp), split = '~'), function(x) x[1]), 
                             Identity = sapply(strsplit(colnames(gene_exp), split = '~'), function(x) x[2])), 
                 col = list(Age = age_colors, Identity = cluster_colors),
               annotation_name_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = 'bold'),
               annotation_legend_param = list(
                 Age = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(gene_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ), 
               Identity = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(gene_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ))
             ), 
         right_annotation = 
   rowAnnotation(foo = anno_block(
    panel_fun = function(index, levels) {
        grid.rect(gp = gpar(fill = NA, col = "black"))
        text <- pathway_list[[levels]]
        #text <- gsub(',', ',\n', text)
        txt = paste(text, collapse = "\n")
        #txt = paste0(txt, "\n", length(index), " rows")
        #txt = paste(index, collapse = ",")
        grid.text(txt, 0.5, 0.5, rot = 0,
            gp = gpar(col = 'black', fontsize = 3))
    },
    width = max_text_width(unlist(pathway_list)) - unit(5, "cm")
))  )

###Fig.5d
features_lvls <- intersect(rownames(gene_exp)[unlist(row_order(gene_exp_p))], unlist(gene_exp_rl_filtered))

gene_exp_filtered <- gene_exp[rownames(gene_exp) %in% unlist(gene_exp_rl_filtered), ]
gene_exp_filtered <- gene_exp_filtered[features_lvls, ]
plot_features <- intersect(rownames(gene_exp_filtered), plot_features)
text_colors <- module_colors[as.character(gene_exp_rl[intersect(features_lvls, plot_features)])]
row_split <- factor(gene_exp_rl[features_lvls], 
                    levels = as.character(sort(unique(gene_exp_rl[features_lvls]))))
p <- Heatmap(gene_exp_filtered,col = colorRamp2(c(-3,-1.5,0,1.5,3), 
                                                c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b')),
             row_title = NULL,
             cluster_columns = F, cluster_rows = F, show_row_names = F,show_column_names = F,show_heatmap_legend=F,
         #row_split = gene_exp_cl, 
        row_split = row_split,
        heatmap_legend_param=list(title='FC', direction = "vertical", 
                                  labels_gp = gpar(fontsize = 13.75),title_gp = gpar(fontsize = 13.75)),
        row_names_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = 'italic'), 
        column_names_centered=T,
       top_annotation = HeatmapAnnotation(show_legend = F,show_annotation_name = T,
                                       
               df=data.frame(Age = sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[1]), 
                             Identity = sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[2])), 
                 col = list(Age = age_colors, Identity = cluster_colors),
                 annotation_name_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = 'bold'),                              
               annotation_legend_param = list(
                 Age = list(
                   #title = "Baaaaaaar",
                   #at = colnames(gene_exp_filtered),
                   title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
                   ncol = 1
                 ), 
               Identity = list(
                   #title = "Baaaaaaar",
                   #at = colnames(gene_exp_filtered),
                   title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
                   ncol = 1
                 ))
             ),              
             right_annotation =  rowAnnotation(foo = anno_block(
                panel_fun = function(index, levels) {
                    grid.rect(gp = gpar(fill = NA, col = "black"))
                    text <- pathway_list[[levels]]
                    text <- gsub('_', ' ', gsub('GO_|KEGG_', '', 
                                                gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β',text)))
                    
                    
                    text <- str_to_sentence_func(text)
                    
                    text <- gsub('Pten', 'PTEN', text)  
                    text <- gsub('NF KAPPAB', 'NF-κB', text, ignore.case = T)
                    text <- gsub('KAPPAB', 'κB', text, ignore.case = T)
                    text <- gsub('TGF beta', 'TGF-β', text, ignore.case = T)
                    text <- gsub('TGF beta', 'TGF-β', text, ignore.case = T)
                    text <- gsub('interferon gamma', 'INF-γ', text, ignore.case = T)
                    text <- gsub('sasp', 'SASP', text)
                    text <- gsub('Mtor', 'mTOR', text)
                    #text <- gsub(',', ',\n', text)
                    txt = paste(text, collapse = "\n")
                    #txt = paste0(txt, "\n", length(index), " rows")
                    #txt = paste(index, collapse = ",")
                    grid.text(txt, 0.03, 0.5, rot = 0,hjust =0,
                        gp = gpar(col = module_colors[levels], lineheight = 0.7,fontsize = 13.75))
                },
                width = max_text_width(unlist(pathway_list)) #- unit(0.1, "cm")
            )), 
            left_annotation = rowAnnotation(' ' = anno_mark(link_width = unit(7.5,'mm'),side = 'left',
                                                                          #at_adj = -0.5,
                                                                          #extra_text_args = list(vjust = 0.5),
            at = which(rownames(gene_exp_filtered) %in% plot_features), 
            labels = plot_features, 
            labels_gp = gpar(fontsize=13.75,col = text_colors,#'black',#text_colors, 
                             fontface='italic')), 
           'foo' = anno_empty(border = FALSE, 
                              width=
                              (max_text_width(plot_features) - unit(2.2,'cm')) )))
options(repr.plot.width=8.7, repr.plot.height=9.5)
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p, heatmap_legend_side="top", annotation_legend_side="right",
                                                      legend_grouping = "original",ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p
ggsave('PB_Tcell_aging_DEGs.svg', width =8.7, height = 9.5)

###Fig.5d legend
Identity <- unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[2]))
lgd1 = Legend(labels =  lbls_func(Identity), title = "Identity", by_row = T,
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 4,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[Identity])))
                          
Age <- age_group_lbls_func(unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[1])))
lgd2 = Legend(labels =  Age, title = "Age", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 4,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(age_colors[Age])))

                     
Module <- sort(unique(gene_exp_rl[features_lvls]))
lgd3 = Legend(labels =  Module, title = "Module", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 4,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(module_colors[as.character(Module)])))     

col_fun =  colorRamp2(c(-3,-1.5,0,1.5,3), c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))
lgd4 = Legend(col_fun = col_fun, title = "Fold change", direction = 'horizontal',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
             )                                          

pd = packLegend(lgd1,lgd2,lgd3,
                lgd4, direction = "horizontal")

options(repr.plot.width=8, repr.plot.height=1.1)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(0, "cm"), 
                                              y = unit(0, "cm"), just = c("left", "bottom"))))
ggsave('PB_Tcell_aging_DEGs_legend.svg', width = 8, height = 1.1) 



load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/TN_aging_clock_plot_df.RData')

load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/cohorts.RData')

load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/250404_final_model.RData')

healthy_TN_exp_mean <- t(
    readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/TN_exp_mean_orig.rds'))

cor_residuals_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/cor_residuals_df.rds')


color_used <- c('#E5D2DD', '#53A85F', '#F1BB72',  '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#F3B1A0',
'#968175',"#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
'#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
    '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
    '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
    '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
    '#968175'                
)

###Fig.5f left
options(repr.plot.width=4, repr.plot.height=4)
ggplot(train_plot_df, 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5))+
labs(x = 'Actual age', y='Predicted age', title = 'Training cohort')
saveRDS(train_plot_df,file = 'TN_aging_clock_train_scatter.rds')
ggsave('TN_aging_clock_train_scatter.svg', width = 4, height = 4)

###Fig.5f right
options(repr.plot.width=4, repr.plot.height=4)
ggplot(healthy_ind_plot_df, 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5),
     plot.margin = margin(r = 2.4, unit = 'mm'))+
labs(x = 'Actual age', y='Predicted age', title = 'Healthy cohort')
saveRDS(healthy_ind_plot_df,file = 'TN_aging_clock_test_scatter.rds')
ggsave('TN_aging_clock_test_scatter.svg', width = 4, height = 4)

###Fig.5g
options(repr.plot.width=2.5, repr.plot.height=4)
tmp_df <- healthy_ind_plot_df
tmp_df$batch <- rownames(tmp_df)
tmp_df$source <- NULL
tmp_df$age_diff <- NULL
tmp_df <- melt(tmp_df, 
               id.vars = 'batch', variable.name = 'class' , value.name = 'Age')
tmp_df$class <- factor(tmp_df$class, levels=c('true','pred', 'pred_adj'))
ggplot(tmp_df[!tmp_df$class %in% c('pred'),], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0.05, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'Healthy cohort')+
theme(plot.title = element_text(hjust= 0.5))
saveRDS(tmp_df[!tmp_df$class %in% c('pred'),], file= 'TN_healthy_boxplot.rds')
ggsave('TN_healthy_boxplot.svg', width = 2.5, height = 4)

###Fig.5h
options(repr.plot.width=2.5, repr.plot.height=4)
tmp_df <- ind_plot_df[rownames(ind_plot_df) %in% c(centenarians_cohort),]
tmp_df$pred_adj <- tmp_df$pred - predict(loess_m, tmp_df$true)
tmp_df$batch <- rownames(tmp_df)
tmp_df$source <- NULL
tmp_df <- melt(tmp_df, 
               id.vars = 'batch', variable.name = 'class' , value.name = 'Age')
tmp_df$class <- factor(tmp_df$class, levels=c('true','pred', 'pred_adj'))
ggplot(tmp_df[!tmp_df$class %in% c('pred'),], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0, size = 13.75/.pt, hjust = 0.25)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'Centenarian cohort')+
theme(plot.title = element_text(hjust= 0.5))+
scale_y_continuous(limits = c(70,112))
saveRDS(tmp_df[!tmp_df$class %in% c('pred'),], file= 'TN_centenarians_boxplot.rds')
ggsave('TN_centenarians_boxplot.svg', width = 2.5, height = 4)

###Fig.5i
plot_features <-  c('LAT','FLI1','CD74','COX5A','FAS',  'UBA52','RPL38','RPS15A','GCSH' ,  'SLC16A10'
                  )
plot_df <- ind_plot_df[rownames(ind_plot_df) %in% centenarians_cohort,]
tmp_exp <- 
healthy_TN_exp_mean[rownames(plot_df), as.character(cor_residuals_df$ensembl)[
                                 cor_residuals_df$symbol %in% plot_features]]
colnames(tmp_exp) <- cor_residuals_df$symbol[match(colnames(tmp_exp), cor_residuals_df$ensembl)]
tmp_exp <- my_scale(tmp_exp, do_scale = 'col')

plot_df <- cbind(plot_df, tmp_exp)

plot_df$pred_adj <- plot_df$pred - predict(loess_m, plot_df$true)
plot_df$batch <- rownames(plot_df)
plot_df <- melt(plot_df, id.vars = c('pred', 'true', 'source', 'pred_adj', 'batch'), 
                variable.name = 'features', value.name = 'exp')
plot_df$AgeDiff <- plot_df$pred_adj - plot_df$true
plot_df$class<- cor_residuals_df$class[match(plot_df$features, cor_residuals_df$symbol)]
plot_df$class <- gsub('no', 'Neg Correlated', plot_df$class)
plot_df$class <- gsub('yes', 'Pos Correlated', plot_df$class)
plot_df$features <- factor(plot_df$features, levels = plot_features)
# ggplot(plot_df[plot_df$batch %in% centenarians_cohort,], aes())+
# geom_point()
tmp_color_used <- color_used[1:length(unique(plot_features))]
names(tmp_color_used) <- plot_features

plot_list <-lapply(c('Neg Correlated','Pos Correlated') , function(i){
    x_loc <- if(i == 'Pos Correlated'){
        0.05
    }else{
        0.5
    }
    tmp_plot_df <- plot_df[plot_df$class==i, ]
    p <- ggplot(tmp_plot_df, aes(pred_adj-true, exp, color = features))+
    stat_smooth(method = 'lm', se = T, alpha = 0.1)+
    geom_point()+
    #facet_wrap(.~class, scales = 'free_y')+

    scale_color_manual(values = tmp_color_used)+
    d_theme_w(size = 13.75, italic_text = '')+
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.text = element_text(face='italic'))+
    labs(y = 'Scaled Expression', x = 'AgeDiff', color =NULL, title = i)+
        ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..),#fontface='bold',
                         show.legend = FALSE, label.x.npc = x_loc)
    p
})
options(repr.plot.width=6.1, repr.plot.height=4)
wrap_plots(plot_list)+plot_layout(nrow =1, guides = 'collect' )
ggsave('TN_centenarians_cor.svg', width = 6.1, height = 4)
saveRDS(plot_df, file  = 'TN_centenarians_cor.rds')

###Fig.5j
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/apr15/centenarian_pcc.GseaPreranked.1744703898917/gsea_report_for_na_pos_1744703898917.tsv.xlsx')
data1$cluster <- 'pos'


data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/apr15/centenarian_pcc.GseaPreranked.1744703898917/gsea_report_for_na_neg_1744703898917.tsv.xlsx')
data2$cluster <- 'neg'


plot_df <- rbind(data1, data2)


pathways <- c('GO_MHC_PROTEIN_BINDING',
'TCR SIGNALING',
'GO_RESPONSE_TO_INTERFERON_GAMMA',
'DNA REPLICATION',
'GO_INTERLEUKIN_2_PRODUCTION',
'KEGG_GLYCEROPHOSPHOLIPID_METABOLISM',
'INNATE IMMUNE SYSTEM',
'GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY',
#'GO_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
'TRANSCRIPTIONAL REGULATION BY RUNX3',
'GO_T_CELL_MEDIATED_IMMUNITY',
'GO_RESPONSE_TO_INTERLEUKIN_1',
'KEGG_OXIDATIVE_PHOSPHORYLATION',              
'GO_T_CELL_ACTIVATION',
'KEGG_RIBOSOME',
'PEPTIDE CHAIN ELONGATION',
'GO_IMMUNOGLOBULIN_RECEPTOR_BINDING',
'GO_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_BINDING',
'GO_COMPLEMENT_ACTIVATION')

#plot_df$NES <- -1*plot_df$NES
plot_df <- plot_df[plot_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
plot_df <- plot_df[order(plot_df$NES, decreasing = T), ]
plot_df$GS.br..follow.link.to.MSigDB <- factor(plot_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = plot_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(plot_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color ='black')+
scale_y_discrete(label = function(x){
   
    x <- gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β', x)
    
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)

    x <- str_to_sentence_func(x)
    x <- gsub('TUMOR NECROSIS FACTOR', 'TNF', x, ignore.case = T)
    x <- gsub('interferon GAMMA', 'INF-γ', x, ignore.case = T)     
    
},position = 'right')+
theme(plot.margin = margin(l = 0, unit = 'cm'),legend.margin = margin(b= -0.2,unit = 'cm'),
      text = element_text(size = 13.75, family = 'sans'), axis.ticks=element_blank(),#legend.box = 'vertical',
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c('neg'='#009639','pos'='#a074ba'), limits = c('neg','pos'), 
                  label = c('Delay aging','Accelerate aging'))+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=5, repr.plot.height=4)
p
ggsave('TN_centenarians_GSEA_bar.svg', width = 5, height = 4)








HIV_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE243905/HIV_ind_plot_df.rds')
HBV_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE283471/HBV_ind_plot_df.rds')

HIV_ind_plot_df$disease <- 'HIV'

HBV_ind_plot_df$disease <- 'HBV'

MG_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE227835/MG_ind_plot_df.rds')
SLE_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE135779/SLE_ind_plot_df.rds')


MG_ind_plot_df$disease <- 'MG'

SLE_ind_plot_df$disease <- 'SLE'

chronic_infection_df <- rbind(HIV_ind_plot_df, HBV_ind_plot_df)

autoimmune_disease_df <- rbind(MG_ind_plot_df, SLE_ind_plot_df)

load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/chronic_disease.RData')

chronic_infection_cor_residuals_df <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/chronic_infection_cor_residuals_df.rds')

autoimmune_disease_cor_residuals_df <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/autoimmune_disease_cor_residuals_df.rds')

chronic_infection_df_d <- dcast(chronic_infection_df, batch~class, value.var = 'Age')
chronic_infection_residuals <- chronic_infection_df_d$pred_adj-chronic_infection_df_d$true

autoimmune_disease_df_d <- dcast(autoimmune_disease_df, batch~class, value.var = 'Age')
autoimmune_disease_residuals <- autoimmune_disease_df_d$pred_adj-autoimmune_disease_df_d$true





###Fig.5k

options(repr.plot.width=2.5, repr.plot.height=4)
ggplot(HIV_ind_plot_df[!HIV_ind_plot_df$class %in%c('pred') ,], aes(class, Age,fill=class))+

geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'HIV cohort')+
theme(plot.title = element_text(hjust= 0.5))+
scale_y_continuous(limits = c(30,100))
saveRDS(HIV_ind_plot_df[!HIV_ind_plot_df$class %in%c('pred') ,], file = 'TN_HIV_boxplot.rds')
ggsave('TN_HIV_boxplot.svg', width = 2.5, height = 4)

###Fig.5l

options(repr.plot.width=2.5, repr.plot.height=4)
ggplot(HBV_ind_plot_df[!HBV_ind_plot_df$class %in%c('pred'),], aes(class, Age,fill=class))+

geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0.1, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'HBV cohort')+
theme(plot.title = element_text(hjust= 0.5))
saveRDS(HBV_ind_plot_df[!HBV_ind_plot_df$class %in%c('pred') ,], file = 'TN_HBV_boxplot.rds')
ggsave('TN_HBV_boxplot.svg', width = 2.5, height = 4)

###Fig.5m
plot_features <- c('ANXA2','FAS','LAG3', 'CD74','COX5A',   'RGS10','SLC16A10','USP51','CHMP3', 'RPL38'
                  )
plot_df <- chronic_infection_df_d
tmp_exp <- 
TN_exp_mean[plot_df$batch,  plot_features]

tmp_exp <- my_scale(tmp_exp, do_scale = 'col')

plot_df <- cbind(plot_df, tmp_exp)
# plot_df$pred_adj <- plot_df$pred - predict(loess_m, plot_df$true)
plot_df <- melt(plot_df, id.vars = c('pred', 'true', 'batch', 'pred_adj'), 
                variable.name = 'features', value.name = 'exp')
plot_df$AgeDiff <- plot_df$pred_adj - plot_df$true
plot_df$class<- chronic_infection_cor_residuals_df$class[match(plot_df$features, 
                                                               chronic_infection_cor_residuals_df$genes)]
plot_df$class <- gsub('no', 'Neg Correlated', plot_df$class)
plot_df$class <- gsub('yes', 'Pos Correlated', plot_df$class)
# ggplot(plot_df[plot_df$batch %in% centenarians_cohort,], aes())+
# geom_point()
tmp_color_used <- color_used[1:length(unique(plot_features))]
names(tmp_color_used) <- plot_features

plot_list <-lapply(c('Neg Correlated', 'Pos Correlated'), function(i){
    x_loc <- if(i == 'Pos Correlated'){
        0.05
    }else{
        0.48
    }
    tmp_plot_df <- plot_df[plot_df$class==i, ]
    p <- ggplot(tmp_plot_df, aes(pred_adj-true, exp, color = features))+
    stat_smooth(method = 'lm', se = T, alpha = 0.1)+
    geom_point()+
    #facet_wrap(.~class, scales = 'free_y')+

    scale_color_manual(values = tmp_color_used)+
    d_theme_w(size = 13.75, italic_text = '')+
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.text = element_text(face='italic'))+
    labs(y = 'Scaled Expression', x = 'AgeDiff', color =NULL, title = i)+
        ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..),
                         label.y.npc = 1,#fontface='bold',
                         show.legend = FALSE, label.x.npc = x_loc)
    p
})

options(repr.plot.width=6.1, repr.plot.height=4)
wrap_plots(plot_list)+plot_layout(nrow =1, guides = 'collect' )
saveRDS(plot_df, file = 'TN_chronic_infection_cor.rds')
ggsave('TN_chronic_infection_cor.svg', width = 6.1, height = 4)

###Fig.5n
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/may18/chronic_infection_pcc.GseaPreranked.1747551333725/gsea_report_for_na_pos_1747551333725.tsv.xlsx')
data1$cluster <- 'pos'


data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/may18/chronic_infection_pcc.GseaPreranked.1747551333725/gsea_report_for_na_neg_1747551333725.tsv.xlsx')
data2$cluster <- 'neg'


plot_df <- rbind(data1, data2)


pathways <- c('INTERFERON GAMMA SIGNALING',
'TRANSCRIPTIONAL REGULATION BY RUNX3',
'TCR SIGNALING',
'GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY',
'GO_RESPONSE_TO_INTERLEUKIN_1',
'GO_REGULATION_OF_INNATE_IMMUNE_RESPONSE',
'GO_T_CELL_ACTIVATION',
'GO_MHC_PROTEIN_COMPLEX',
'GO_INTERLEUKIN_2_PRODUCTION',
'DNA REPLICATION',
'GO_T_CELL_MEDIATED_CYTOTOXICITY',
'GO_OXIDATIVE_PHOSPHORYLATION',
'KEGG_RIBOSOME',
'PEPTIDE CHAIN ELONGATION',
'EUKARYOTIC TRANSLATION INITIATION',
'GO_CYTOPLASMIC_TRANSLATION')

#plot_df$NES <- -1*plot_df$NES
plot_df <- plot_df[abs(plot_df$NES) > 1 & plot_df$NOM.p.val < 0.05 & plot_df$FDR.q.val < 0.25, ]
plot_df <- plot_df[plot_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
plot_df <- plot_df[order(plot_df$NES, decreasing = T), ]
plot_df$GS.br..follow.link.to.MSigDB <- factor(plot_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = plot_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(plot_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color ='black')+
scale_y_discrete(label = function(x){
    x <- gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β', x)
    
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)

    x <- str_to_sentence_func(x)
    x <- gsub('TUMOR NECROSIS FACTOR', 'TNF', x, ignore.case = T)
    x <- gsub('interferon GAMMA', 'INF-γ', x, ignore.case = T)    
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), axis.ticks=element_blank(),#legend.box = 'vertical',
plot.margin = margin(l = 0.2, unit = 'cm'),legend.margin = margin(b= -0.2,unit = 'cm'),
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c('neg'='#009639','pos'='#a074ba'), limits = c('neg','pos'), 
                  label = c('Delay aging','Accelerate aging'))+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=5, repr.plot.height=4)
p
ggsave('TN_chronic_infection_GSEA_bar.svg', width = 5, height = 4)





###Fig.5o
options(repr.plot.width=2.5, repr.plot.height=4)
ggplot(SLE_ind_plot_df[!SLE_ind_plot_df$class %in%c('pred') ,], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'SLE cohort')+
theme(plot.title = element_text(hjust= 0.5))+

scale_y_continuous(limits = c(5,85))
saveRDS(SLE_ind_plot_df[!SLE_ind_plot_df$class %in%c('pred') ,], file = 'TN_SLE_boxplot.rds')
ggsave('TN_SLE_boxplot.svg', width = 2.5, height = 4)

###Fig.5p
options(repr.plot.width=2.5, repr.plot.height=4)
ggplot(MG_ind_plot_df[!MG_ind_plot_df$class %in%c('pred') ,], aes(class, Age,fill=class))+

geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0.15, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'MG cohort')+
theme(plot.title = element_text(hjust= 0.5))
saveRDS(MG_ind_plot_df[!MG_ind_plot_df$class %in%c('pred') ,], file = 'TN_MG_boxplot.rds')
ggsave('TN_MG_boxplot.svg', width = 2.5, height = 4)





###Fig.5q
plot_features <- c('ANXA2','LAG3','COX5A','CD74','LAT',    'RGS10','SLC16A10','UBA52', 'USP51','GCSH'
                  )
plot_df <- autoimmune_disease_df_d
tmp_exp <- 
TN_exp_mean[autoimmune_disease_df_d$batch,  plot_features]

tmp_exp <- my_scale(tmp_exp, do_scale = 'col')

plot_df <- cbind(plot_df, tmp_exp)
# plot_df$pred_adj <- plot_df$pred - predict(loess_m, plot_df$true)
plot_df <- melt(plot_df, id.vars = c('pred', 'true', 'batch', 'pred_adj'), 
                variable.name = 'features', value.name = 'exp')
plot_df$AgeDiff <- plot_df$pred_adj - plot_df$true

plot_df$class<- autoimmune_disease_cor_residuals_df$class[match(plot_df$features, 
                                                               autoimmune_disease_cor_residuals_df$genes)]
plot_df$class <- gsub('no', 'Neg Correlated', plot_df$class)
plot_df$class <- gsub('yes', 'Pos Correlated', plot_df$class)
# ggplot(plot_df[plot_df$batch %in% centenarians_cohort,], aes())+
# geom_point()
tmp_color_used <- color_used[1:length(unique(plot_features))]
names(tmp_color_used) <- plot_features

plot_list <-lapply(c('Neg Correlated', 'Pos Correlated'), function(i){
    x_loc <- if(i == 'Pos Correlated'){
        0.05
    }else{
        0.48
    }
    tmp_plot_df <- plot_df[plot_df$class==i, ]
    p <- ggplot(tmp_plot_df, aes(pred_adj-true, exp, color = features))+
    stat_smooth(method = 'lm', se = T, alpha = 0.1)+
    geom_point()+
    #facet_wrap(.~class, scales = 'free_y')+

    scale_color_manual(values = tmp_color_used)+
    d_theme_w(size = 13.75, italic_text = '')+
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.text = element_text(face='italic'))+
    labs(y = 'Scaled Expression', x = 'AgeDiff', color =NULL, title = i)+
        ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..),
                         label.y.npc = 1,#fontface='bold',
                         show.legend = FALSE, label.x.npc = x_loc)
    p
})

options(repr.plot.width=6.1, repr.plot.height=4)
wrap_plots(plot_list)+plot_layout(nrow =1, guides = 'collect' )
saveRDS(plot_df, file = 'TN_autoimmune_disease_cor.rds')
ggsave('TN_autoimmune_disease_cor.svg', width = 6.1, height = 4)



###Fig.5r
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/may18/autoimmune_disease_pcc.GseaPreranked.1747551356470/gsea_report_for_na_pos_1747551356470.tsv.xlsx')
data1$cluster <- 'pos'


data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/may18/autoimmune_disease_pcc.GseaPreranked.1747551356470/gsea_report_for_na_neg_1747551356470.tsv.xlsx')
data2$cluster <- 'neg'


plot_df <- rbind(data1, data2)


pathways <- c('INNATE IMMUNE SYSTEM',
'GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR',
'INTERLEUKIN-1 SIGNALING',
'GO_T_CELL_MEDIATED_IMMUNITY',
'GO_INTERLEUKIN_2_PRODUCTION',
'GO_T_CELL_ACTIVATION',
'KEGG_GLYCEROPHOSPHOLIPID_METABOLISM',
'GO_RESPONSE_TO_INTERFERON_GAMMA',
'GO_OXIDATIVE_PHOSPHORYLATION',
'GO_T_CELL_PROLIFERATION',
'TRANSCRIPTIONAL REGULATION BY RUNX3',
'GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
'GO_MHC_PROTEIN_BINDING',
'PEPTIDE CHAIN ELONGATION',
'KEGG_RIBOSOME',
'INITIAL TRIGGERING OF COMPLEMENT',
'EUKARYOTIC TRANSLATION INITIATION',
'GO_IMMUNOGLOBULIN_RECEPTOR_BINDING')

#plot_df$NES <- -1*plot_df$NES
plot_df <- plot_df[abs(plot_df$NES) > 1 & plot_df$NOM.p.val < 0.05 & plot_df$FDR.q.val < 0.25, ]
plot_df <- plot_df[plot_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
plot_df <- plot_df[order(plot_df$NES, decreasing = T), ]
plot_df$GS.br..follow.link.to.MSigDB <- factor(plot_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = plot_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(plot_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color ='black')+
scale_y_discrete(label = function(x){ 
    x <- gsub('TRANSFORMING_GROWTH_FACTOR_BETA', 'TGF-β', x)
    
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)

    x <- str_to_sentence_func(x)
    x <- gsub('TUMOR NECROSIS FACTOR', 'TNF', x, ignore.case = T)
    x <- gsub('interferon GAMMA', 'INF-γ', x, ignore.case = T)      
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), axis.ticks=element_blank(),#legend.box = 'vertical',
      plot.margin = margin(l = 0.2, unit = 'cm'),legend.margin = margin(b= -0.2,unit = 'cm'),
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c('neg'='#009639','pos'='#a074ba'), limits = c('neg','pos'), 
                  label = c('Delay aging','Accelerate aging'))+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=5, repr.plot.height=4)
p
ggsave('TN_autoimmune_disease_GSEA_bar.svg', width = 5, height = 4)




load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/T_aging_clock_plot_df.RData')

###Extended Data Fig.6f top
options(repr.plot.width=5, repr.plot.height=5)
ggplot(train_plot_df, 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5))+
labs(x = 'Actual age', y='Predicted age', title = 'Training cohort')
saveRDS(train_plot_df, file = 'T_aging_clock_train_scatter.rds')
ggsave('T_aging_clock_train_scatter.svg', width = 5, height = 5)

###Extended Data Fig.6f bottom
options(repr.plot.width=5, repr.plot.height=5)
ggplot(healthy_ind_plot_df[ ], 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5),
     plot.margin = margin(r = 2.4, unit = 'mm'))+
labs(x = 'Actual age', y='Predicted age', title = 'Healthy cohort')
saveRDS(healthy_ind_plot_df, file = 'T_aging_clock_test_scatter.rds')
ggsave('T_aging_clock_test_scatter.svg', width = 5, height = 5)










