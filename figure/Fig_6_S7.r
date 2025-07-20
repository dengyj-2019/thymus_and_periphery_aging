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




plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig6_FigS7'
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

age_colors <- c('<=12'="#b7d9ee",'Prepuberal' = '#b7d9ee', '13-39'="#3ca88e",'Adult' = '#3ca88e', 
                '>=40'="#ec8f46", '40-99'="#ec8f46",'Aged' = '#ec8f46', '>=100' = '#c04e5c', 'Centenarian' = '#c04e5c')

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

library(ggraph)
library(tidygraph)

load('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_Tcell_dev.RData')



###Extended Data Fig.7a

PB_CD4T_FDG_df <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD4T_FDG_df.csv', 
                         row.names= 1)


options(repr.plot.width=7.5, repr.plot.height=7.5)
connection <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD4T_paga_con_df.csv', 
                       row.names = 1, check.names = F)
pos <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD4T_paga_pos_df.csv',
                row.names = 1)
colnames(pos) <- c('x', 'y')
connection$Identity <- rownames(connection)
plot_graph <- melt(connection, id.vars = 'Identity', variable.name = 'Identity2', value.name = 'Connectivity')
plot_graph$threshold <- 'no'
plot_graph$threshold[plot_graph$Connectivity > 0.13] <- 'yes' #0.13, 0.11

graph <- as_tbl_graph(plot_graph, directed = F) %>%
  mutate(size = as.data.frame(table(Idents(PB_CD4T)), row.names = 1)[name, ])


options(repr.plot.width=4.5, repr.plot.height=4.5)
p <- ggraph(graph, layout = pos)+ 
  geom_edge_link(aes(width = Connectivity, colour = threshold, linetype = threshold), show.legend = F)+ 
  scale_edge_width_continuous(range = c(0.25, 2)) + 
scale_edge_color_manual(values = c('yes' = '#000000', 'no' = 'lightgrey')) + 
  scale_edge_linetype_manual(values = c('yes' = 'solid','no' = 'dashed')) + 
  geom_node_point(aes(#size = size, 
                      fill = name), shape=21, color = 'black', size = 8) + 
  scale_size(range=c(5, 10))  + scale_fill_manual(values = cluster_colors) + 
scale_color_manual(values = cluster_colors) + 
  theme(panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title=element_text(size=13.75, color ='black'),
       #legend.text = element_text(size=17), 
       #legend.title = element_text(size=17)
       ) + 
labs(x = 'FA_x', y = 'FA_y')+
guides(size = F, fill = F, color = F)+
  labs(colour = 'Identity') +
geom_point(data = PB_CD4T_FDG_df,mapping = aes(x,y, color = identity),size = 0.5)

p$layers <- c(p$layers[[length(p$layers)]], p$layers[1:(length(p$layers)-1)])

p

my_saveplot_fun(p, 'PB_CD4_T_PAGA', width = 4.5,height = 4.5)




###Extended Data Fig.7a legend
CD4_dev_lbls <- c('CD4+RTE', 'CD4+TN', 'CD4+TCM', 'CD4+TEFF', 'CD4+TEM', 'CD4+CTL','Treg')
lgd1 = Legend(labels = lbls_func(CD4_dev_lbls) , title = "Identity", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 19,size = unit(0.2, "npc"),
              ncol = 1,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[ CD4_dev_lbls])))


pd = packLegend(lgd1, direction = "horizontal")


options(repr.plot.width=2.5, repr.plot.height=3)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(1, "cm"), 
                                              y = unit(1, "cm"), just = c("left", "bottom"))))
ggsave('PB_CD4_T_PAGA_legend.svg', width =2.5, height = 3)



###Extended Data Fig.7b
PB_CD8T_FDG_df <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD8T_FDG_df.csv', 
                         row.names= 1)


options(repr.plot.width=7.5, repr.plot.height=7.5)
connection <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD8T_paga_con_df.csv', 
                       row.names = 1, check.names = F)
pos <- read.csv('/data1/02.private/dengyj/analysis/thymus/Development/PB_dev/PB_CD8T_paga_pos_df.csv',
                row.names = 1)
colnames(pos) <- c('x', 'y')
connection$Identity <- rownames(connection)
plot_graph <- melt(connection, id.vars = 'Identity', variable.name = 'Identity2', value.name = 'Connectivity')
plot_graph$threshold <- 'no'
plot_graph$threshold[plot_graph$Connectivity > 0.13] <- 'yes' #0.13, 0.11
graph <- as_tbl_graph(plot_graph, directed = F) %>%
  mutate(size = as.data.frame(table(Idents(PB_CD8T)), row.names = 1)[name, ])


options(repr.plot.width=4.5, repr.plot.height=4.5)
p <- ggraph(graph, layout = pos)+ 
  geom_edge_link(aes(width = Connectivity, colour = threshold, linetype = threshold), show.legend = F)+ 
  scale_edge_width_continuous(range = c(0.25, 2)) + 
scale_edge_color_manual(values = c('yes' = '#000000', 'no' = 'lightgrey')) + 
  scale_edge_linetype_manual(values = c('yes' = 'solid','no' = 'dashed')) + 
  geom_node_point(aes(#size = size, 
                      fill = name), shape=21, color = 'black', size = 8) + 
  scale_size(range=c(5, 10))  + scale_fill_manual(values = cluster_colors) + 
scale_color_manual(values = cluster_colors) + 
  theme(panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title=element_text(size=13.75, color ='black'),
       #legend.text = element_text(size=17), 
       #legend.title = element_text(size=17)
       ) + 
labs(x = 'FA_x', y = 'FA_y')+
guides(size = F, fill = F, color = F)+
  labs(colour = 'Identity') +
geom_point(data = PB_CD8T_FDG_df,mapping = aes(x,y, color = identity),size = 0.5)

p$layers <- c(p$layers[[length(p$layers)]], p$layers[1:(length(p$layers)-1)])

p
my_saveplot_fun(p, 'PB_CD8_T_PAGA', width = 4.5,height = 4.5)




###Extended Data Fig.7b legend
CD8_dev_lbls <- c('CD8+RTE','CD8+TN', 'CD8+TCM', 'CD8+TEM-K', 'CD8+TEM-B')
lgd1 = Legend(labels =  lbls_func(CD8_dev_lbls), title = "Identity", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 16,size = unit(0.2, "npc"),
              ncol = 1,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[CD8_dev_lbls])))


pd = packLegend(lgd1, direction = "horizontal")

options(repr.plot.width=2.5, repr.plot.height=2.5)
ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(1, "cm"), y = unit(1, "cm"), 
                                              just = c("left", "bottom"))))
ggsave('PB_CD8_T_PAGA_legend.svg', width =2.5, height = 2.5)



load('/data1/02.private/dengyj/analysis//thymus/Integrated/thymus_PB/combined_Naive.RData')

###Fig.6a
plot_df <- cbind(FetchData(CD4_combined_Naive, vars = 'identity'), 
                 Embeddings(CD4_combined_Naive, 'pca')[, c(2:3)])

cutoff <- quantile(plot_df$PC_2, probs = c(0.01, 0.99))

plot_df_sub <- plot_df[plot_df$PC_2 > cutoff[1] & plot_df$PC_2 < cutoff[2],  ]
loess.data <- stats::loess(PC_3 ~ PC_2, data = plot_df_sub,
                           span = 0.5)

loess.predict <- predict(loess.data, se = F)
loess.df <- data.frame(fit = loess.predict,# se = loess.predict$se.fit, 
                       PC_2 = plot_df_sub$PC_2, PC_3 = plot_df_sub$PC_3)

options(repr.plot.width=4, repr.plot.height=4)
CD4_Naive_p <- ggplot()+
geom_point(data = plot_df, mapping = aes(PC_2, PC_3, color=identity), size = 0.25)+
scale_color_manual(values = cluster_colors, label = function(x){
    x <- gsub('Naive', 'SP', x)
    x <- gsub('thymus', '', x)
    x <- gsub('PB', '', x)
    x <- gsub('CD.\\+', ' ', x)
}, limits = c('CD4+SOX4+Naive' , 'CD4+RTE', 'CD4+TN'))+
scale_x_reverse()+
theme(panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA), 
      legend.spacing.x = unit(0, 'cm'),legend.key.size = unit(0, 'cm'),
      axis.ticks=element_blank(),
      axis.text=element_blank(),

      axis.title = element_blank(),
      legend.title = element_blank(),legend.key=element_rect(fill =NA),
      legend.background=element_rect(fill =NA),
      legend.text = element_text(size = 13.75, color = 'black'),
      legend.position = c(0.24, 0.12)
     )+
guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_line(data = loess.df, aes(x = PC_2, y = fit), linetype='dashed', size = 0.8,
            arrow = arrow(length=unit(0.50,"cm"), ends="last", type = 'closed')) 


plot_df <- 
cbind(FetchData(CD4_combined_Naive, vars = c('SOX4', 'TOX', 'KLF2',  
                                                'S1PR1', 'SELL'), 
                    slot = 'data'), Embeddings(CD4_combined_Naive, reduction = 'pca'))



plot_list <- lapply(c('SOX4', 'TOX', 'KLF2',   'S1PR1', 'SELL'), function(feature){
    p <- ggplot(plot_df, aes_string('PC_2', 'PC_3', color =feature))+
    geom_point(size = 0.25)+
    scale_color_gradientn(colours = c('lightgrey', 'red'), label = function(x){round(x, 1)},
                          breaks = breaks_func)+
    theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.title=element_blank(),#element_blank(),#element_text(size = 25, face = 'italic', family = 'sans'),
        legend.text=element_text(size =13.75, angle=45,hjust=1,family = 'sans', colour = 'black'),
      axis.ticks=element_blank(),
          axis.title=element_blank(),
          #axis.title=element_text(size = 17.5),
      axis.text=element_blank(),
          axis.line=element_blank(), legend.key.size = unit(0.025, 'npc'),
       legend.position = c(0.20, 0.1), legend.title=element_blank(),legend.background =element_blank(),
       legend.direction = 'horizontal')+
    scale_x_reverse()+
    #guides(color = F)+
    annotation_custom(grid::textGrob(label = feature,hjust = 0,
                                   x=grid::unit(0.75,"npc") ,                                
                                   y=grid::unit(0.92,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75, fontface='italic')))#+
#    scale_x_reverse()
    return(p)
})
plot_list <- c(list(CD4_Naive_p), plot_list)


options(repr.plot.width=9, repr.plot.height=6)
CD4_P <- wrap_plots(plot_list)+plot_layout(ncol = 3)
CD4_P
my_saveplot_fun(CD4_P, 'thymus_PB_CD4_Naive_DEGs', width = 9, height = 6)

###Fig.6b
plot_df <- cbind(FetchData(CD8_combined_Naive, vars = 'identity'), 
                 Embeddings(CD8_combined_Naive, 'pca')[, c(2:3)])

cutoff <- quantile(plot_df$PC_2, probs = c(0.01, 0.99))

plot_df_sub <- plot_df[plot_df$PC_2 > cutoff[1] & plot_df$PC_2 < cutoff[2],  ]
loess.data <- stats::loess(PC_3 ~ PC_2, data = plot_df_sub,
                           span = 0.5)

loess.predict <- predict(loess.data, se = F)
loess.df <- data.frame(fit = loess.predict,# se = loess.predict$se.fit, 
                       PC_2 = plot_df_sub$PC_2, PC_3 = plot_df_sub$PC_3)


options(repr.plot.width=6, repr.plot.height=6)
CD8_Naive_p <- ggplot()+
geom_point(data = plot_df, mapping = aes(PC_2, PC_3, color=identity),
           size = 0.25)+scale_color_manual(values = cluster_colors, label = function(x){
    x <- gsub('Naive', 'SP', x)
    x <- gsub('thymus', '', x)
    x <- gsub('PB', '', x)
    x <- gsub('CD.\\+', ' ', x)
}, limits =  c('CD8+SOX4+Naive' , 'CD8+RTE', 'CD8+TN'))+
theme(panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA), 
      legend.spacing.x = unit(0, 'cm'),legend.key.size = unit(0, 'cm'),
      axis.ticks=element_blank(),
      axis.text=element_blank(),

      axis.title = element_blank(),
      legend.title = element_blank(),legend.key=element_rect(fill =NA),
      legend.background=element_rect(fill =NA),
      legend.text = element_text(size = 13.75, color = 'black'),
      legend.position = c(0.25, 0.15)
     )+
guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_line(data = loess.df, aes(x = PC_2, y = fit), linetype = 'dashed',size = 0.8,
            arrow = arrow(length=unit(0.50,"cm"), ends="last", type = 'closed')) #+
#scale_x_reverse()


plot_df <- 
cbind(FetchData(CD8_combined_Naive, vars = c('SOX4', 'TOX', 
                                               'KLF2', 'S1PR1',  'SELL'), 
                    slot = 'data'), Embeddings(CD8_combined_Naive, reduction = 'pca'))



plot_list <- lapply(c('SOX4', 'TOX', 
                                               'KLF2', 'S1PR1',  'SELL'), function(feature){
    p <- ggplot(plot_df, aes_string('PC_2', 'PC_3', color =feature))+
    geom_point(size = 0.25)+
    scale_color_gradientn(colours = c('lightgrey', 'red'), label = function(x){round(x, 1)},
                      breaks = breaks_func)+
    theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.title=element_blank(),#element_blank(),#element_text(size = 25, face = 'italic', family = 'sans'),
        legend.text=element_text(size =13.75, angle=45,hjust=1,family = 'sans', colour = 'black'),
      axis.ticks=element_blank(),
          axis.title=element_blank(),
          #axis.title=element_text(size = 13.75),
      axis.text=element_blank(),
          axis.line=element_blank(), legend.key.size = unit(0.025, 'npc'),
       legend.position = c(0.27, 0.1), legend.title=element_blank(),legend.background =element_blank(),
       legend.direction = 'horizontal')+
    #guides(color = F)+
    annotation_custom(grid::textGrob(label = feature,hjust = 0,
                                   x=grid::unit(0.75,"npc") ,                                
                                   y=grid::unit(0.92,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75, fontface='italic')))#+
#    scale_x_reverse()
    return(p)
})
plot_list <- c(list(CD8_Naive_p), plot_list)


options(repr.plot.width=9, repr.plot.height=6)
CD8_P <- wrap_plots(plot_list)+plot_layout(ncol = 3)
CD8_P
my_saveplot_fun(CD8_P, 'thymus_PB_CD8_Naive_DEGs', width = 9, height = 6)



###Fig.6c
features<-c('SOX4','ITM2A','LRRN3','TOX','STMN1','SATB1','CD38','PECAM1','ANXA1','GPR183','PASK')

p <- DotPlot(CD4_combined_Naive[, grep('RTE|TN', CD4_combined_Naive$identity)], group.by = 'identity',
             features = features)

g <- ggplot(p$data, aes(id, features.plot,fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_x_discrete(limits =c('CD4+RTE', 'CD4+TN'), label = function(x){
    x <- gsub('CD\\d\\+', ' ', x)
    x <- gsub('PB', '', x)
})+ 
scale_size(range=c(0.5, 6.5))+
theme(axis.title = element_blank(), 
      legend.spacing.y = unit(0.1,'mm'),
      legend.margin = margin(0,0,0,0, unit="mm"),
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.y = element_text(hjust = 1, 
                                 size = 13.75, family = 'sans',  face = 'italic'),
     axis.text.x = element_text(size = 13.75, family = 'sans', face = 'plain',angle =45,hjust=1), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", order=1,override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression")) 
options(repr.plot.width=3.5, repr.plot.height=3.5)
g
ggsave('thymus_PB_CD4_Naive_dotplot_DEGs.svg', width = 3.5, height = 3.5)




###Fig.6d
features<-c('SOX4','FGFBP2','TOX','STMN1','LRRN3','CD9','CD38','MARCKSL1','ANXA1','FXYD5',
            'RIC3', 'SELENOM', 
           'SERPINB6')
options(repr.plot.width=4, repr.plot.height=4.5)
p <- DotPlot(CD8_combined_Naive[, grep('RTE|TN', CD8_combined_Naive$identity)], group.by = 'identity',
             features = features)

g <- ggplot(p$data, aes(id, features.plot,fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black')+ 
scale_x_discrete(limits = c('CD8+RTE', 'CD8+TN'), label = function(x){
    x <- gsub('CD\\d\\+', ' ', x)
    x <- gsub('PB', '', x)    
})+ 
scale_size(range=c(0.5, 6.5))+
theme(axis.title = element_blank(), 
      legend.spacing.y = unit(0.1,'mm'),
      legend.margin = margin(0,0,0,0, unit="mm"),      
      legend.text = element_text(size = 13.75, family = 'sans'), 
      legend.title = element_text(size = 13.75, family = 'sans'), 
      axis.text.y = element_text(hjust = 1, 
                                 size = 13.75, family = 'sans',  face = 'italic'),
     axis.text.x = element_text(size = 13.75, angle =45,hjust = 1,family = 'sans', face = 'plain'), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
scale_fill_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))+
guides(size = guide_legend(title = "% Cells", order=1,override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression")) 
options(repr.plot.width=3.5, repr.plot.height=3.5)
g
ggsave('thymus_PB_CD8_Naive_dotplot_DEGs.svg', width = 3.5, height = 3.5)






matureSP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/matureSP.rds')

Tcell <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell.rds')

matureSP$Identity <- as.character(Idents(matureSP))
Tcell$Identity <- as.character(Idents(Tcell))


matureSP$Identity <- paste0('thymus_', matureSP$Identity)

Tcell$Identity <- paste0('PB_', Tcell$Identity)



CD4_combined_T <- merge(matureSP[,grepl('CD4.*Naive', matureSP$Identity)], 
                       Tcell[, grepl('CD4', Tcell$Identity)])


Idents(CD4_combined_T) <- CD4_combined_T$Identity

CD4_genes_mean <- rowMeans(expm1(CD4_combined_T$RNA@data))

CD4_exp <- AverageExpression(CD4_combined_T, features = names(CD4_genes_mean)[CD4_genes_mean > 0.5], 
                             verbose = F,assays = 'RNA')$RNA

CD4_cor_result <- cor(CD4_exp, method = 'spearman')

CD4_cor_result

###Fig.6e
CD4_cor_result_plot_df <- as.data.frame(as.table(CD4_cor_result))
colnames(CD4_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- c(
    
'thymus_CD4+SOX4+Naive',#'thymus_CD4+SOX4-Naive',
    'PB_CD4+RTE','PB_CD4+TN','PB_CD4+TCM','PB_CD4+TEFF',
'PB_CD4+TEM','PB_CD4+CTL')


p <- ggplot(CD4_cor_result_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),na.value="lightgrey",
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         #limits = c(0.7,1)
                        )+
    theme(#plot.margin = margin(b = 0.9,l = 0.9, unit = 'cm'),
          strip.background = element_blank(),strip.text.x = element_text(size = 13.75),
          axis.text.x = element_text(size = 13.75, family = 'sans', angle =45, hjust = 1), 
          panel.background = element_rect(color = NA, fill ='white'),
          panel.border = element_rect(color = 'black', fill =NA),
          legend.box = 'vertical',
          legend.position = 'right',
          axis.ticks = element_blank(), axis.title = element_blank(), 
         axis.text.y = element_text(size = 13.75, family = 'sans'),
          legend.text = element_text(size = 13.75, family = 'sans'),
          legend.title = element_text(size = 13.75, family = 'sans'),
          plot.title = element_text(size = 13.75, family = 'sans', hjust =0.5),


          axis.line =element_blank())+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=5, repr.plot.height=3.5)
p

ggsave('Spearman_cor_CD4_T.svg', width = 5, height = 3.5)



CD8_combined_T <- merge(matureSP[,grepl('CD8.*Naive', matureSP$Identity)], 
                       Tcell[, grepl('CD8', Tcell$Identity)])


Idents(CD8_combined_T) <- CD8_combined_T$Identity

CD8_genes_mean <- rowMeans(expm1(CD8_combined_T$RNA@data))

CD8_exp <- AverageExpression(CD8_combined_T, features = names(CD8_genes_mean)[CD8_genes_mean > 0.5], 
                             verbose = F,assays = 'RNA')$RNA

CD8_cor_result <- cor(CD8_exp, method = 'spearman')

CD8_cor_result

###Fig.6f
CD8_cor_result_plot_df <- as.data.frame(as.table(CD8_cor_result))
colnames(CD8_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- c(
'thymus_CD8+SOX4+Naive',#'thymus_CD8+SOX4-Naive',
         'PB_CD8+RTE','PB_CD8+TN','PB_CD8+TCM','PB_CD8+TEM-K', 'PB_CD8+TEM-B')



p <- ggplot(CD8_cor_result_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
   #geom_text(mapping = aes(label = round(Cor,2)))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),na.value="lightgrey",
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         #limits = c(0.6,1)
                        )+
    theme(#plot.margin = margin(b = 0.9,l = 0.9, unit = 'cm'),
          strip.background = element_blank(),strip.text.x = element_text(size = 13.75),
          axis.text.x = element_text(size = 13.75, family = 'sans', angle =45, hjust = 1),
          panel.background = element_rect(color = NA, fill ='white'),
          panel.border = element_rect(color = 'black', fill =NA),
          legend.box = 'vertical',
          legend.position = 'right',
          axis.ticks = element_blank(), axis.title = element_blank(), 
         axis.text.y = element_text(size = 13.75, family = 'sans'),
          legend.text = element_text(size = 13.75, family = 'sans'),
          legend.title = element_text(size = 13.75, family = 'sans'),
          plot.title = element_text(size = 13.75, family = 'sans', hjust =0.5),
          #strip.background.y = element_blank(), 
         #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

          axis.line =element_blank())+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=5, repr.plot.height=3.5)
p

ggsave('Spearman_cor_CD8_T.svg', width = 5, height = 3.5)





###Extended Data Fig.7c
CD4_cor_result_plot_df <- as.data.frame(as.table(CD4_cor_result))
colnames(CD4_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- c(
'thymus_CD4+SOX4+Naive','thymus_CD4+SOX4-Naive',
         'PB_CD4+RTE','PB_CD4+TN')



p <- ggplot(CD4_cor_result_plot_df[CD4_cor_result_plot_df$cluster1 %in% cluster_used & 
                                   CD4_cor_result_plot_df$cluster2 %in% cluster_used, ], aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
   #geom_text(mapping = aes(label = round(Cor,2)))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),na.value="lightgrey",
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         breaks = seq(0.85, 1.00, 0.05),
                         #limits = c(0.6,1)
                        )+
    theme(#plot.margin = margin(b = 0.9,l = 0.9, unit = 'cm'),
          strip.background = element_blank(),strip.text.x = element_text(size = 17.75),
          axis.text.x = element_text(size = 17.75, family = 'sans', angle =45, hjust = 1),
          panel.background = element_rect(color = NA, fill ='white'),
          panel.border = element_rect(color = 'black', fill =NA),
          legend.box = 'vertical',
          legend.position = 'right',
          axis.ticks = element_blank(), axis.title = element_blank(), 
         axis.text.y = element_text(size = 17.75, family = 'sans'),
          legend.text = element_text(size = 17.75, family = 'sans'),
          legend.title = element_text(size = 17.75, family = 'sans'),
          plot.title = element_text(size = 17.75, family = 'sans', hjust =0.5),


          axis.line =element_blank())+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=8, repr.plot.height=6.5)
p

ggsave('Spearman_cor_CD4_SP_RTE_TN.svg', width = 8, height = 6.5)

###Extended Data Fig.7d
CD8_cor_result_plot_df <- as.data.frame(as.table(CD8_cor_result))
colnames(CD8_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- c(
'thymus_CD8+SOX4+Naive','thymus_CD8+SOX4-Naive',
         'PB_CD8+RTE','PB_CD8+TN')



p <- ggplot(CD8_cor_result_plot_df[CD8_cor_result_plot_df$cluster1 %in% cluster_used & 
                                   CD8_cor_result_plot_df$cluster2 %in% cluster_used, ], aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
   #geom_text(mapping = aes(label = round(Cor,2)))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),
                         na.value="lightgrey",breaks = seq(0.85, 1.00, 0.05),
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         #limits = c(0.6,1)
                        )+
    theme(#plot.margin = margin(b = 0.9,l = 0.9, unit = 'cm'),
          strip.background = element_blank(),strip.text.x = element_text(size = 17.75),
          axis.text.x = element_text(size = 17.75, family = 'sans', angle =45, hjust = 1),
          panel.background = element_rect(color = NA, fill ='white'),
          panel.border = element_rect(color = 'black', fill =NA),
          legend.box = 'vertical',
          legend.position = 'right',
          axis.ticks = element_blank(), axis.title = element_blank(), 
         axis.text.y = element_text(size = 17.75, family = 'sans'),
          legend.text = element_text(size = 17.75, family = 'sans'),
          legend.title = element_text(size = 17.75, family = 'sans'),
          plot.title = element_text(size = 17.75, family = 'sans', hjust =0.5),


          axis.line =element_blank())+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Naive', 'SP', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=8, repr.plot.height=6.5)
p

ggsave('Spearman_cor_CD8_SP_RTE_TN.svg', width = 8, height = 6.5)





combined_Naive <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/combined_Naive.rds')

###Extended Data Fig.7e


p <- DotPlot(combined_Naive[, grep('CD4', Idents(combined_Naive))], cols = c('lightgrey', 'red'),
        features = unique(c('CCR9','CCR4','CD38',
'TOX2','SOX4','ID3','TOX','SATB1','BCL11A','BACH2','DUSP1',
'CD69','GPR183','FOS','JUN','ZFP36L2',
'ZFP36L1','JUNB','STAT4','STMN1','MARCKSL1','LRRN3',
'IL2RG','SOX4','ETS1','TCF12','TOX2','IKZF2','BCL11B','TOX',
'GIMAP5','ANXA1','S1PR4','CD55','S100A9',
'KLF3','JUND','KLF2', 'S1PR1', 'STAT4', 'LAT', 'S100A4')))


options(repr.plot.width=8, repr.plot.height=10.5)
CD4_dotplot <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black', key_glyph=draw_key_point)+ 
theme(axis.text.y =element_text(size = 17.75, face = 'italic'),
     axis.text.x =element_text(angle = 45, hjust = 1, size = 17.75), 
     axis.title = element_blank(), axis.line = element_blank(), 
     panel.background = element_rect(fill = 'white', colour = NA), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     legend.position = 'right', legend.justification='center', legend.text = element_text( size =17.75),
     legend.title = element_text( size = 17.75))+ 
scale_size_continuous(range=c(0, 7.5), limits = c(0,100),breaks = c(0,25,75,100)) +coord_flip()+
scale_y_discrete(limits = c('CD4+SOX4+Naive','CD4+SOX4-Naive','CD4+RTE','CD4+TN'), label = function(x){
    x <- gsub('Naive', 'SP', gsub('CD.\\+','', x))
})+
scale_fill_gradientn(colours = c('white', 'red'), limits = c(-1.5, 1.5))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression"))
CD4_dotplot
ggsave('CD4_Naive_DEGs.svg', width = 8, height = 10.5)

###Extended Data Fig.7f



p <- DotPlot(combined_Naive[, grep('CD8', Idents(combined_Naive))], cols = c('lightgrey', 'red'),
        features = unique(c('CCR9','CD38',
'TOX2','SOX4','ID3','TOX','SATB1','BCL11A','BACH2','DUSP1',
'CD69','GPR183','FOS','JUN','ZFP36L2',
'ZFP36L1','JUNB','STAT4','STMN1','MARCKSL1','LRRN3',
'IL2RG','SOX4','ETS1','TCF12','TOX2','IKZF2','TOX',
'GIMAP5','ANXA1','S1PR4','CD55','S100A9',
'KLF3','JUND','KLF2', 'S1PR1', 'STAT4', 'LAT', 'S100A4', 'S100A6')))


options(repr.plot.width=8, repr.plot.height=10.5)
CD8_dotplot <- ggplot(p$data, aes(features.plot, id, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21, color = 'black', key_glyph=draw_key_point)+ 
theme(axis.text.y =element_text(size = 17.75, face = 'italic'),
     axis.text.x =element_text(angle = 45, hjust = 1, size = 17.75), 
     axis.title = element_blank(), axis.line = element_blank(), 
     panel.background = element_rect(fill = 'white', colour = NA), 
     panel.border = element_rect(fill = NA, colour = 'black'),
     legend.position = 'right', legend.justification='center', legend.text = element_text( size =17.75),
     legend.title = element_text( size = 17.75))+ 
scale_size_continuous(range=c(0, 7.5), limits = c(0,100),breaks = c(0,25,75,100)) +coord_flip()+
scale_y_discrete(limits = c('CD8+SOX4+Naive','CD8+SOX4-Naive','CD8+RTE','CD8+TN'), label = function(x){
    x <- gsub('Naive', 'SP', gsub('CD.\\+','', x))
})+
scale_fill_gradientn(colours = c('white', 'red'), limits = c(-1.5, 1.5))+
guides(size = guide_legend(title = "% Cells", override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Expression"))
CD8_dotplot
ggsave('CD8_Naive_DEGs.svg', width = 8, height = 10.5)





###combined analysis of our data and CITE-seq data
combined_PB <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE164378/combined_PB.rds')

###project a represents CITE-seq data objtained from GEO with accession GSE164378
sub_T1 <- combined_PB[, combined_PB$project=='a']
DefaultAssay(sub_T1) <- 'ADT'
sub_T1 <- NormalizeData(sub_T1,normalization.method = 'CLR')
DefaultAssay(sub_T1) <- 'RNA'

###project b represents our data
sub_T2 <- combined_PB[, combined_PB$project == 'b']


###Fig.6g
DefaultAssay(sub_T1) <- 'ADT'
plot_df <- cbind(FetchData(sub_T1, vars = c(#'CD45RA', 'CD45RO', 
                                            'CD4-2', 'CD8a', 'CD31','CD38-1','CD38-2', 'SOX4')), 
                 Embeddings(sub_T1, reduction = 'umap'))

colnames(plot_df)[colnames(plot_df) %in% c('rna_SOX4')] <- 'SOX4'

plot_list <- lapply(c(
                      'CD4-2', 'CD8a', 'CD31','CD38-1','CD38-2', 'SOX4'), function(feature){
    if(feature == 'SOX4'){
        colors <- c('lightgrey', 'red')
        fontface='italic'
        
    }else{
        colors <- c('#d1cdd5', '#2100ff')
        fontface='plain'
    }
    plot_df_iter <- plot_df[, c('UMAP_1', 'UMAP_2', feature)]
    min_cutoff <- if(feature == 'SOX4'){
        quantile(plot_df_iter[, feature], 0.05)
    }else{
        quantile(plot_df_iter[, feature], 0.05)
    }
    max_cutoff <- if(feature == 'SOX4'){
        quantile(plot_df_iter[, feature], 0.95)
    }else{
        quantile(plot_df_iter[, feature], 0.925)
    }    
    plot_df_iter[, feature][plot_df_iter[, feature] <= min_cutoff] <- min_cutoff
    plot_df_iter[, feature][plot_df_iter[, feature] >= max_cutoff] <- max_cutoff
    p <- ggplot(plot_df_iter, aes_string('UMAP_1', 'UMAP_2'))+
    geom_point(size = 1, shape = 19, stroke = 0, mapping = aes_string(color =paste0('`', feature, '`')))+ 

    theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        legend.text=element_text(size =13.75,  family = 'sans', color ='black', angle = 45, hjust=1),
      axis.ticks=element_blank(),axis.title=element_blank(),#element_text(size = 15),
      axis.text=element_blank(),axis.line=element_blank(), legend.key.height = unit(0.0125, 'npc'),
          legend.key.width = unit(0.03, 'npc'),
       legend.position = c(0.78, 0.8), legend.title=element_blank(),legend.background =element_blank(),
       legend.direction = 'horizontal')+
    #labs(title = feature)
    annotation_custom(grid::textGrob(label = feature,hjust = 0,
                                   x=grid::unit(0.72,"npc") ,                                
                                   y=grid::unit(0.95,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 13.75, fontface=fontface)))+

    scale_color_gradientn(colours = colors, label = function(x){round(x, 1)},
                      breaks = breaks_func)
    return(p)
})
options(repr.plot.width=6, repr.plot.height=9)
p <- wrap_plots(plot_list)+plot_layout(ncol =2)
p
my_saveplot_fun(p, 'combined_PB_DEGs', width = 6, height = 9)













###Fig.7j top-right
clu_used <- c('CD4+CTL', 'CD4+TEFF', 'CD4+TN','CD4+TEM', 'Treg', 
                         'CD4+TCM',
                        'CD4+RTE')
plot_df <- data.frame(pct=apply(table(Idents(sub_T2), sub_T2$batch)[clu_used,], 2, function(x){
    sum(x[grep('RTE', names(x))])/sum(x) * 100
}))
plot_df$age_group <- sub_T2$age_3[match(rownames(plot_df), sub_T2$batch)]

plot_df <- plot_df[!plot_df$age_group %in% c('>=100'),]

stat_df <- rstatix::pairwise_t_test(plot_df, pct~age_group, pool.sd = F)
stat_df$p.signif <- 'ns'
stat_df$p.signif[stat_df$p < 0.05] <- '*'
stat_df$p.signif[stat_df$p < 0.01] <- "**"
stat_df$p.signif[stat_df$p < 0.001] <- "***"
stat_df$p.signif[stat_df$p < 0.0001] <- "****"
logi <- stat_df$p.signif == 'ns'
stat_df$p.signif[logi] <- stat_df$p[logi]
max_val <- max(plot_df$pct)
len <- nrow(stat_df)
by_val <- 0.2* max_val
stat_df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = len)
p <- ggplot(plot_df, aes(age_group, pct))+

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


d_theme_w(size = 15, italic_text = '')+
theme(panel.border = element_rect(fill = NA, color = NA),axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),plot.margin = margin(t = 0, unit = 'cm'),
      axis.text = element_text(hjust = 0.5,size = 15))+      
guides(color = F,shape = F)+
labs(x = NULL, color = 'Age group', y = '% CD4 RTE of CD4 T-cells')+
geom_point(mapping = aes(color = age_group, shape = age_group))+
scale_shape_manual(values = c('<=12' = 16, '13-39' = 15, '40-99' = 17))+
scale_y_continuous(expand = expansion(mult = .1))
options(repr.plot.width=3, repr.plot.height=3.25)
p

###Fig.7j bottom-right
clu_used <- c('CD8+TEM-K', 'CD8+TCM', 'CD8+TEM-B', 'CD8+TN', 'CD8+RTE')
plot_df <- data.frame(pct=apply(table(Idents(sub_T2), sub_T2$batch)[clu_used,], 2, function(x){
    sum(x[grep('RTE', names(x))])/sum(x) * 100
}))
plot_df$age_group <- sub_T2$age_3[match(rownames(plot_df), sub_T2$batch)]

plot_df <- plot_df[!plot_df$age_group %in% c('>=100'),]

stat_df <- rstatix::pairwise_t_test(plot_df, pct~age_group, pool.sd = F)
stat_df$p.signif <- 'ns'
stat_df$p.signif[stat_df$p < 0.05] <- '*'
stat_df$p.signif[stat_df$p < 0.01] <- "**"
stat_df$p.signif[stat_df$p < 0.001] <- "***"
stat_df$p.signif[stat_df$p < 0.0001] <- "****"
logi <- stat_df$p.signif == 'ns'
stat_df$p.signif[logi] <- stat_df$p[logi]
max_val <- max(plot_df$pct)
len <- nrow(stat_df)
by_val <- 0.2* max_val
stat_df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = len)
p <- ggplot(plot_df, aes(age_group, pct))+

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


d_theme_w(size = 15, italic_text = '')+
theme(panel.border = element_rect(fill = NA, color = NA),axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),plot.margin = margin(t = 0, unit = 'cm'),
      axis.text = element_text(hjust = 0.5,size = 15))+      
guides(color = F,shape = F)+
labs(x = NULL, color = 'Age group', y = '% CD8 RTE of CD8 T-cells')+
geom_point(mapping = aes(color = age_group, shape = age_group))+
scale_shape_manual(values = c('<=12' = 16, '13-39' = 15, '40-99' = 17))+
scale_y_continuous(expand = expansion(mult = .1))
p

###Fig.7j top-left
clu_used <- c('CD4+RTE', 'CD4+TN')
plot_df <- data.frame(pct=apply(table(Idents(sub_T2), sub_T2$batch)[clu_used,], 2, function(x){
    x[grep('RTE', names(x))]/sum(x) * 100
}))
plot_df$age_group <- sub_T2$age_3[match(rownames(plot_df), sub_T2$batch)]

plot_df <- plot_df[!plot_df$age_group %in% c('>=100'),]

stat_df <- rstatix::pairwise_t_test(plot_df, pct~age_group, pool.sd = F)
stat_df$p.signif <- 'ns'
stat_df$p.signif[stat_df$p < 0.05] <- '*'
stat_df$p.signif[stat_df$p < 0.01] <- "**"
stat_df$p.signif[stat_df$p < 0.001] <- "***"
stat_df$p.signif[stat_df$p < 0.0001] <- "****"
logi <- stat_df$p.signif == 'ns'
stat_df$p.signif[logi] <- stat_df$p[logi]
max_val <- max(plot_df$pct)
len <- nrow(stat_df)
by_val <- 0.2* max_val
stat_df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = len)
p <- ggplot(plot_df, aes(age_group, pct))+

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


d_theme_w(size = 15, italic_text = '')+
theme(panel.border = element_rect(fill = NA, color = NA),axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),plot.margin = margin(t = 0, unit = 'cm'),
      axis.text = element_text(hjust = 0.5,size = 15))+      
guides(color = F,shape = F)+
labs(x = NULL, color = 'Age group', y = '% CD4 RTE of CD4 naive T-cells')+
geom_point(mapping = aes(color = age_group, shape = age_group))+
scale_shape_manual(values = c('<=12' = 16, '13-39' = 15, '40-99' = 17))+
scale_y_continuous(expand = expansion(mult = .1))
p



###Fig.7j bottom-left
clu_used <- c('CD8+RTE', 'CD8+TN')
plot_df <- data.frame(pct=apply(table(Idents(sub_T2), sub_T2$batch)[clu_used,], 2, function(x){
    x[grep('RTE', names(x))]/sum(x) * 100
}))
plot_df$age_group <- sub_T2$age_3[match(rownames(plot_df), sub_T2$batch)]

plot_df <- plot_df[!plot_df$age_group %in% c('>=100'),]

stat_df <- rstatix::pairwise_t_test(plot_df, pct~age_group, pool.sd = F)
stat_df$p.signif <- 'ns'
stat_df$p.signif[stat_df$p < 0.05] <- '*'
stat_df$p.signif[stat_df$p < 0.01] <- "**"
stat_df$p.signif[stat_df$p < 0.001] <- "***"
stat_df$p.signif[stat_df$p < 0.0001] <- "****"
logi <- stat_df$p.signif == 'ns'
stat_df$p.signif[logi] <- stat_df$p[logi]
max_val <- max(plot_df$pct)
len <- nrow(stat_df)
by_val <- 0.2* max_val
stat_df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = len)
p <- ggplot(plot_df, aes(age_group, pct))+

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


d_theme_w(size = 15, italic_text = '')+
theme(panel.border = element_rect(fill = NA, color = NA),axis.line = element_line(),
    plot.title = element_text(hjust = 0.5),plot.margin = margin(t = 0, unit = 'cm'),
      axis.text = element_text(hjust = 0.5,size = 15))+      
guides(color = F,shape = F)+
labs(x = NULL, color = 'Age group', y = '% CD8 RTE of naive T-cells')+
geom_point(mapping = aes(color = age_group, shape = age_group))+
scale_shape_manual(values = c('<=12' = 16, '13-39' = 15, '40-99' = 17))+
scale_y_continuous(expand = expansion(mult = .1))
p




