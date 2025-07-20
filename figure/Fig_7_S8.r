library(ggplot2)
library(Seurat)

library(rstatix)
library(data.table)

library(ggpubr)

library(openxlsx)

library(patchwork)

library(plyr)
library(reshape2)

library(treemapify)


library(Matrix)

library(circlize)

library(ComplexHeatmap)
library(future)
library(aplot)

plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig7_FigS8'

setwd(plot_dir)



T_dev_3 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/T_dev/old/250429/T_dev_3.rds')
combined_PB <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE164378/combined_PB.rds')
Tcell <- combined_PB[, combined_PB$project == 'b']

T_dev_3$Identity <- as.character(Idents(T_dev_3))
Tcell$Identity <- as.character(Idents(Tcell))

T_dev_3$Identity[grepl('CD4.*Naive', T_dev_3$Identity)] <- 'CD4+Naive'
T_dev_3$Identity[grepl('CD8.*Naive', T_dev_3$Identity)] <- 'CD8+Naive'

T_dev_3$Identity <- paste0('thymus_', T_dev_3$Identity)

Tcell$Identity <- paste0('PB_', Tcell$Identity)

T_dev_3$age_3 <- ''
T_dev_3$age_3[T_dev_3$age <=12] <- '<=12'
T_dev_3$age_3[T_dev_3$age >= 13 & T_dev_3$age <= 39 ] <- '13-39'
T_dev_3$age_3[T_dev_3$age >=40] <- '40-99'

integer_breaks <- function(n = 5, neg = F, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    if(neg == T){
      breaks <- breaks * -1  
    }      
    names(breaks) <- attr(breaks, "labels")
    breaks
  
  }
  return(fxn)
}


source('/data1/02.private/dengyj/analysis/thymus/TCR/TCR_caculation.R')

cluster_colors <- 
 c('pre_Treg1'='#E5D2DD','CD8+Naive'='#53A85F','CD4+Memory_like'='#808040',
                    'CD8+Memory_like'='#F3B1A0',
                             'CD4+Naive'='#D6E7A3','IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3','Agonist_2'='#a83117','Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'Agonist_1' = '#a6c050', 'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'pre_Treg2' = '#9566ac' ,'Tex_like' = '#9566ac' , 'None' = 'lightgrey'       
                            )

PB_cluster_colors <-  c('CD8+TEM-B'='#cb50a0',
                    'CD8+TEM-K'='#db72b4',#
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
                    'CD8+TCM'='#E59CC4')

cluster_colors <- c(cluster_colors, PB_cluster_colors)

lbls_func <- 
function(x){
    x <- gsub('CD4.*Naive', 'mature CD4 T', x)
    x <- gsub('CD8.*Naive', 'mature CD8 T', x)
    x <- gsub('CD4.*pre_T', 'pre-CD4 T', x)
    x <- gsub('CD8.*pre_T', 'pre-CD8 T', x)
    x <- gsub('CD4.*Memory_like', 'memory-like CD4 T', x)
    x <- gsub('CD8.*Memory_like', 'memory-like CD8 T', x)
    #x <- gsub('Agonist', 'agonist', x)
    x <- gsub('thymus_', 'thymus ', x)
    x <- gsub('PB_', 'PB  ', x)
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x <- gsub('thymus', 'Thy', x)
    x
}

tree_color <- c("#9E6BAB", "#CFE3C6", "#E3CEE4", "#F3746C", "#86382A", "#ABA300", "#D6D4EB"
, "#B9A96B", "#FF68A1", "#EAA944", "#7AAF93", "#7A8DBB", "#7673AE", "#396F68"
, "#ECAFCF", "#EACF68", "#F7DDD4", "#EBDBE4", "#66CEF6", "#F8F4A8", "#C35338"
, "#EF5276", "#A0D7C9", "#63B472", "#F9DBE5", "#0CB702", "#F48930", "#6B6A6B"
, "#27BDCF", "#F8BFAF", "#F5C785", "#DEEAB1", "#217CB8", "#31FEB3", "#74517B"
, "#588198", "#CAA57D", "#9C99C4", "#2D563D", "#FF77AB", "#9F8CFF", "#D5E7F7"
, "#22A3FF", "#00E8F7", "#BB4A94", "#69B4CE", "#C9BDB2")

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


thymus_identity <- 
as.character(Idents(T_dev_3))
names(thymus_identity) <- names(Idents(T_dev_3))
thymus_T_idents <- c('αβ_Entry','Agonist_1','CD4+pre_T', 'CD8+pre_T', 'Agonist_2', 
                    'CD4+Naive', 'CD8+Naive', 'IFN_response','pre_Treg1', 'pre_Treg2','Treg','CD4+Memory_like',
                    'CD8+Memory_like')
thymus_T_idents2 <- c('DP_P', 'DP_Q', 'HSP_DP', 
                      'αβ_Entry','Agonist_1','CD4+pre_T', 'CD8+pre_T', 'Agonist_2', 
                    'CD4+Naive', 'CD8+Naive', 'IFN_response','pre_Treg1','pre_Treg2', 'Treg', 'CD4+Memory_like',
                    'CD8+Memory_like')

PB_identity <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell_tcr_idents.rds')
##############只包含abT
PB_T_idents <- c('CD4+RTE','CD4+TN','CD4+TCM','CD4+TEFF','CD4+TEM','CD4+CTL','Treg',
                 'CD8+RTE','CD8+TN','CD8+TCM',
'CD8+TEM-K','CD8+TEM-B', 'Tex','IFN_T', 'T_un')

setwd('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/data/')
files <- dir()
files <- files[-grep('metadata.txt', files)]
thymus_tcr <- list()
thymus_tcr$meta <- read.table('metadata.txt', header = T)  
thymus_tcr$data <- 
lapply(files, function(i){
    df_iter <- read.csv(i)
    df_iter$Sample <- strsplit(i, split = '\\.')[[1]][1]
    df_iter$age <- thymus_tcr$meta[match(df_iter$Sample, thymus_tcr$meta$Sample), 'age']
    df_iter$age_3 <- thymus_tcr$meta[match(df_iter$Sample, thymus_tcr$meta$Sample), 'age_3']
    df_iter$identity <- thymus_identity[match(df_iter$barcode, names(thymus_identity))]
    df_iter$tissue <- 'thymus'
    df_iter$tissue_identity <- paste0( df_iter$tissue, '_', df_iter$identity)
    df_iter <- df_iter[df_iter$identity %in% thymus_T_idents, ]

    df_iter <- df_iter[df_iter$productive == 'True', ]
    tbl=table(df_iter$barcode, df_iter$chain)
    logi <- apply(tbl, 1, function(x){
        x['TRA'] > 0 & x['TRB'] > 0 
    })  
    df_iter <- df_iter[df_iter$barcode %in% rownames(tbl)[logi], ]
    return(df_iter)
})
names(thymus_tcr$data) <- sapply(strsplit(files, split = '\\.'), function(x) x[1])                          
setwd(plot_dir)                           

thymus_tcr$meta$age_3 <- factor(thymus_tcr$meta$age_3, levels =  c('<=12', '13-39', '>=40'))

setwd('/data1/02.private/dengyj/analysis/thymus/TCR/PB/data/')
columns_selected <- c('barcode','is_cell','contig_id','high_confidence','length','chain','v_gene',
'd_gene','j_gene','c_gene','full_length','productive','cdr3','cdr3_nt',
'reads','umis','raw_clonotype_id','raw_consensus_id')
files <- dir()
files <- files[-grep('metadata.txt', files)]
PB_tcr <- list()
PB_tcr$meta <- read.table('metadata.txt', header = T)  
PB_tcr$data <- 
lapply(files, function(i){
    df_iter <- read.csv(i)
    df_iter <- df_iter[, columns_selected]
    df_iter$Sample <- strsplit(i, split = '\\.')[[1]][1]
    df_iter$age <- PB_tcr$meta[match(df_iter$Sample, PB_tcr$meta$Sample), 'age']
    df_iter$age_3 <- PB_tcr$meta[match(df_iter$Sample, PB_tcr$meta$Sample), 'age_3']
    df_iter$identity <- PB_identity[match(df_iter$barcode, names(PB_identity))]
    df_iter$tissue <- 'PB'
    df_iter$tissue_identity <- paste0( df_iter$tissue, '_', df_iter$identity)
    df_iter <- df_iter[df_iter$identity %in% PB_T_idents, ]

    df_iter <- df_iter[grepl('True|true', df_iter$productive), ]####true是cellranger 6的结果
    tbl=table(df_iter$barcode, df_iter$chain)
    logi <- apply(tbl, 1, function(x){
        x['TRA'] > 0 & x['TRB'] > 0 
    })  
    df_iter <- df_iter[df_iter$barcode %in% rownames(tbl)[logi], ]    
    
    return(df_iter)
})
names(PB_tcr$data) <- sapply(strsplit(files, split = '\\.'), function(x) x[1])                          
setwd(plot_dir)   
                              
PB_tcr$meta$age_3 <- factor(PB_tcr$meta$age_3, levels =  c('<=12', '13-39', '40-99', '>=100'))

###Fig.7a
options(repr.plot.width=7, repr.plot.height=7)

for(i in 1:length(thymus_tcr$data)){
    df_iter <- thymus_tcr$data[[i]]
    df_iter <- df_iter[!duplicated(df_iter$barcode), ]  
    df_iter <- data.frame(table(df_iter$raw_clonotype_id))
    set.seed(123)
    df_iter$color <- sample(tree_color, nrow(df_iter), replace=TRUE)
    name <- names(thymus_tcr$data)[i]
    age <- thymus_tcr$meta[thymus_tcr$meta$Sample %in% name, 'age']
    p1 <- ggplot(data=df_iter,aes(area=Freq,fill = color))+
    geom_treemap(colour = 'black')+guides(fill = F)+
    scale_fill_manual(values = tree_color) #+ 
    print(p1+labs(title = paste0(gsub('-10X5', '', name), '_', age, 'Y')))
    ggsave(paste0(gsub('-10X5', '', name), '_', age, 'Y', '.png'), width = 7, height = 7,plot = p1)
}

###Fig.7b
div_by_cells_age <- lapply(thymus_tcr$data, 
                           function(x) get_div_by_cells(df = x, group_var = 'age_3', method = 'pielou')
                          )
div_by_cells_age <- do.call(rbind, div_by_cells_age)
div_by_cells_age$age_3 <- factor(div_by_cells_age$age_3, levels =  c('<=12', '13-39', '>=40'))
stat_test <- wilcox_test(
 Diversity ~ age_3, data = div_by_cells_age,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.3,fun = 'mean')
stat_test$y.position <- seq(min(stat_test$y.position), max(stat_test$y.position), 
                            length.out = length(stat_test$y.position))
p <- ggplot(div_by_cells_age, aes(age_3, Diversity))+
geom_boxplot(mapping = aes(fill = age_3))+
scale_x_discrete(#limits = c('<=12', '13-39', '40-99', '>=100'), 
                 labels= age_group_lbls_func)+
geom_point( mapping = aes(x=age_3, y=Diversity), 
           position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors,labels= c('Prepuberal', 'Adult', 'Aged'))+
labs(x = NULL, y = "Pielou's evenness index", fill = 'Age')+
guides(fill = F)+
scale_y_continuous(expand = c(0, 0.005))+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt)
options(repr.plot.width=3, repr.plot.height=4.5)
p
saveRDS(div_by_cells_age, file = 'thymus_TCR_div.rds')
ggsave('thymus_TCR_div.svg', width = 3, height =4.5)
#my_saveplot_fun(p, 'thymus_TCR_div', width = 3, height =4.5)

###Fig.7c
options(repr.plot.width=7, repr.plot.height=7)
for(i in 1:length(PB_tcr$data)){
    df_iter <- PB_tcr$data[[i]]
    df_iter <- df_iter[!duplicated(df_iter$barcode), ]  
    df_iter <- data.frame(table(df_iter$raw_clonotype_id))
    set.seed(1234)
    df_iter$color <- sample(tree_color, nrow(df_iter), replace=TRUE)
    name <- names(PB_tcr$data)[i]
    age <- PB_tcr$meta[PB_tcr$meta$Sample %in% name, 'age']
    p1 <- ggplot(data=df_iter,aes(area=Freq,fill = color))+
        geom_treemap(colour = 'black')+guides(fill = F)+
        scale_fill_manual(values = tree_color) #+ 
        #ggtitle(paste0(gsub('-10X5', '', name), '_', age, 'Y'))
    print(p1+labs(title = paste0(gsub('-10X5', '', name), '_', age, 'Y')))
    ggsave(paste0(gsub('-10X5', '', name), '_', age, 'Y', '.png'), width = 7, height = 7,plot = p1)
}

###Fig.7d
div_by_cells_age <- lapply(PB_tcr$data, 
                           function(x) get_div_by_cells(df = x, group_var = 'age_3', method = 'pielou')
                          )
div_by_cells_age <- do.call(rbind, div_by_cells_age)
div_by_cells_age$age_3 <- factor(div_by_cells_age$age_3, levels =  c('<=12', '13-39', '40-99', '>=100'))

stat_test <- wilcox_test(
 Diversity ~ age_3, data = div_by_cells_age,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.25,fun = 'mean')
stat_test$y.position <- seq(min(stat_test$y.position), max(stat_test$y.position), 
                            length.out = length(stat_test$y.position))


p <- ggplot(div_by_cells_age, aes(age_3, Diversity))+
geom_boxplot(mapping = aes(fill = age_3))+
scale_x_discrete(#limits = c('<=12', '13-39', '40-99', '>=100'), 
                 labels=age_group_lbls_func)+
geom_point( mapping = aes(x=age_3, y=Diversity), 
           position = position_dodge(0.9), size = 0.25)+
scale_fill_manual(values = age_colors,labels= c('Prepuberal', 'Adult', 'Aged', 'Centenarian'))+

labs(x = NULL, y = "Pielou's evenness index", fill = 'Age')+
guides(fill = F)+
scale_y_continuous(expand = c(0, 0.1))+
theme(#plot.margin = margin(0, 0, 0, 0, 'lines'),
    panel.background = element_rect(fill = 'white', color = NA), 
    panel.border = element_rect(fill = NA, color = 'black'),
      axis.title = element_text(size=13.75*15/13),axis.text.y = element_text(size = 13.75*15/13),
      axis.text.x = element_text(size = 13.75*15/13, angle = 45, hjust=1),
      plot.title = element_text(size=13.75*15/13,hjust=0.5),
      #axis.title.y = element_blank(),
     legend.title = element_text(size = 13.75*15/13),legend.text = element_text(size = 13.75*15/13))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75*15/13/.pt)
options(repr.plot.width=3, repr.plot.height=4.5)
p
saveRDS(div_by_cells_age, file = 'PB_TCR_div.rds')
ggsave('PB_TCR_div.svg', width = 3, height = 4.5)
#my_saveplot_fun(p, 'PB_TCR_div', width = 3, height = 4.5)





###Extended Data Fig.8a-c

tcr_clu_age_test <- readRDS('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/tcr_clu_age_test.rds')

classic_T_idents <- c('DN','DP_P', 'DP_Q',  
              'αβ_Entry','CD4+pre_T', 'CD8+pre_T', 
                    'CD4+Naive', 'CD8+Naive')

trav_plot_df <- lapply(tcr_clu_age_test[['trav_clu_age_test']], function(x){x$gene <- rownames(x);x})
trav_plot_df <- do.call(rbind, trav_plot_df)
#p_val_df <- p_val_df[, -grep('p_val$', colnames(p_val_df))]
trav_plot_df <- melt(trav_plot_df, id.vars = c('cluster', 'gene'))
trav_plot_df$age_group <- gsub('_.*$', '', trav_plot_df$variable)
trav_plot_df$variable <- gsub('^.*\\d_', '', trav_plot_df$variable)
trav_plot_df <- dcast(trav_plot_df, cluster+gene+age_group~variable, value.var = 'value')

trav_plot_df <- ddply(trav_plot_df, c('gene'), transform, useage = scale(useage))

age_lbls <- c('<=12' = 'Prepubertal', '13-39' = 'Adult', '>=40' = 'Aged')
trav_plot_df$cluster <- lbls_func(as.character(trav_plot_df$cluster))
trav_plot_df$cluster <- factor(trav_plot_df$cluster ,levels = 
                                  lbls_func(setdiff(classic_T_idents, 'DN')))
trav_plot_df$age_group <- factor(age_lbls[trav_plot_df$age_group], levels = 
                                c('Prepubertal', 'Adult', 'Aged'))

trav_plot_df <- trav_plot_df[!trav_plot_df$gene %in% trav_plot_df$gene[is.na(trav_plot_df$p_val_adj)],]



val=unique(trav_plot_df$gene)
gene_order=data.frame(gene = val, loc = gsub('TRAV', '',val))
gene_order$loc1 <- as.numeric(gsub('-.*$|/DV.*$', '', gene_order$loc))
gene_order$loc2 <- as.numeric(gsub('^.*-|^.*/DV', '', gene_order$loc))
gene_order <- gene_order[order(gene_order$loc1, gene_order$loc2),]


options(repr.plot.width=12, repr.plot.height=2.5)
ggplot(trav_plot_df[trav_plot_df$gene %in% unique(trav_plot_df$gene[trav_plot_df$p_val_adj < 0.05]), ],
       aes(age_group , gene, fill = useage,, size = -log10(p_val_adj)
                        ))+
#geom_tile()+
geom_point(shape = 21)+
scale_fill_gradientn(colours = viridis::viridis(5))+
# scale_y_discrete(limits = gene_order$gene)+
labs(x = NULL, y = NULL, size = '-log10 P.adj')+
# scale_x_discrete(limits = )+
facet_wrap(.~cluster, nrow = 1)+
guides(size = guide_legend(ncol = 2, override.aes = list(fill = 'black')), 
       fill = guide_colorbar(order = 1,direction = 'horizontal', title.position = 'top',
                             title = 'scaled usage'))+
scale_size_continuous(range = c(0,5))+
theme(axis.text.x = element_text(size= 13.75, angle = 45, vjust=1,hjust = 1), strip.background = element_blank(),
      legend.title = element_text(size=13.75),
      strip.text.x = element_text(size= 13.75),
     legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), # axis.ticks=element_blank(),
      axis.text.y=element_text(size= 13.75),
      axis.line=element_blank(), axis.title=element_text(size = 17.5),
      legend.text = element_text(size=13.75))+
geom_tile(data = trav_plot_df[trav_plot_df$p_val_adj < 0.05, ], color = 'red',fill = NA, size = 0.5)+
scale_x_discrete(limits = c('Prepubertal', 'Adult', 'Aged'), label = age_group_lbls_func)
ggsave('TRAV.svg', width = 12, height = 2.5)

traj_plot_df <- lapply(tcr_clu_age_test[['traj_clu_age_test']], function(x){x$gene <- rownames(x);x})
traj_plot_df <- do.call(rbind, traj_plot_df)
#p_val_df <- p_val_df[, -grep('p_val$', colnames(p_val_df))]
traj_plot_df <- melt(traj_plot_df, id.vars = c('cluster', 'gene'))
traj_plot_df$age_group <- gsub('_.*$', '', traj_plot_df$variable)
traj_plot_df$variable <- gsub('^.*\\d_', '', traj_plot_df$variable)
traj_plot_df <- dcast(traj_plot_df, cluster+gene+age_group~variable, value.var = 'value')

traj_plot_df <- ddply(traj_plot_df, c('gene'), transform, useage = scale(useage))

age_lbls <- c('<=12' = 'Prepubertal', '13-39' = 'Adult', '>=40' = 'Aged')
traj_plot_df$cluster <- lbls_func(as.character(traj_plot_df$cluster))
traj_plot_df$cluster <- factor(traj_plot_df$cluster ,levels = 
                                  lbls_func(setdiff(classic_T_idents, 'DN')))
traj_plot_df$age_group <- factor(age_lbls[traj_plot_df$age_group], levels = 
                                c('Prepubertal', 'Adult', 'Aged'))

traj_plot_df <- traj_plot_df[!traj_plot_df$gene %in% traj_plot_df$gene[is.na(traj_plot_df$p_val_adj)],]



val=unique(traj_plot_df$gene)
gene_order=data.frame(gene = val, loc = gsub('traj', '',val))
gene_order$loc1 <- as.numeric(gsub('-.*$|/DV.*$', '', gene_order$loc))
gene_order$loc2 <- as.numeric(gsub('^.*-|^.*/DV', '', gene_order$loc))
gene_order <- gene_order[order(gene_order$loc1, gene_order$loc2),]


options(repr.plot.width=12, repr.plot.height=2.5)
ggplot(traj_plot_df[traj_plot_df$gene %in% unique(traj_plot_df$gene[traj_plot_df$p_val_adj < 0.05]), ],
       aes(age_group , gene, fill = useage,, size = -log10(p_val_adj)
                        ))+
#geom_tile()+
geom_point(shape = 21)+
scale_fill_gradientn(colours = viridis::viridis(5))+
# scale_y_discrete(limits = gene_order$gene)+
labs(x = NULL, y = NULL, size = '-log10 P.adj')+
# scale_x_discrete(limits = )+
facet_wrap(.~cluster, nrow = 1)+
guides(size = guide_legend(ncol = 2, override.aes = list(fill = 'black')), 
       fill = guide_colorbar(order = 1,direction = 'horizontal', title.position = 'top',
                             title = 'scaled usage'))+
scale_size_continuous(range = c(0,5))+
theme(axis.text.x = element_text(size= 13.75, angle = 45, vjust=1,hjust = 1), strip.background = element_blank(),
      legend.title = element_text(size=13.75),
      strip.text.x = element_text(size= 13.75),
     legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), # axis.ticks=element_blank(),
      axis.text.y=element_text(size= 13.75),
      axis.line=element_blank(), axis.title=element_text(size = 17.5),
      legend.text = element_text(size=13.75))+
geom_tile(data = traj_plot_df[traj_plot_df$p_val_adj < 0.05, ], color = 'red',fill = NA, size = 0.5)+
scale_x_discrete(limits = c('Prepubertal', 'Adult', 'Aged'), label = age_group_lbls_func)
ggsave('TRAJ.svg', width = 12, height = 2.5)

trbv_plot_df <- lapply(tcr_clu_age_test[['trbv_clu_age_test']], function(x){x$gene <- rownames(x);x})
trbv_plot_df <- do.call(rbind, trbv_plot_df)
#p_val_df <- p_val_df[, -grep('p_val$', colnames(p_val_df))]
trbv_plot_df <- melt(trbv_plot_df, id.vars = c('cluster', 'gene'))
trbv_plot_df$age_group <- gsub('_.*$', '', trbv_plot_df$variable)
trbv_plot_df$variable <- gsub('^.*\\d_', '', trbv_plot_df$variable)
trbv_plot_df <- dcast(trbv_plot_df, cluster+gene+age_group~variable, value.var = 'value')

trbv_plot_df <- ddply(trbv_plot_df, c('gene'), transform, useage = scale(useage))

age_lbls <- c('<=12' = 'Prepubertal', '13-39' = 'Adult', '>=40' = 'Aged')
trbv_plot_df$cluster <- lbls_func(as.character(trbv_plot_df$cluster))
trbv_plot_df$cluster <- factor(trbv_plot_df$cluster ,levels = 
                                  lbls_func(classic_T_idents))
trbv_plot_df$age_group <- factor(age_lbls[trbv_plot_df$age_group], levels = 
                                c('Prepubertal', 'Adult', 'Aged'))

trbv_plot_df <- trbv_plot_df[!trbv_plot_df$gene %in% trbv_plot_df$gene[is.na(trbv_plot_df$p_val_adj)],]



val=unique(trbv_plot_df$gene)
gene_order=data.frame(gene = val, loc = gsub('trbv', '',val))
gene_order$loc1 <- as.numeric(gsub('-.*$|/DV.*$', '', gene_order$loc))
gene_order$loc2 <- as.numeric(gsub('^.*-|^.*/DV', '', gene_order$loc))
gene_order <- gene_order[order(gene_order$loc1, gene_order$loc2),]


options(repr.plot.width=13, repr.plot.height=2.5)
ggplot(trbv_plot_df[trbv_plot_df$gene %in% unique(trbv_plot_df$gene[trbv_plot_df$p_val_adj < 0.05]), ],
       aes(age_group , gene, fill = useage,, size = -log10(p_val_adj)
                        ))+
#geom_tile()+
geom_point(shape = 21)+
scale_fill_gradientn(colours = viridis::viridis(5))+
# scale_y_discrete(limits = gene_order$gene)+
labs(x = NULL, y = NULL, size = '-log10 P.adj')+
# scale_x_discrete(limits = )+
facet_wrap(.~cluster, nrow = 1)+
guides(size = guide_legend(ncol = 2, override.aes = list(fill = 'black')), 
       fill = guide_colorbar(order = 1,direction = 'horizontal', title.position = 'top',
                             title = 'scaled usage'))+
scale_size_continuous(range = c(0,5))+
theme(axis.text.x = element_text(size= 13.75, angle = 90, vjust=0.5,hjust = 1), strip.background = element_blank(),
      legend.title = element_text(size=13.75),
      strip.text.x = element_text(size= 13.75),
     legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), # axis.ticks=element_blank(),
      axis.text.y=element_text(size= 13.75),
      axis.line=element_blank(), axis.title=element_text(size = 17.5),
      legend.text = element_text(size=13.75))+
geom_tile(data = trbv_plot_df[trbv_plot_df$p_val_adj < 0.05, ], color = 'red',fill = NA, size = 0.5)+
scale_x_discrete(limits = c('Prepubertal', 'Adult', 'Aged'), label = age_group_lbls_func)

trbj_plot_df <- lapply(tcr_clu_age_test[['trbj_clu_age_test']], function(x){x$gene <- rownames(x);x})
trbj_plot_df <- do.call(rbind, trbj_plot_df)
#p_val_df <- p_val_df[, -grep('p_val$', colnames(p_val_df))]
trbj_plot_df <- melt(trbj_plot_df, id.vars = c('cluster', 'gene'))
trbj_plot_df$age_group <- gsub('_.*$', '', trbj_plot_df$variable)
trbj_plot_df$variable <- gsub('^.*\\d_', '', trbj_plot_df$variable)
trbj_plot_df <- dcast(trbj_plot_df, cluster+gene+age_group~variable, value.var = 'value')

trbj_plot_df <- ddply(trbj_plot_df, c('gene'), transform, useage = scale(useage))

age_lbls <- c('<=12' = 'Prepubertal', '13-39' = 'Adult', '>=40' = 'Aged')
trbj_plot_df$cluster <- lbls_func(as.character(trbj_plot_df$cluster))
trbj_plot_df$cluster <- factor(trbj_plot_df$cluster ,levels = 
                                  lbls_func(classic_T_idents))
trbj_plot_df$age_group <- factor(age_lbls[trbj_plot_df$age_group], levels = 
                                c('Prepubertal', 'Adult', 'Aged'))

trbj_plot_df <- trbj_plot_df[!trbj_plot_df$gene %in% trbj_plot_df$gene[is.na(trbj_plot_df$p_val_adj)],]



val=unique(trbj_plot_df$gene)
gene_order=data.frame(gene = val, loc = gsub('trbj', '',val))
gene_order$loc1 <- as.numeric(gsub('-.*$|/DV.*$', '', gene_order$loc))
gene_order$loc2 <- as.numeric(gsub('^.*-|^.*/DV', '', gene_order$loc))
gene_order <- gene_order[order(gene_order$loc1, gene_order$loc2),]


options(repr.plot.width=12, repr.plot.height=2.5)
ggplot(trbj_plot_df[trbj_plot_df$gene %in% unique(trbj_plot_df$gene[trbj_plot_df$p_val_adj < 0.05]), ],
       aes(age_group , gene, fill = useage,, size = -log10(p_val_adj)
                        ))+
#geom_tile()+
geom_point(shape = 21)+
scale_fill_gradientn(colours = viridis::viridis(5))+
# scale_y_discrete(limits = gene_order$gene)+
labs(x = NULL, y = NULL, size = '-log10 P.adj')+
# scale_x_discrete(limits = )+
facet_wrap(.~cluster, nrow = 1)+
guides(size = guide_legend(ncol = 2, override.aes = list(fill = 'black')), 
       fill = guide_colorbar(order = 1,direction = 'horizontal', title.position = 'top',
                             title = 'scaled usage'))+
scale_size_continuous(range = c(0,5))+
theme(axis.text.x = element_text(size= 13.75, angle = 45, vjust=1,hjust = 1), strip.background = element_blank(),
      legend.title = element_text(size=13.75),
      strip.text.x = element_text(size= 13.75),
     legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), # axis.ticks=element_blank(),
      axis.text.y=element_text(size= 13.75),
      axis.line=element_blank(), axis.title=element_text(size = 17.5),
      legend.text = element_text(size=13.75))+
geom_tile(data = trbj_plot_df[trbj_plot_df$p_val_adj < 0.05, ], color = 'red',fill = NA, size = 0.5)+
scale_x_discrete(limits = c('Prepubertal', 'Adult', 'Aged'), label = age_group_lbls_func)
ggsave('TRBJ.svg', width = 12, height = 2.5)

###Extended Data Fig.8d top
div_by_cells_idents <- lapply(thymus_tcr$data, 
                           function(x) get_div_by_cells(df = x, group_var = 'identity', method = 'pielou')
                          )
div_by_cells_idents <- do.call(rbind, div_by_cells_idents)
clusters_used <- thymus_T_idents
div_by_cells_idents$identity <- factor(div_by_cells_idents$identity, levels = rev(clusters_used))



p1 <- ggplot(div_by_cells_idents[div_by_cells_idents$identity %in% clusters_used, ], 
       aes(Diversity, identity, fill = identity))+ geom_boxplot()+
geom_point(position = position_dodge(0.9), size = 0.25)+
stat_compare_means( method="kruskal.test",  size = 13.75/.pt, vjust = -0.3, label.x.npc = 0.3)+
theme(axis.text.y = element_text(angle = 0, hjust=1, size = 13.75), 
     axis.text.x = element_text(size = 13.75), 
     axis.title.x = element_text(size = 13.75),
     legend.title = element_text(size = 13.75),legend.text = element_text(size = 13.75), 
      panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA),
      plot.title = element_text(size = 13.75, hjust = 0.5)
     )+
guides(fill = guide_legend(override.aes = list(size=0.1), keywidth = 1, keyheight = 1))+
scale_fill_manual(values = cluster_colors,limits=clusters_used,
                 label = lbls_func)+
labs(x = "Pielou's evenness index", fill = 'Identity', y = NULL, title = 'thymus')+
scale_x_continuous(expand = c(0, 0.02))+
scale_y_discrete(label = lbls_func)+
guides(fill = F)
options(repr.plot.width=5, repr.plot.height=3.75)
p1
saveRDS(p1$data, file = 'thymus_idents_TCR_div.rds')
#ggsave('thymus_idents_TCR_div.svg', width = 8, height = 6)
#my_saveplot_fun(p, 'thymus_idents_TCR_div', width = 6, height = 4.5)


###Extended Data Fig.8d bottom
div_by_cells_idents <- lapply(PB_tcr$data, 
                           function(x) get_div_by_cells(df = x, group_var = 'identity', method = 'pielou')
                          )
div_by_cells_idents <- do.call(rbind, div_by_cells_idents)


clusters_used <- PB_T_idents
div_by_cells_idents$identity <- factor(div_by_cells_idents$identity, levels = rev(clusters_used))


####Pielous evenness was NA when there is only one clone in a cell cluster from a certain donor
p2 <- ggplot(div_by_cells_idents[div_by_cells_idents$identity %in% clusters_used, ], 
       aes(Diversity, identity, fill = identity))+ geom_boxplot()+
geom_point(position = position_dodge(0.9), size = 0.25)+
stat_compare_means( method="kruskal.test", label.x.npc = 0.3, size = 13.75/.pt, vjust = -0.4)+
theme(axis.text.y = element_text(angle = 0, hjust=1, size = 13.75), 
     axis.text.x = element_text(size = 13.75), 
     axis.title.x = element_text(size = 13.75),
     legend.title = element_text(size = 13.75),legend.text = element_text(size = 13.75), 
      panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA),
      plot.title = element_text(size = 13.75, hjust = 0.5)
     )+
guides(fill = guide_legend(override.aes = list(size=0.1), keywidth = 1, keyheight = 1))+
scale_fill_manual(values = cluster_colors,limits=clusters_used, label = lbls_func)+
scale_y_discrete(label = lbls_func)+
labs(x = "Pielou's evenness index", fill = 'Identity', y = NULL, title = 'PB')+
scale_x_continuous(expand = c(0, 0.075))+
guides(fill = F)
options(repr.plot.width=5, repr.plot.height=3.5)
p2
saveRDS(p2$data[!is.na(p2$data$Diversity), ], file = 'PB_idents_TCR_div.rds')
#ggsave('PB_idents_TCR_div.svg', width = 8, height = 6)
#my_saveplot_fun(p, 'PB_idents_TCR_div', width = 6, height = 4.5)

options(warn = -1)
plot_list <- get_gg_align(list(p1, p2))
options(repr.plot.width=5, repr.plot.height=8)
wrap_plots(list(p1, p2))+plot_layout(ncol = 1)
ggsave('thymus_PB_idents_TCR_div.svg', width = 5, height = 8)
options(warn = 0)



combined_cluster_colors <- c('PB_CD8+TEM-B'='#cb50a0','PB_CD8+TEM-K'='#db72b4','PB_CD8+RTE'='#276d6d',
'PB_CD4+TN'='#7bacd9','PB_CD4+CTL'='#7c75ce','PB_CD56_Dim_NK'='#c8af1e',
'PB_CD56_Bright_NK'='#d99333','PB_CD27lo_Vδ1T'='#93c88e','PB_CD8+TN'='#DDC4DA',
'PB_IFN_response'='#00E8F7','PB_CD4+TEM'='#868ae3','PB_Treg'='#6e7cad','PB_CD4+TCM'='#80a2e6',
'PB_MAIT'='#b1b772','PB_CD4+RTE'='#66c5b9','PB_CD27lo_Vδ2Vγ9T'='#A0D7C9','PB_proliferating'='#bb9400',
'PB_NKT'='#e3cd82','PB_CD27+Vδ1T'='#9cb39b','PB_T_un'='#747474','PB_Tex'='#A1B4BB','PB_CD4+TEFF'='#55b3db',
  'PB_CD8+TCM'='#E59CC4', 
'thymus_pre_Treg1'='#E5D2DD','thymus_CD4+Memory_like'='#808040','thymus_CD8+Memory_like'='#F3B1A0',
  'thymus_CD4+GZMK+T'='#57C3F3',
'thymus_IFN_response'='#476D87','thymus_Treg'='#E95C59','thymus_CD4+pre_T'='#83c3bf',
'thymus_DP_P'='#AB3282','thymus_DP_Q'='#A8B0E3','thymus_Agonist_2'='#a83117','thymus_Exhausted_T'='#8C549C',
'thymus_αβ_Entry'='#c28578','thymus_HSP_DP'='#9FA3A8','thymus_IFN_response_2'='#CCE0F5',
  'thymus_CD8+pre_T'='#625D9E',
'thymus_Agonist_1'='#a6c050','thymus_HSP_SP'='#d98ed6','thymus_pre_Treg2'='#9566ac',
  'thymus_CD8+Naive'='#53A85F','thymus_CD4+Naive'='#D6E7A3')

combined_lvls <- c('thymus_αβ_Entry','thymus_Agonist_1','thymus_CD4+pre_T','thymus_CD8+pre_T','thymus_Agonist_2',
'thymus_CD4+Naive','thymus_CD8+Naive','thymus_IFN_response','thymus_pre_Treg1','thymus_pre_Treg2',
                   'thymus_Treg',
'thymus_CD4+Memory_like','thymus_CD8+Memory_like',
         'PB_CD4+RTE','PB_CD4+TN','PB_CD4+TCM','PB_CD4+TEFF','PB_CD4+TEM','PB_CD4+CTL','PB_Treg',

          'PB_CD8+RTE','PB_CD8+TN','PB_CD8+TCM','PB_CD8+TEM-K','PB_CD8+TEM-B','PB_Tex','PB_T_un')

thymus_tcr$meta$age_3 <- as.character(thymus_tcr$meta$age_3)
thymus_tcr$meta$age_3[thymus_tcr$meta$age_3 == '>=40'] <- '40-99'
for(i in 1:length(thymus_tcr$data)){
    df_iter <- thymus_tcr$data[[i]]
    df_iter$age_3[df_iter$age_3 == '>=40'] <- '40-99'
    thymus_tcr$data[[i]] <- df_iter
}
thymus_tcr$meta$age_3 <- factor(thymus_tcr$meta$age_3, levels =  c('<=12', '13-39', '40-99'))

combined_tcr <- 
list(
    data = c(thymus_tcr$data, PB_tcr$data),
    meta = rbind(thymus_tcr$meta, PB_tcr$meta)
)

sample_used <- c('T20041304F-10X5','T20041305F-10X5','T20041307F-10X5', 
                    'T20041501F-10X5','T20041502M-10X5','T20041503M-10X5','T20041504F-10X5',
                    'T20041601M-10X5','T20041603M-10X5','T20041604F-10X5',
                   'P20041304F-10X5','P20041305F-10X5','P20041307F-10X5',
                   'P20041501F-10X5','P20041502M-10X5','P20041503M-10X5','P20041504F-10X5',
                   'P20041601M-10X5','P20041603M-10X5','P20041604F-10X5')
combined_tcr$data <- combined_tcr$data[names(combined_tcr$data) %in% sample_used]
combined_tcr$meta <- combined_tcr$meta[combined_tcr$meta$Sample %in% sample_used, ]

age_order <- c('4','9', '22', '24', '35', '35', '46', '53', '56', '69')

sample_list <- c('20041504F-10X5', '20041502M-10X5', '20041503M-10X5',
        '20041307F-10X5', '20041305F-10X5', '20041501F-10X5', '20041604F-10X5', '20041603M-10X5', 
         '20041304F-10X5', '20041601M-10X5')

combined_clonotype_list <- 
lapply(sample_list, function(i){
    t_df = combined_tcr$data[[paste0('T', i)]];#t_df$tissue <- 'thymus'
    p_df = combined_tcr$data[[paste0('P', i)]];#p_df$tissue <- 'PB'
    df <- rbind(t_df, p_df)
    df <- df[df$cdr3_nt != 'None', ]
    df <- df[order(df$barcode, df$cdr3_nt), ]
    df_cp <- df[!duplicated(df$barcode), ]
    nt_df <- ddply(df, 'barcode', summarise, 
                   CDR3.nt = paste(cdr3_nt, collapse = ';'))
    nt_df$v_gene <- sapply(1:nrow(nt_df), function(j){
        paste(df[df$barcode %in% nt_df$barcode[j], 'v_gene'], collapse = ';')})
    nt_df$d_gene <- sapply(1:nrow(nt_df), function(j){
        paste(df[df$barcode %in% nt_df$barcode[j], 'd_gene'], collapse = ';')})
    nt_df$j_gene <- sapply(1:nrow(nt_df), function(j){
        paste(df[df$barcode %in% nt_df$barcode[j], 'j_gene'], collapse = ';')})
    nt_df$c_gene <- sapply(1:nrow(nt_df), function(j){
        paste(df[df$barcode %in% nt_df$barcode[j], 'c_gene'], collapse = ';')})

    nt_df$tissue_identity <- df_cp[match(nt_df$barcode, df_cp$barcode), 'tissue_identity']
    nt_df$tissue <- df_cp[match(nt_df$barcode, df_cp$barcode), 'tissue']
    nt_df$age <- df_cp[match(nt_df$barcode, df_cp$barcode), 'age']
    nt_df$age_3 <- df_cp[match(nt_df$barcode, df_cp$barcode), 'age_3']
    nt_df$Sample <- df_cp[match(nt_df$barcode, df_cp$barcode), 'Sample']
    nt_df$Gene <- paste(nt_df$v_gene, nt_df$d_gene, nt_df$j_gene, nt_df$c_gene, sep = ';')
    nt_df$clonotype <- paste(nt_df$Gene, nt_df$CDR3.nt, sep = ';')
    return(nt_df)
})

names(combined_clonotype_list) <- sample_list



options(repr.plot.width=4, repr.plot.height=4.1)

clu_list <- lapply(1:length(sample_list), function(i){
    nt_df <-  combined_clonotype_list[[i]]
    t_nt_df <- nt_df[nt_df$tissue == 'thymus', ]
    p_nt_df <- nt_df[nt_df$tissue == 'PB', ]
    p_nt_df_2 <- p_nt_df[p_nt_df$clonotype %in% t_nt_df$clonotype,]
    t_nt_df_2 <- t_nt_df[t_nt_df$clonotype %in% p_nt_df$clonotype,]
    t_max_val <- max(table(t_nt_df_2$clonotype))
    p_max_val <- max(table(p_nt_df_2$clonotype))
    max_val <- max(c(t_max_val,p_max_val))
    stats <- table(p_nt_df_2$clonotype)
    stats <- stats[order(stats, decreasing = T)]
    p1 <- ggplot(p_nt_df_2, aes(clonotype, fill = tissue_identity))+geom_bar(color = 'black', size = 0.3)+
        theme(axis.text.y=element_blank(),
             axis.line = element_blank(), axis.ticks = element_blank(),
             panel.background = element_rect(fill = 'white', colour = NA),
             panel.border = element_rect(fill = NA, colour = 'black'), 
             axis.text.x = element_text(size=13.75),
             axis.title = element_text(size=13.75),
             legend.text =  element_text(size=13.75),
             legend.title =  element_text(size=13.75),
             plot.title =  element_text(size=13.75, hjust=0.5))+
        labs(x = 'Clonotype', y='Counts', fill = 'PB_identity', title = 'PB')+
        scale_x_discrete(limits = rev(names(stats)))+
        scale_fill_manual(values = combined_cluster_colors)+
        scale_y_reverse(breaks = integer_breaks(), limits = c(max_val,0)) +
        #ggtitle(paste0('Age: ', age_order[i], 'Y'))+
    guides(fill = F)+coord_flip()

    p2 <- ggplot(t_nt_df_2, aes(clonotype, fill = tissue_identity))+geom_bar(color = 'black', size = 0.3)+
        theme(axis.text.y=element_blank(),
             axis.line = element_blank(), axis.ticks = element_blank(),
             panel.background = element_rect(fill = 'white', colour = NA),
             panel.border = element_rect(fill = NA, colour = 'black'), 
             axis.text.x = element_text(size=13.75),
             axis.title = element_text(size=13.75),
             legend.text =  element_text(size=13.75),
             legend.title =  element_text(size=13.75),
             plot.title =  element_blank())+
        labs(x = NULL, y='Counts', fill = 'thymus_identity', title = 'thymus')+
        scale_x_discrete(limits = rev(names(stats)))+
        scale_fill_manual(values = combined_cluster_colors)+
        scale_y_continuous(breaks = integer_breaks(), limits = c(0, max_val)) +
    guides(fill = F)+coord_flip()
        #ggtitle(paste0('Age: ', age_order[i], 'Y'))+
    p <-wrap_plots(list(p1,p2))+plot_layout(ncol = 2,guides ='collect')+
    plot_annotation(title = paste0('Age: ', age_order[i], 'Y'))&
    theme(plot.margin = unit(rep(0,4), "lines"), plot.title = element_text(hjust = 0.5, size = 13.75))
    print(p)
    return(p)
}
)
names(clu_list) <- sample_list

###Fig.7e
sample_used2 <- c('20041504F-10X5','20041502M-10X5','20041307F-10X5',
'20041305F-10X5','20041604F-10X5','20041601M-10X5')
options(repr.plot.width=18, repr.plot.height=3)
p <- wrap_plots(clu_list[sample_used2])+plot_layout(ncol = 6)
p
ggsave('clonotype_shared.svg', width = 18, height = 3)


###Fig.7e legend
clu_selected <- unique(unlist(lapply(sample_used2, function(clu){
    clu1 <- unique(clu_list[[clu]]$data$tissue_identity)
    clu2 <- unique(clu_list[[clu]]$patches$plots[[1]]$data$tissue_identity)
    union(clu1, clu2)
})))
clu_selected <- intersect(combined_lvls,clu_selected)

lgd = Legend(labels = lbls_func(clu_selected), title = "Identity", nrow=2,by_row = T,title_position = 'leftcenter',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
              legend_gp = gpar(fill = combined_cluster_colors[clu_selected]))

pd = packLegend(list = list(lgd), row_gap = unit(1, "cm"))
options(repr.plot.width=13.75, repr.plot.height=1)
ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd)))
ggsave('thymus_PB_clonotype_shared_legend.svg', width =13.75, height = 1)

get_pairs <- function(char){
    char <- unique(char)
    result_x <- lapply(char, function(x){
        if(which(char == x) != length(char)){
            lapply(char[(which(char == x)+ 1):length(char)], function(y){
                return(c(x, y))
            })            
        }
    })
    result_x <- unlist(result_x, recursive = F)
    result_x <- result_x[!sapply(result_x, is.null)]
    return(result_x)
}


total_TRSS <- lapply(1:length(sample_list), function(i){
    nt_df <-  combined_clonotype_list[[i]]
    #nt_df$tissue_identity[nt_df$tissue_identity %in% c('CD8+RTE','CD8+TN')] <- 'CD8+TN'
    
    #nt_df$tissue_tissue_identity <- paste0(nt_df$tissue, '_',nt_df$tissue_identity)
    pairs <- get_pairs(unique(nt_df$tissue_identity))
    pairs <- c(pairs, lapply(pairs, function(x) rev(x)))
    clu_val <- sapply(pairs, function(pair){
        clu_1_df <- nt_df[nt_df$tissue_identity %in% pair[1], ]
        clu_2_df <- nt_df[nt_df$tissue_identity %in% pair[2], ]
        clu_1_df_2 <- clu_1_df[clu_1_df$clonotype %in% clu_2_df$clonotype,]
        clu_2_df_2 <- clu_2_df[clu_2_df$clonotype %in% clu_1_df$clonotype,]  
        clu_1_stats <- table(clu_1_df_2$clonotype)
        clu_2_stats <- table(clu_2_df_2$clonotype)    
        val <- sapply(names(clu_1_stats), function(clonotype){
            sqrt(clu_1_stats[clonotype]/nrow(clu_1_df) * clu_2_stats[clonotype]/nrow(clu_2_df))
        })
        val <- sum(unlist(val))              
    })    
    cluster1 <- sapply(pairs, function(x) x[1])
    cluster2 <- sapply(pairs, function(x) x[2])                   
    df <- data.frame(cluster1 = cluster1, cluster2 = cluster2, TRSS = clu_val, 
                     sample = sample_list[i], age = age_order[i])        
      

    return(df)
})
total_TRSS <- do.call(rbind, total_TRSS)

total_TRSS$cluster1 <- factor(total_TRSS$cluster1, levels = intersect(combined_lvls, total_TRSS$cluster1))

total_TRSS$cluster2 <- factor(total_TRSS$cluster2, levels = intersect(combined_lvls, total_TRSS$cluster2))

total_TRSS$class <- paste0(total_TRSS$cluster1, '~', total_TRSS$cluster2)
TRSS_sum <- ddply(total_TRSS, `.variables` = c('class'), summarise, TRSS = mean(TRSS))
TRSS_sum$cluster1 <- sapply(strsplit(TRSS_sum$class, split = '~'), function(x){x[1]})
TRSS_sum$cluster2 <- sapply(strsplit(TRSS_sum$class, split = '~'), function(x){x[2]})


total_TRSS_plot_fun <- function(clu, font_size = 15){
    clu <- as.character(clu)
    plot_df <- total_TRSS[total_TRSS$cluster1 %in% clu, ]
    plot_df$tissue <- sapply(strsplit(as.character(plot_df$cluster2), split = '_'), function(x){x[1]})
    plot_df$sample <- factor(plot_df$sample, levels = sample_list)
    p <- ggplot(plot_df, 
                aes(sample, cluster2, fill = TRSS, size = TRSS))+
    geom_point(shape = 21)+
    scale_fill_gradientn(colours = c('lightgrey', 'red'))+#viridis::viridis(20)
    scale_size_continuous(range = c(0,3.5))+

    guides(fill = guide_legend(shape = 19, color = 'black', reverse = T), 
          size = guide_legend(reverse = T))+
        labs(x = "Age")+
    facet_wrap(.~cluster1, labeller = function(x){
        x <- gsub("CD4.*Naive", "mature CD4 T", as.character(clu))
        x <- gsub("CD8.*Naive", "mature CD8 T", x)
        x <- gsub("CD4.*pre_T", "pre-CD4 T", x)
        x <- gsub("CD8.*pre_T", "pre-CD8 T", x)
        x <- gsub("CD4.*Memory_like", "mem-like CD4", x)
        x <- gsub("CD8.*Memory_like", "mem-like CD8", x)
        x <- gsub("thymus_", "Thy ", x)
        x <- gsub("PB_", "PB  ", x)
        x <- gsub("_", "-", x)
        x <- gsub("\\+", " ", x)
        
    })+    
    theme(text = element_text(size = font_size), 
              axis.title.x = element_text(size=font_size),axis.ticks = element_blank(), 
              axis.title.y = element_blank(),
             axis.text.x = element_text(size=font_size, hjust=1, angle = 90, vjust = 0.5), axis.text.y = element_blank(),
              axis.line  = element_blank(), 
             panel.border=element_rect(colour = "black", fill =NA),
              panel.background = element_rect(colour = NA, fill ='white'),
              legend.key.width = unit(2, "lines"),
           legend.key.height = unit(2, "lines"), legend.title = element_text(size=font_size),
              plot.title = element_text(size=font_size, hjust=0.5),
              legend.text = element_text(size=font_size),
          strip.text.x = element_text(size=font_size),
         strip.background.x=element_rect(fill = combined_cluster_colors[clu], color = 'black'))+
        scale_x_discrete(label = function(x){
            age_order[match(x, sample_list)]
        })
    left_label_1 <- ggplot(plot_df, aes(1, y=cluster2, fill=cluster2)) + geom_tile() +
      scale_fill_manual(values = combined_cluster_colors)+
    guides(fill=F)+
    theme_void()
    left_label_2 <- ggplot(plot_df, aes(1, y=cluster2, fill=tissue)) + geom_tile() +
      scale_fill_manual(values = c('PB' = '#E87D72', 'thymus' = '#56BCC2'))+
    guides(fill=F)+
    theme_void()
    new_p <- p %>% insert_left(left_label_1, width=0.05) %>% insert_left(left_label_2, width=0.05)
    #print(left_label_2)
    return(new_p)
    
}

options(repr.plot.width=4, repr.plot.height=5)
selected_cluster <- c('thymus_CD8+Memory_like', 'thymus_CD4+Memory_like',  'thymus_Treg', 
                     'PB_CD4+TCM', 'PB_CD4+TEFF', 'PB_CD4+TEM', 'PB_CD4+CTL', 'PB_Treg', 'PB_CD8+TCM', 
                     'PB_CD8+TEM-K','PB_CD8+TEM-B')

plot_list <- lapply(selected_cluster, function(cluster1){
    p <- print(total_TRSS_plot_fun(cluster1, font_size = 13.75))
    return(p)
})


names(plot_list) <- selected_cluster


###Extended Data Fig.8e
options(repr.plot.width=18.5, repr.plot.height=6)
p <- wrap_plots(plot_list)+plot_layout(ncol = 6)
p
ggsave(width = 18.5, height = 6, file = 'thymus_PB_shared_TCR2.svg')
#my_saveplot_fun(p, width = 18.5, height = 6, file = 'thymus_PB_shared_TCR2')

tmp_func <- function(x){
    x <- gsub("CD4.*Naive", "mature CD4 T", as.character(x))
    x <- gsub("CD8.*Naive", "mature CD8 T", x)
    x <- gsub("CD4.*pre_T", "pre-CD4 T", x)
    x <- gsub("CD8.*pre_T", "pre-CD8 T", x)
    x <- gsub("CD4.*Memory_like", "mem-like CD4", x)
    x <- gsub("CD8.*Memory_like", "mem-like CD8", x)
    x <- gsub("thymus_", "thymus ", x)
    x <- gsub("PB_", "PB  ", x)
    x <- gsub("_", "-", x)
    x <- gsub("\\+", " ", x)
    x <- gsub("thymus", "Thy", x)    
}
clu_selected <- levels(total_TRSS$cluster2)
lgd1 = Legend(labels = tmp_func(clu_selected), title = "Identity", nrow= 4,
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
              legend_gp = gpar(fill = combined_cluster_colors[clu_selected]))

lgd2 = Legend(labels = c('thymus', 'PB'), title = "Tissue", nrow = 2,
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
              legend_gp = gpar(fill = c('PB' = '#E87D72', 'thymus' = '#56BCC2')[c('thymus', 'PB')]))

pd = packLegend(list = list(lgd2, lgd1), column_gap = unit(0.5, "cm"), direction = 'horizontal')


###Extended Data Fig.8e legend
options(repr.plot.width=16, repr.plot.height=1.5)
ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd)))
ggsave('thymus_PB_TRSS_legend.svg', width =16, height = 1.5)

###Fig.7f
cluster_used <- c('thymus_CD8+Memory_like','PB_CD8+RTE',
'PB_CD8+TN','PB_CD8+TCM','PB_CD8+TEM-K','PB_CD8+TEM-B')
plot_df <- TRSS_sum[TRSS_sum$cluster1 %in% cluster_used & 
                          TRSS_sum$cluster2 %in% cluster_used, ]

min_val <- min(plot_df$TRSS)
max_val <- max(plot_df$TRSS)
tmp_breaks <- round(seq(min_val, max_val, length.out = 5),2)

p <- ggplot(plot_df, aes(cluster1, cluster2))+
geom_point(shape = 21, mapping = aes(fill = TRSS, size = TRSS))+
geom_tile(fill =NA, color = 'black', size = 0.1)+
scale_fill_gradientn(colours = c('white', '#db72b4'), breaks = tmp_breaks,#limits = rev,
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                    )+
scale_size_continuous(range = c(0,8.5), , breaks = tmp_breaks
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                     )+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(cluster_used))+
theme(#plot.margin = margin(b = 0.7,l = 0.7, unit = 'cm'),
      axis.title.x = element_text(size=13.75),axis.ticks = element_blank(),strip.background = element_blank(),
      legend.position = 'right',
          axis.title.y = element_blank(),
         axis.text.x = element_text(size=13.75, hjust=1, angle = 45), 
      axis.text.y = element_text(size=13.75),
          axis.line  = element_blank(), 
         panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'),
          legend.key.width = unit(0, "lines"),
       legend.key.height = unit(0, "lines"), legend.title = element_text(size=13.75),
          plot.title = element_text(size=13.75, hjust=0.5),
          legend.text = element_text(size=13.75),
      strip.text.x = element_text(size=13.75))+
labs(x = NULL, y = NULL)+
guides(fill = guide_legend(shape = 19, color = 'black', ncol = 1, reverse = T), 
      size = guide_legend(reverse = T))

options(repr.plot.width=5, repr.plot.height=3.5)
p
ggsave('CD8_total_TRSS.svg', width = 5, height = 3.5)
#my_saveplot_fun(p, 'CD8_total_TRSS', width = 5, height = 3.5)

###Extended Data Fig.8f
cluster_used <- c(
'thymus_CD4+Memory_like',
         'PB_CD4+RTE','PB_CD4+TN','PB_CD4+TCM','PB_CD4+TEFF','PB_CD4+TEM','PB_CD4+CTL')
plot_df <- TRSS_sum[TRSS_sum$cluster1 %in% cluster_used & 
                          TRSS_sum$cluster2 %in% cluster_used, ]

min_val <- min(plot_df$TRSS)
max_val <- max(plot_df$TRSS)
tmp_breaks <- round(seq(min_val, max_val, length.out = 5),2)

p <- ggplot(plot_df, aes(cluster1, cluster2))+
geom_point(shape = 21, mapping = aes(fill = TRSS, size = TRSS))+
geom_tile(fill =NA, color = 'black', size = 0.1)+
scale_fill_gradientn(colours = c('white', '#909ee7'), breaks = tmp_breaks
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                    )+
scale_size_continuous(range = c(0,8.5), , breaks = tmp_breaks
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                     )+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(cluster_used))+
theme(#plot.margin = margin(b = 0.7,l = 0.8, unit = 'cm'),
      axis.title.x = element_text(size=13.75),axis.ticks = element_blank(),strip.background = element_blank(),
      legend.position = 'right',
          axis.title.y = element_blank(),
         axis.text.x = element_text(size=13.75, hjust=1, angle = 45), 
      axis.text.y = element_text(size=13.75),
          axis.line  = element_blank(), 
         panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'),
          legend.key.width = unit(0, "lines"),
       legend.key.height = unit(0, "lines"), legend.title = element_text(size=13.75),
          plot.title = element_text(size=13.75, hjust=0.5),
          legend.text = element_text(size=13.75),
      strip.text.x = element_text(size=13.75))+
labs(x = NULL, y = NULL)+
guides(fill = guide_legend(shape = 19, color = 'black', ncol = 1, reverse = T), 
      size = guide_legend(reverse = T))

options(repr.plot.width=5, repr.plot.height=3.5)
p
ggsave('CD4_total_TRSS.svg', width = 5, height = 3.5)
#my_saveplot_fun(p, 'CD4_total_TRSS', width = 5, height = 3.5)



CD8_combined_T <- merge(T_dev_3[,grepl('CD8', T_dev_3$Identity)], 
                       Tcell[, grepl('CD8', Tcell$Identity)])

Idents(CD8_combined_T) <- CD8_combined_T$Identity



CD8_genes_mean <- rowMeans(expm1(CD8_combined_T$RNA@data))

CD8_exp <- AverageExpression(CD8_combined_T, features = names(CD8_genes_mean)[CD8_genes_mean > 0.5], 
                             verbose = F,assays = 'RNA')$RNA

CD8_cor_result <- cor(CD8_exp, method = 'spearman')


###Fig.7g
CD8_cor_result_plot_df <- as.data.frame(as.table(CD8_cor_result))
colnames(CD8_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- c(
'thymus_CD8+Memory_like',
         'PB_CD8+RTE','PB_CD8+TN','PB_CD8+TCM','PB_CD8+TEM-K','PB_CD8+TEM-B')


p <- ggplot(CD8_cor_result_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),na.value="lightgrey",
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         limits = round(range(unlist(CD8_cor_result)), 2))+
    theme(#plot.margin = margin(b = 0.7,l = 0.8, unit = 'cm'),
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
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=5, repr.plot.height=3.5)
p
ggsave('Spearman_cor_CD8_T.svg', width = 5, height = 3.5)
#my_saveplot_fun(p, 'Spearman_cor_CD8_T', width = 5, height = 3.5)

CD4_combined_T <- merge(T_dev_3[,grepl('CD4', T_dev_3$Identity)], 
                       Tcell[, grepl('CD4', Tcell$Identity)])

Idents(CD4_combined_T) <- CD4_combined_T$Identity



CD4_genes_mean <- rowMeans(expm1(CD4_combined_T$RNA@data))

CD4_exp <- AverageExpression(CD4_combined_T, features = names(CD4_genes_mean)[CD4_genes_mean > 0.5], 
                             verbose = F,assays = 'RNA')$RNA

CD4_cor_result <- cor(CD4_exp, method = 'spearman')



###Extended Data Fig.8g
CD4_cor_result_plot_df <- as.data.frame(as.table(CD4_cor_result))
colnames(CD4_cor_result_plot_df) <- c('cluster1', 'cluster2', 'Cor')

cluster_used <- cluster_used <- c(
'thymus_CD4+Memory_like',
         'PB_CD4+RTE','PB_CD4+TN','PB_CD4+TCM','PB_CD4+TEFF','PB_CD4+TEM','PB_CD4+CTL')


p <- ggplot(CD4_cor_result_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),na.value="lightgrey",
                         oob = scales::squish,#c("navy", "white", "firebrick3"), 
                         limits = round(range(unlist(CD4_cor_result)), 2))+
    theme(#plot.margin = margin(b = 0.7,l = 0.8, unit = 'cm'),
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
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(cluster_used))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(cluster_used))+
labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=5, repr.plot.height=3.5)
p
ggsave('Spearman_cor_CD4_T.svg', width = 5, height = 3.5)
#my_saveplot_fun(p, 'Spearman_cor_CD4_T', width = 5, height = 3.5)









###Fig.7h
CD8_combined_T_sub <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/CD8_combined_T_sub.rds')
levels(CD8_combined_T_sub) <- c('PB_CD8+TCM', 'thymus_CD8+Memory_like', 'PB_CD8+TEM-K')

CD8_combined_T_sub_markers <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/CD8_combined_T_sub_markers.rds')

top_markers <- lapply(unique(CD8_combined_T_sub_markers$cluster), function(clu){
    df <- CD8_combined_T_sub_markers[CD8_combined_T_sub_markers$p_val_adj < 0.05 &
                                     CD8_combined_T_sub_markers$avg_logFC > log(1.3) &
                                     CD8_combined_T_sub_markers$cluster == clu, ]
    df <- df[order(df$p_val_adj), ]
    genes <- df$gene[1:200]
    genes <- genes[!is.na(genes)]
    genes
})
top_markers <- unique(unlist(top_markers))
scale_exp <- my_scale(my_group_mean(expm1(CD8_combined_T_sub$RNA@data[top_markers, ]), Idents(CD8_combined_T_sub)),do_scale = 'row')
scale_exp <- order_heatmap2(scale_exp)
scale_exp <- t(scale_exp)
tmp_colors <- combined_cluster_colors[rownames(scale_exp)]
names(tmp_colors) <- rownames(scale_exp) <- 
gsub('Memory like', 'memory-like T',
     gsub('_', ' ', gsub('CD8\\+', '', gsub('thymus|PB', '', rownames(scale_exp)))))
plot_features <- c(
'ACTN1','MAL','LEF1','CCR7','TCF7','SELL','LTB','IL7R','NOSIP','KLF2',
    'CXCR6','JUN','FOS','ITGA1','GZMK','PDCD1','CXCR3','CCR4','ITGAE','CCR5','CD69','FGFBP2','CX3CR1',
    'GZMH','NKG7','KLRD1','FCGR3A','PRF1','CST7','CCL5','GZMA')
text_colors <- rev(
    tmp_colors[
        rownames(scale_exp)[apply(scale_exp[, plot_features,drop = F], 2, which.max)]
    ])
p <- Heatmap(scale_exp,cluster_rows=F, cluster_columns = F,show_column_names = F,
             column_names_rot = 90,column_names_centered = F,column_names_side = 'bottom',
             #row_names_gp = gpar(fontsize = 0, fontfamily = 'sans'),
             #row_names_side = c("left"), 
             #width = unit(13.75, 'cm'),
             col = circlize::colorRamp2(c(-1.2, 0, 1.2), c('#00add9', '#040600', '#fdfc00')),
             rect_gp = gpar(col = NA),
             heatmap_legend_param = list(
                 title_position = c('leftcenter'),
               direction = "horizontal",                 
               title = 'Expression',
               title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
               labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans')),
             show_heatmap_legend=T, show_row_names = F, 
            bottom_annotation =  HeatmapAnnotation(' ' = anno_mark(link_height = unit(7.5,'mm'),side = 'bottom',
                                                                          #at_adj = -0.5,
                                                                          #extra_text_args = list(vjust = 0.5),
            at = which(colnames(scale_exp) %in% plot_features),
            labels = colnames(scale_exp)[colnames(scale_exp) %in% plot_features], 
            labels_gp = gpar(fontsize=13.75,col = text_colors, fontface='italic')), 
           'foo' = anno_empty(border = FALSE, 
                              height = max_text_height(colnames(scale_exp)[colnames(scale_exp) %in% plot_features]))                           
                                         ) ,        
             left_annotation = rowAnnotation(
                 #annotation_label = 
               df=data.frame(Identity = rownames(scale_exp)), 
                 col = list(Identity = tmp_colors[rownames(scale_exp)]),
               show_annotation_name=F,
               annotation_legend_param = list(
                 Identity = list(
                      title_position = c('leftcenter'),
                   #legend_direction = "horizontal",
                   #title = "Baaaaaaar",
                   at = rownames(scale_exp),
                   title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
                   nrow = 1
                 ))
             )
            )
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p,merge_legend = T,  heatmap_legend_side="top",ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
options(repr.plot.width=7.5, repr.plot.height=3.5)
gg_p
ggsave('CD8_combined_T_sub_DEGs.svg', width = 7.5, height = 3.5)
#my_saveplot_fun(gg_p,'CD8_combined_T_sub_DEGs', width = 7.5, height = 3.5)



###Extended Data Fig.7j left
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TCM.Gsea.1687096756411/gsea_report_for_thymus_CD8_Memory_like_1687096756411.tsv.xlsx')
data1$cluster <- 'CD8+Memory_like'





data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TCM.Gsea.1687096756411/gsea_report_for_PB_CD8_TCM_1687096756411.tsv.xlsx')
data2$cluster <- 'CD8+TCM'

data3 <- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TCM.Gsea.1687099406217/gsea_report_for_thymus_CD8_Memory_like_1687099406217.tsv.xlsx')
data3$cluster <- 'CD8+Memory_like'


CD8_memory_GSEA_df <- rbind(data1, data2, data3)


pathways <- c('TRM_SIGNATURE','GO_T_CELL_CHEMOTAXIS','GO_T_CELL_MEDIATED_CYTOTOXICITY',
'GO_REGULATION_OF_T_CELL_CYTOKINE_PRODUCTION','PHOSPHORYLATION OF CD3 AND TCR ZETA CHAINS',
'PD-1 SIGNALING','GO_POSITIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION',
'GO_MHC_PROTEIN_COMPLEX','GO_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION',
'GO_T_CELL_RECEPTOR_COMPLEX','GO_RIBOSOME','GO_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS',
'TRANSLATION','GO_RESPONSE_TO_INTERLEUKIN_7', 
             'GO_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM')

CD8_memory_GSEA_df$NES <- -1*CD8_memory_GSEA_df$NES
CD8_memory_GSEA_df <- CD8_memory_GSEA_df[CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
CD8_memory_GSEA_df <- CD8_memory_GSEA_df[order(CD8_memory_GSEA_df$NES, decreasing = T), ]
CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB <- factor(CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(CD8_memory_GSEA_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color = 'black')+
scale_y_discrete(label = function(x){
    x <- gsub('TRM_SIGNATURE', 'TISSUE RESIDENT MEMORY T CELL SIGNATURE', x)
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)
    x <- str_to_sentence_func(x)
    x <- gsub('INTERFERON GAMMA', 'IFN-γ', x, ignore.case = T)
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), 
      plot.margin = margin(l = 1, unit = 'cm'),
      legend.box.margin=margin(b = -0.3, unit = 'cm'),
      axis.ticks=element_blank(),#legend.box = 'vertical',
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c(cluster_colors),
                  label = c('Thy memory-like T', 'PB TCM'))+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=5.6, repr.plot.height=3.5)
p
ggsave('thymus_CD8_Memory_VS_PB_TCM_GSEA.svg', width = 9, height = 4)
#my_saveplot_fun(p, 'thymus_CD8_Memory_VS_PB_TCM_GSEA', width = 9, height = 4)




###Extended Data Fig.7j right
data1<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TEM.Gsea.1687081667359/gsea_report_for_thymus_CD8_Memory_like_1687081667359.tsv.xlsx')
data1$cluster <- 'CD8+Memory_like'
data2<- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TEM.Gsea.1687081667359/gsea_report_for_PB_CD8_TEM_1687081667359.tsv.xlsx')
data2$cluster <- 'CD8+TEM-K'

data3 <- read.xlsx('/data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/jun18/thymus_CD8_Memory_like_vs_PB_CD8_TEM.Gsea.1687089341786/gsea_report_for_thymus_CD8_Memory_like_1687089341786.tsv.xlsx')
data3$cluster <- 'CD8+Memory_like'

CD8_memory_GSEA_df <- rbind(data1, data2, data3)


pathways <- c('TRM_SIGNATURE',
    'GO_T_CELL_CHEMOTAXIS',
'GO_ACTIVATED_T_CELL_PROLIFERATION',
              'GO_TOLERANCE_INDUCTION',
'GO_NEGATIVE_REGULATION_OF_T_CELL_APOPTOTIC_PROCESS',
'FOXO-MEDIATED TRANSCRIPTION OF CELL CYCLE GENES',
'GO_T_CELL_HOMEOSTASIS',
'GO_NEGATIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION',
'REGULATION OF RUNX3 EXPRESSION AND ACTIVITY',
#'TCF DEPENDENT SIGNALING IN RESPONSE TO WNT',
'APOPTOSIS',
'PHOSPHORYLATION OF CD3 AND TCR ZETA CHAINS',
'GO_MHC_PROTEIN_COMPLEX',
'INTERFERON GAMMA SIGNALING',
'GO_T_CELL_MIGRATION',
'GO_LEUKOCYTE_MEDIATED_CYTOTOXICITY',
             'RRNA PROCESSING IN THE NUCLEUS AND CYTOSOL')

CD8_memory_GSEA_df$NES <- -1*CD8_memory_GSEA_df$NES
CD8_memory_GSEA_df <- CD8_memory_GSEA_df[CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB %in% pathways, ]
CD8_memory_GSEA_df <- CD8_memory_GSEA_df[order(CD8_memory_GSEA_df$NES, decreasing = T), ]
CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB <- factor(CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB, 
                                                         levels = CD8_memory_GSEA_df$GS.br..follow.link.to.MSigDB)


p <- ggplot(CD8_memory_GSEA_df, aes(NES, `GS.br..follow.link.to.MSigDB`, fill = cluster))+
geom_bar(stat='identity', color = 'black')+
scale_y_discrete(label = function(x){
    x <- gsub('TRM_SIGNATURE', 'TISSUE RESIDENT MEMORY T CELL SIGNATURE', x)
    x <- gsub("KEGG_|GO_", '', x)
    x <- gsub('_', ' ', x)
    x <- str_to_sentence_func(x) 
    x <- gsub('INTERFERON GAMMA', 'IFN-γ', x, ignore.case = T)
    
},position = 'right')+
theme(text = element_text(size = 13.75, family = 'sans'), 
      plot.margin = margin(l = 1, unit = 'cm'),
      legend.box.margin=margin(b = -0.4, unit = 'cm'),
      axis.ticks=element_blank(),#legend.box = 'vertical',
      axis.text=element_text(size = 13.75, family = 'sans'), 
      legend.text=element_text(size = 13.75, family = 'sans'), 
      axis.title=element_text(size = 13.75, family = 'sans'), 
      panel.background = element_rect(color = NA, fill ='white'),
      axis.line =element_blank(),
      legend.position='top',
      panel.border = element_rect(color = 'black', fill =NA),)+

scale_fill_manual(values = c(cluster_colors),
                  label = c('Thy memory-like T', 'PB TEM-K'))+
labs(y=NULL,fill=NULL)+
guides(fill = guide_legend(ncol = 1))
options(repr.plot.width=5.6, repr.plot.height=3.5)
p
ggsave('thymus_PB_CD8_Memory_GSEA.svg', width = 9, height = 4)
#my_saveplot_fun(p, 'thymus_PB_CD8_Memory_GSEA', width = 9, height = 4)


total_TRSS$age <- as.numeric(total_TRSS$age)

total_TRSS$age_3 <- ''
total_TRSS$age_3[total_TRSS$age <=12] <- 'Group I'
total_TRSS$age_3[total_TRSS$age >= 13 & total_TRSS$age <= 39 ] <- 'Group II'
total_TRSS$age_3[total_TRSS$age >=40] <- 'Group III'

TRSS_by_age_group <- ddply(total_TRSS, `.variables` = c('class', 'age_3'), summarise, TRSS = mean(TRSS))
TRSS_by_age_group$cluster1 <- sapply(strsplit(TRSS_by_age_group$class, split = '~'), function(x){x[1]})
TRSS_by_age_group$cluster2 <- sapply(strsplit(TRSS_by_age_group$class, split = '~'), function(x){x[2]})

###Fig.7i
cluster_used <- c('thymus_CD8+Memory_like', 'PB_CD8+TCM', 
                     'PB_CD8+TEM-K','PB_CD8+TEM-B')
plot_df <- TRSS_by_age_group[TRSS_by_age_group$cluster1 %in% cluster_used & 
                          TRSS_by_age_group$cluster2 %in% cluster_used, ]
plot_df$age_3 <-factor(plot_df$age_3, levels = c('Group I', 'Group II', 'Group III'))

min_val <- min(plot_df$TRSS)
max_val <- max(plot_df$TRSS)
tmp_breaks <- round(seq(min_val, max_val, length.out = 5),2)

p <- ggplot(plot_df, aes(cluster1, cluster2))+
geom_point(shape = 21, mapping = aes(fill = TRSS, size = TRSS))+
geom_tile(fill =NA, color = 'black', size = 0.1)+
scale_fill_gradientn(colours = c('white', '#db72b4'), 
                     breaks = tmp_breaks)+
scale_size_continuous(range = c(0,10), 
                     breaks = tmp_breaks)+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(c('thymus_CD8+Memory_like', 'PB_CD8+TCM', 
                     'PB_CD8+TEM-K','PB_CD8+TEM-B')))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(c('thymus_CD8+Memory_like', 'PB_CD8+TCM', 
                     'PB_CD8+TEM-K','PB_CD8+TEM-B')))+
theme(#plot.margin = margin(b = 0.6,l = 0.8, unit = 'cm'),
      axis.title.x = element_text(size=13.75),axis.ticks = element_blank(),strip.background = element_blank(),
      legend.position = 'right',
          axis.title.y = element_blank(),
         axis.text.x = element_text(size=13.75, hjust=1, angle = 45), 
      axis.text.y = element_text(size=13.75),
          axis.line  = element_blank(), 
         panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'),
          legend.key.width = unit(2, "lines"),
       legend.key.height = unit(2, "lines"), legend.title = element_text(size=13.75),
          plot.title = element_text(size=13.75, hjust=0.5),
          legend.text = element_text(size=13.75),
      strip.text.x = element_text(size=13.75))+
labs(x = NULL, y = NULL)+
facet_wrap(.~age_3, ncol = 3)+
guides(fill = guide_legend(shape = 19, color = 'black', ncol = 1, reverse = T), 
      size = guide_legend(reverse = T))



# coord_cartesian(ylim = c(1, 4),expand = 1, clip="off", xlim = c(1, 4)) 
options(repr.plot.width=9, repr.plot.height=3.5)
p
ggsave('CD8_total_TRSS_by_age.svg', width = 9, height = 3.5)
#my_saveplot_fun(p, 'CD8_total_TRSS_by_age', width = 9, height = 3.5)


###Extended Data Fig.8h
cluster_used <-  c('thymus_CD4+Memory_like', 
                     'PB_CD4+TCM', 'PB_CD4+TEFF', 'PB_CD4+TEM', 'PB_CD4+CTL')
plot_df <- TRSS_by_age_group[TRSS_by_age_group$cluster1 %in% cluster_used & 
                          TRSS_by_age_group$cluster2 %in% cluster_used, ]
plot_df$age_3 <-factor(plot_df$age_3, levels = c('Group I', 'Group II', 'Group III'))

min_val <- min(plot_df$TRSS)
max_val <- max(plot_df$TRSS)
tmp_breaks <- round(seq(min_val, max_val, length.out = 5),2)


p <- ggplot(plot_df, aes(cluster1, cluster2))+
geom_point(shape = 21, mapping = aes(fill = TRSS, size = TRSS))+
geom_tile(fill =NA, color = 'black', size = 0.1)+
scale_fill_gradientn(colours = c('white', '#909ee7'), 
                     breaks = tmp_breaks)+
scale_size_continuous(range = c(0,10), 
                     breaks = tmp_breaks)+
scale_x_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                 limits = rev(c('thymus_CD4+Memory_like','PB_CD4+TCM', 'PB_CD4+TEFF', 'PB_CD4+TEM', 'PB_CD4+CTL')))+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
}, 
                limits = rev(c('thymus_CD4+Memory_like','PB_CD4+TCM', 'PB_CD4+TEFF', 'PB_CD4+TEM', 'PB_CD4+CTL')))+
theme(#plot.margin = margin(b = 0.6,l = 0.8 , unit = 'cm'),
      axis.title.x = element_text(size=13.75),axis.ticks = element_blank(),strip.background = element_blank(),
      legend.position = 'right',
          axis.title.y = element_blank(),
         axis.text.x = element_text(size=13.75, hjust=1, angle = 45), 
      axis.text.y = element_text(size=13.75),
          axis.line  = element_blank(), 
         panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'),
          legend.key.width = unit(2, "lines"),
       legend.key.height = unit(2, "lines"), legend.title = element_text(size=13.75),
          plot.title = element_text(size=13.75, hjust=0.5),
          legend.text = element_text(size=13.75),
      strip.text.x = element_text(size=13.75))+
labs(x = NULL, y = NULL)+
facet_wrap(.~age_3, ncol = 3)+
guides(fill = guide_legend(shape = 19, color = 'black', ncol = 1, reverse = T), 
      size = guide_legend(reverse = T))


options(repr.plot.width=9, repr.plot.height=3.5)
p
ggsave('CD4_total_TRSS_by_age.svg', width = 9, height = 3.5)
#my_saveplot_fun(p, 'CD4_total_TRSS_by_age', width = 9, height = 3.5)




CD8_memory <- merge(T_dev_3[,grepl('CD8', T_dev_3$Identity) & !grepl('Naive|pre_T', T_dev_3$Identity)], 
                       Tcell[, grepl('CD8', Tcell$Identity) & !grepl('TN|RTE', Tcell$Identity)])

Idents(CD8_memory) <- CD8_memory$Identity

levels(CD8_memory) <- c('thymus_CD8+Memory_like', 'PB_CD8+TCM', 'PB_CD8+TEM-K','PB_CD8+TEM-B')

CD8_memory_genes_mean <- rowMeans(expm1(CD8_memory$RNA@data))



CD8_age_cor_result <- lapply(unique(CD8_memory$age_3), function(age){
    CD8_exp <- AverageExpression(CD8_memory[, CD8_memory$age_3 == age], 
                                 features = names(CD8_memory_genes_mean)[CD8_memory_genes_mean > 0.5], 
                                 verbose = F,assays = 'RNA')$RNA   
    result <- cor(CD8_exp, method = 'spearman')
    result
})#, cor(CD8_exp, method = 'spearman')

names(CD8_age_cor_result) <- unique(CD8_memory$age_3)

###Fig.7j
CD8_age_cor_plot_df <- lapply(names(CD8_age_cor_result), function(i){
    result <- CD8_age_cor_result[[i]]
    colnames(result) <- rownames(result) <- gsub('thymus_|PB_', '',colnames(result))
    colnames(result) <- rownames(result) <- gsub('CD.\\+', '',colnames(result))
    result <- result[rev(rownames(result)), rev(colnames(result))]
    mapping <- c('<=12' = 'Group I', '13-39' = 'Group II', '40-99' = 'Group III', '>=100' = 'Centenarian')
    plot_df <- as.data.frame(as.table(result))
    colnames(plot_df) <- c('cluster1', 'cluster2', 'Cor')
    plot_df$age <- mapping[[i]]
    plot_df
})
CD8_age_cor_plot_df <- do.call(rbind, CD8_age_cor_plot_df)
CD8_age_cor_plot_df <- CD8_age_cor_plot_df[!CD8_age_cor_plot_df$age %in% c('Centenarian'), ]
CD8_age_cor_plot_df$age <- factor(CD8_age_cor_plot_df$age, levels = c('Group I', 'Group II', 'Group III'))



p <- ggplot(CD8_age_cor_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
facet_wrap(.~age)+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),#c('white', '#db72b4'),
                         oob = scales::squish,
                         na.value = 'lightgrey',#c("navy", "white", "firebrick3"), 
                         limits = round(range(unlist(CD8_age_cor_result)), 2))+
    theme(#plot.margin = margin(b = 0.6,l = 0.8, unit = 'cm'),
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
    x <- gsub('Memory_like', 'memory-like T', x)
})+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
})+


labs(fill  = 'Spearman \nCorrelation')

options(repr.plot.width=9, repr.plot.height=3.5)
p
ggsave('Spearman_cor_CD8_T_age.svg', width = 9, height = 3.5)
#my_saveplot_fun(p, 'Spearman_cor_CD8_T_age', width = 9, height = 3.5)



CD4_memory <- merge(T_dev_3[,grepl('CD4', T_dev_3$Identity) & !grepl('Naive|pre_T', T_dev_3$Identity)], 
                       Tcell[, grepl('CD4', Tcell$Identity) & !grepl('TN|RTE', Tcell$Identity)])

Idents(CD4_memory) <- CD4_memory$Identity



levels(CD4_memory) <- c('thymus_CD4+Memory_like', 
                        'PB_CD4+TCM', 'PB_CD4+TEM','PB_CD4+TEFF', 'PB_CD4+CTL')

CD4_memory_genes_mean <- rowMeans(expm1(CD4_memory$RNA@data))

CD4_age_cor_result <- lapply(unique(CD4_memory$age_3), function(age){
    CD4_exp <- AverageExpression(CD4_memory[, CD4_memory$age_3 == age], 
                                 features = names(CD4_memory_genes_mean)[CD4_memory_genes_mean > 0.5], 
                                 verbose = F,assays = 'RNA')$RNA   
    result <- cor(CD4_exp, method = 'spearman')
    result
})#, cor(CD4_exp, method = 'spearman')

names(CD4_age_cor_result) <- unique(CD4_memory$age_3)

###Extended Data Fig.8i
CD4_age_cor_plot_df <- lapply(names(CD4_age_cor_result), function(i){
    result <- CD4_age_cor_result[[i]]
    colnames(result) <- rownames(result) <- gsub('thymus_|PB_', '',colnames(result))
    colnames(result) <- rownames(result) <- gsub('CD.\\+', '',colnames(result))
    result <- result[rev(rownames(result)), rev(colnames(result))]
    mapping <- c('<=12' = 'Group I', '13-39' = 'Group II', '40-99' = 'Group III', '>=100' = 'Centenarian')
    plot_df <- as.data.frame(as.table(result))
    colnames(plot_df) <- c('cluster1', 'cluster2', 'Cor')
    plot_df$age <- mapping[[i]]
    plot_df
})
CD4_age_cor_plot_df <- do.call(rbind, CD4_age_cor_plot_df)
CD4_age_cor_plot_df <- CD4_age_cor_plot_df[!CD4_age_cor_plot_df$age %in% c('Centenarian'), ]
CD4_age_cor_plot_df$age <- factor(CD4_age_cor_plot_df$age, levels = c('Group I', 'Group II', 'Group III'))



p <- ggplot(CD4_age_cor_plot_df, aes(cluster1, cluster2))+
    geom_tile(color = 'black', mapping = aes(fill = Cor))+
facet_wrap(.~age)+
    scale_fill_gradientn(colours = c("midnightblue","dodgerblue3","goldenrod1","darkorange2"),#c('white', '#db72b4'),
                         oob = scales::squish,
                         na.value = 'lightgrey',#c("navy", "white", "firebrick3"), 
                         limits = round(range(unlist(CD4_age_cor_result)), 2))+
    theme(#plot.margin = margin(b = 0.6,l = 0.8, unit = 'cm'),
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
    x <- gsub('Memory_like', 'memory-like T', x)
})+
scale_y_discrete(label = function(x){
    x <- gsub('thymus_|PB_', '', x)
    x <- gsub('CD.\\+', '', x)
    x <- gsub('Memory_like', 'memory-like T', x)
})+
labs(fill  = 'Spearman \nCorrelation')



options(repr.plot.width=9, repr.plot.height=3.5)
p
ggsave('Spearman_cor_CD4_T_age.svg', width = 9, height = 3.5)
#my_saveplot_fun(p, 'Spearman_cor_CD4_T_age', width = 9, height = 3.5)



TCR_DB <- read.xlsx('/data1/02.private/dengyj/analysis/thymus/TCR/TCR_database/TCR-repertoire-database-main/TCR-db.xlsx')

PB_tcr_df <- do.call(rbind, PB_tcr$data)


PB_match_TRB_df<-  PB_tcr_df[
    PB_tcr_df$cdr3 %in% TCR_DB$cdr3.beta & PB_tcr_df$chain == 'TRB',
    c('barcode', 'cdr3', 'identity')]
colnames(PB_match_TRB_df)[colnames(PB_match_TRB_df)== 'cdr3'] <- 'cdr3.beta'
PB_match_TRB_df <- lapply(1:nrow(PB_match_TRB_df), function(row){
    tmp_df <- PB_match_TRB_df[row, ]
    ID <- TCR_DB[TCR_DB$cdr3.beta %in% tmp_df$cdr3.beta, 'ID']
    new_df <- data.frame(barcode = tmp_df$barcode,  cdr3.beta = tmp_df$cdr3.beta, 
                         identity = tmp_df$identity,
                        TCR_DB_ID = ID)
    new_df
})

PB_match_TRB_df <- do.call(rbind, PB_match_TRB_df)


PB_match_TRB_df$antigen.species <- TCR_DB[match(PB_match_TRB_df$TCR_DB_ID, TCR_DB$ID), 'antigen.species']

antigen.species <- unique(PB_match_TRB_df[!is.na(PB_match_TRB_df$cdr3.beta),'antigen.species'])
PB_clone_cell_df <- 
lapply(antigen.species, 
       function(antigen.specie){
           df <- PB_match_TRB_df[!is.na(PB_match_TRB_df$cdr3.beta) & 
                                   PB_match_TRB_df$antigen.species %in% antigen.specie,]
           df <- df[!duplicated(df$barcode), c('barcode', 'identity')]
           setdiff_df <- data.frame(barcode = setdiff(colnames(Tcell),df$barcode), identity = 'None')
           combined_df <- rbind(df, setdiff_df)
           rownames(combined_df) <- combined_df$barcode
           combined_df <- combined_df[colnames(Tcell), ]
           combined_df <- combined_df[, c('identity'), drop = F]
           colnames(combined_df)[colnames(combined_df) == 'identity'] <- antigen.specie
           combined_df
       })

PB_clone_cell_df <- do.call(cbind, PB_clone_cell_df)

Tcell <- AddMetaData(Tcell, metadata = PB_clone_cell_df[, 'CMV', drop = F])

options(repr.plot.width=16, repr.plot.height=16)
plot_df <- cbind(FetchData(Tcell, vars = c('CMV', 'age_3')), 
                 Embeddings(Tcell, reduction = 'umap'))
plot_df$age_3 <- factor(plot_df$age_3, levels = c('<=12', '13-39', '40-99', '>=100'))
levels(plot_df$age_3)[levels(plot_df$age_3) == '<=12'] <- 'Group I'
levels(plot_df$age_3)[levels(plot_df$age_3) == '13-39'] <- 'Group II'
levels(plot_df$age_3)[levels(plot_df$age_3) == '40-99'] <- 'Group III'
levels(plot_df$age_3)[levels(plot_df$age_3) == '>=100'] <- 'Centenarian'


###Fig.7K left
idents_order <-c(PB_T_idents, 'None')
antigen.specie <- 'CMV'

p <- ggplot(mapping = aes_string('UMAP_1', 'UMAP_2', color = antigen.specie))+
geom_point(data = plot_df[plot_df[, antigen.specie] == 'None', ], size = 0.1)+
geom_point(data = plot_df[plot_df[, antigen.specie] != 'None', ], size = 1.5)+
facet_wrap(.~age_3, nrow =1)+
scale_color_manual(values =cluster_colors, 
                   limits = intersect(idents_order, unique(plot_df[, antigen.specie])), 
                  label = lbls_func)  +
theme(panel.background = element_rect(fill = 'white', color = NA), 
      strip.background.x = element_rect(fill = 'white'),
      panel.border = element_rect(fill = NA, color = 'black'), 
      strip.text.x = element_text(size=20),
      axis.title = element_text(size=20),axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=20,hjust=0.5),
     legend.title = element_text(size = 20),legend.text = element_text(size = 20), 
     plot.margin = unit(c(0.5, 0, 0, 0), 'cm'))+
labs(title = antigen.specie, color = NULL)+
guides(color = guide_legend(override.aes = list(size= 5))) +
theme(text=element_text(color='black'),
     strip.text=element_text(color='black'), legend.text=element_text(color='black'))
options(repr.plot.width=18.5, repr.plot.height=5)
p
my_saveplot_fun(p, paste0('PB_', antigen.specie, '_TCR'), width = 18.5, height = 5)


###Fig.7K right
options(repr.plot.width=5, repr.plot.height=3.5)
plot_df <- Tcell@meta.data[, c('CMV', 'batch', 'age_3'), 'drop'  = F]
plot_df$barcode <- rownames(plot_df)
plot_df$CMV_count <- ifelse(!plot_df$CMV %in% c('CD8+TCM', 'CD4+TEM', 'CD8+TEM-K','CD4+CTL', 'CD8+TEM-B'), 0, 1)

plot_df <- ddply(plot_df, c('batch'), summarise,  CMV_count = sum(CMV_count), age_3 = unique(age_3))


plot_df_sum <- ddply(PB_tcr_df[!duplicated(PB_tcr_df$barcode),], 'Sample',summarise, count = length(Sample))
colnames(plot_df_sum)[colnames(plot_df_sum) =='Sample'] <- 'batch'
plot_df <- merge(plot_df, plot_df_sum)
plot_df$Freq <- plot_df$CMV_count/plot_df$count * 100
plot_df$age_3 <- factor(plot_df$age_3, levels =rev(c('<=12', '13-39', '40-99', '>=100')))

stat_test <- wilcox_test(
 Freq ~ age_3, data = plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.3,fun = 'mean')
stat_test$y.position <- seq(min(stat_test$y.position), max(stat_test$y.position), 
                            length.out = length(stat_test$y.position))+0.005

p1 <- ggplot(plot_df, aes(age_3, Freq))+
geom_boxplot(mapping = aes(fill = age_3))+
geom_point(position = 'jitter')+
guides(fill = F)+
scale_fill_manual(values = age_colors)+
    theme(panel.background = element_rect(fill = 'white', color = NA),plot.margin = margin(r = 1, unit = 'mm'), 
          panel.border = element_rect(fill = NA, color = 'black'), 
          axis.title.x = element_text(size=13.75),axis.text.x = element_text(size=13.75, angle = 45, hjust = 1),
          axis.title.y = element_text(size=13.75, angle = 0, vjust = 0.5),
          axis.ticks = element_blank(),axis.text.y =  element_text(size=13.75),
          plot.title = element_text(size=13.75,hjust=0.5),
          #axis.title.y = element_blank(),
         legend.title = element_text(size = 13.75),legend.text = element_text(size = 13.75), )+
labs(y = 'CMV-specific TCR in effector/memory cells %', x = NULL, fill = 'Age')+
coord_flip()+scale_x_discrete(labels =  age_group_lbls_func)+
scale_y_continuous(expand = expansion(mult = c(0.02, 0.035)))+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif', hjust = 0.5,coord.flip= T,
                   size = 13.75*15/13/.pt)

p1
saveRDS(plot_df, file = 'exp_T_PB_CMV_count_age.rds')
#ggsave(file = 'exp_T_PB_CMV_count_age.svg', width = 5, height = 3.5)
#my_saveplot_fun(p, file = 'exp_T_PB_CMV_count_age', width = 5, height = 3.5)





setwd('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/data/')
files <- dir()
files <- files[-grep('metadata.txt', files)]
thymus_tcr2 <- list()
thymus_tcr2$meta <- read.table('metadata.txt', header = T)  
thymus_tcr2$data <- 
lapply(files, function(i){
    df_iter <- read.csv(i)
    df_iter$Sample <- strsplit(i, split = '\\.')[[1]][1]
    df_iter$age <- thymus_tcr2$meta[match(df_iter$Sample, thymus_tcr2$meta$Sample), 'age']
    df_iter$age_3 <- thymus_tcr2$meta[match(df_iter$Sample, thymus_tcr2$meta$Sample), 'age_3']
    df_iter$identity <- thymus_identity[match(df_iter$barcode, names(thymus_identity))]
    df_iter$tissue <- 'thymus'
    df_iter$tissue_identity <- paste0( df_iter$tissue, '_', df_iter$identity)
    ####contain DP cells 
    df_iter <- df_iter[df_iter$identity %in% thymus_T_idents2, ]

    df_iter <- df_iter[df_iter$productive == 'True', ]
    tbl=table(df_iter$barcode, df_iter$chain)
    logi <- apply(tbl, 1, function(x){
        x['TRA'] > 0 & x['TRB'] > 0 
    })  
    df_iter <- df_iter[df_iter$barcode %in% rownames(tbl)[logi], ]
    return(df_iter)
})
names(thymus_tcr2$data) <- sapply(strsplit(files, split = '\\.'), function(x) x[1])                          
setwd(plot_dir)                          

thymus_tcr2$meta$age_3 <- factor(thymus_tcr2$meta$age_3, levels =  c('<=12', '13-39', '>=40'))

thymus_tcr_df <- do.call(rbind, thymus_tcr2$data)


thymus_match_TRB_df<-  thymus_tcr_df[
    thymus_tcr_df$cdr3 %in% TCR_DB$cdr3.beta & thymus_tcr_df$chain == 'TRB',
    c('barcode', 'cdr3', 'identity')]
colnames(thymus_match_TRB_df)[colnames(thymus_match_TRB_df)== 'cdr3'] <- 'cdr3.beta'
thymus_match_TRB_df <- lapply(1:nrow(thymus_match_TRB_df), function(row){
    tmp_df <- thymus_match_TRB_df[row, ]
    ID <- TCR_DB[TCR_DB$cdr3.beta %in% tmp_df$cdr3.beta, 'ID']
    new_df <- data.frame(barcode = tmp_df$barcode,  cdr3.beta = tmp_df$cdr3.beta, 
                         identity = tmp_df$identity,
                        TCR_DB_ID = ID)
    new_df
})

thymus_match_TRB_df <- do.call(rbind, thymus_match_TRB_df)


thymus_match_TRB_df$antigen.species <- TCR_DB[match(thymus_match_TRB_df$TCR_DB_ID, TCR_DB$ID), 'antigen.species']

antigen.species <- unique(thymus_match_TRB_df[!is.na(thymus_match_TRB_df$cdr3.beta),'antigen.species'])
thymus_clone_cell_df <- 
lapply(antigen.species, 
       function(antigen.specie){
           df <- thymus_match_TRB_df[!is.na(thymus_match_TRB_df$cdr3.beta) & 
                                   thymus_match_TRB_df$antigen.species %in% antigen.specie,]
           df <- df[!duplicated(df$barcode), c('barcode', 'identity')]
           setdiff_df <- data.frame(barcode = setdiff(colnames(T_dev_3),df$barcode), identity = 'None')
           combined_df <- rbind(df, setdiff_df)
           rownames(combined_df) <- combined_df$barcode
           combined_df <- combined_df[colnames(T_dev_3), ]
           combined_df <- combined_df[, c('identity'), drop = F]
           colnames(combined_df)[colnames(combined_df) == 'identity'] <- antigen.specie
           combined_df
       })

thymus_clone_cell_df <- do.call(cbind, thymus_clone_cell_df)

T_dev_3 <- AddMetaData(T_dev_3, metadata = thymus_clone_cell_df[, 'CMV', drop = F])

options(repr.plot.width=16, repr.plot.height=16)
plot_df <- cbind(FetchData(T_dev_3, vars = c('CMV', 'age_3')), 
                 Embeddings(T_dev_3, reduction = 'umap'))
plot_df$age_3 <- factor(plot_df$age_3, levels = c('<=12', '13-39', '40-99'))
levels(plot_df$age_3)[levels(plot_df$age_3) == '<=12'] <- 'Group I'
levels(plot_df$age_3)[levels(plot_df$age_3) == '13-39'] <- 'Group II'
levels(plot_df$age_3)[levels(plot_df$age_3) == '40-99'] <- 'Group III'

###Fig.7l left
idents_order <- c(thymus_T_idents2, 'None')
antigen.specie <- 'CMV'

p <- ggplot(mapping = aes_string('UMAP_1', 'UMAP_2', color = antigen.specie))+
geom_point(data = plot_df[plot_df[, antigen.specie] == 'None', ], size = 0.1)+
geom_point(data = plot_df[plot_df[, antigen.specie] != 'None', ], size = 1.5)+
facet_wrap(.~age_3, nrow =1)+
scale_color_manual(values =cluster_colors, 
                   limits = intersect(idents_order, unique(plot_df[, antigen.specie])), 
                  label = lbls_func)  +
theme(panel.background = element_rect(fill = 'white', color = NA),
      strip.background.x = element_rect(fill = 'white'),
      panel.border = element_rect(fill = NA, color = 'black'), 
      strip.text.x = element_text(size=20),
      axis.title = element_text(size=20),axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=20,hjust=0.5),
     legend.title = element_text(size = 20),legend.text = element_text(size = 20), 
     plot.margin = unit(c(0.5, 0, 0, 0), 'cm'))+
labs(title = antigen.specie, color = NULL)+
guides(color = guide_legend(override.aes = list(size= 5))) +
theme(text=element_text(color='black'),
     strip.text=element_text(color='black'), legend.text=element_text(color='black'))
options(repr.plot.width=15.4, repr.plot.height=5)
p
my_saveplot_fun(p, paste0('thymus_', antigen.specie, '_TCR'), width = 15.4, height = 5)


###Fig.7l right
options(repr.plot.width=5, repr.plot.height=3.5)
plot_df <- T_dev_3@meta.data[, c('CMV', 'batch', 'age_3'), 'drop'  = F]
plot_df$barcode <- rownames(plot_df)
plot_df$CMV_count <- ifelse(!plot_df$CMV %in% c('CD8+Memory_like', 'CD4+Memory_like'), 0, 1)

plot_df <- ddply(plot_df, c('batch'), summarise,  CMV_count = sum(CMV_count), age_3 = unique(age_3))


plot_df_sum <- ddply(thymus_tcr_df[!duplicated(thymus_tcr_df$barcode),], 'Sample',summarise, count = length(Sample))
colnames(plot_df_sum)[colnames(plot_df_sum) =='Sample'] <- 'batch'
plot_df <- merge(plot_df, plot_df_sum)
plot_df$Freq <- plot_df$CMV_count/plot_df$count * 100
plot_df$age_3 <- factor(plot_df$age_3, levels =rev(c( '<=12','13-39', '40-99')))

stat_test <- wilcox_test(
 Freq ~ age_3, data = plot_df,p.adjust.method = 'fdr',
)%>% #filter(p.adj.signif!='ns')%>% 
add_xy_position(step.increase = 0.2,fun = 'mean')
stat_test$y.position <- seq(min(stat_test$y.position), max(stat_test$y.position), 
                            length.out = length(stat_test$y.position))+0.5


p2 <- ggplot(plot_df, aes(age_3, Freq))+
geom_boxplot(mapping = aes(fill = age_3))+
geom_point(position = 'jitter')+
guides(fill = F)+
scale_fill_manual(values = age_colors)+
    theme(panel.background = element_rect(fill = 'white', color = NA),plot.margin = margin(0, 0, 0, 0, 'mm'), 
          panel.border = element_rect(fill = NA, color = 'black'), 
          axis.title.x = element_text(size=13.75),axis.text.x = element_text(size=13.75, angle = 45, hjust = 1),
          axis.title.y = element_text(size=13.75, angle = 0, vjust = 0.5),
          axis.ticks = element_blank(),axis.text.y =  element_text(size=13.75),
          plot.title = element_text(size=13.75,hjust=0.5),
          #axis.title.y = element_blank(),
         legend.title = element_text(size = 13.75),legend.text = element_text(size = 13.75), )+
labs(y = 'CMV-specific TCR in effector/memory cells %', x = NULL, fill = 'Age')+
coord_flip()+
scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))+
scale_x_discrete(label = age_group_lbls_func)+
stat_pvalue_manual(data = stat_test, label = 'p.adj.signif', hjust = 0.5,coord.flip= T,
                   size = 13.75*15/13/.pt)

p2
saveRDS(plot_df, file = 'exp_T_thymus_CMV_count_age.rds')
#my_saveplot_fun(p, file = 'exp_T_thymus_CMV_count_age', width = 5, height = 3.5)

options(repr.plot.width=5.1, repr.plot.height=3.5)
plot_list <- get_gg_align(list(p1,p2))
plot_list[[1]]
ggsave(file = 'exp_T_PB_CMV_count_age.svg', width = 5, height = 3.5)
plot_list[[2]]
ggsave(file = 'exp_T_thymus_CMV_count_age.svg', width = 5, height = 3.5)








