library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(treemapify)
library(circlize)
library(dplyr)
library(plyr)
library(vegan)
library(Rmisc)
library(rcompanion)
library(ineq)
library(ggsignif)
library(ggpubr)
setwd('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/')
library(ComplexHeatmap)
source('/data1/02.private/dengyj/analysis/thymus/TCR/TCR_caculation.R')
library(patchwork)
library(rlang)
library(ggalluvial)

library(reshape2)

library(rstatix)

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


ETP_DN_identity <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_DN_tcr_idents.rds')

T_dev_identity <- 
as.character(Idents(T_dev_3))
names(T_dev_identity) <- names(Idents(T_dev_3))

thymus_identity <- c(T_dev_identity, ETP_DN_identity)

thymus_identity[grep('DN', thymus_identity)] <- 'DN'

lbls_func <- function(x){
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
}

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

age_colors <- c('<=12'="#b7d9ee",'Prepuberal' = '#b7d9ee', '13-39'="#3ca88e",'Adult' = '#3ca88e', 
                '>=40'="#ec8f46",'Aged' = '#ec8f46')

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

T_idents <- c('DN','DP_P', 'DP_Q',
              'HSP_DP',
              'αβ_Entry','Agonist_1','CD4+pre_T', 'CD8+pre_T', 'Agonist_2', 
                    'CD4+Naive', 'CD8+Naive', 'IFN_response','pre_Treg1', 'pre_Treg2','Treg', 'CD4+Memory_like',
                    'CD8+Memory_like')
other_idents <- setdiff(unique(thymus_identity), T_idents)

clusters_used <- T_idents

setwd('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/data/')
files <- dir()
files <- files[-grep('metadata.txt', files)]
tcr_raw <- list()
tcr_raw$meta <- read.table('metadata.txt', header = T)  
tcr_raw$data <- 
lapply(files, function(i){
    df_iter <- read.csv(i)
    df_iter$Sample <- strsplit(i, split = '\\.')[[1]][1]
    df_iter$age <- tcr_raw$meta[match(df_iter$Sample, tcr_raw$meta$Sample), 'age']
    df_iter$age_3 <- tcr_raw$meta[match(df_iter$Sample, tcr_raw$meta$Sample), 'age_3']
    df_iter$identity <- thymus_identity[match(df_iter$barcode, names(thymus_identity))]
    df_iter <- df_iter[df_iter$identity %in% T_idents, ]
    df_iter$type <- 'other'


    df_iter <- df_iter[df_iter$productive == 'True', ]
    tbl=table(df_iter$barcode, df_iter$chain)
    logi <- apply(tbl, 1, function(x){
        x['TRA'] > 0 & x['TRB'] > 0 
    })  
    logi2 <- apply(tbl, 1, function(x){
        x['TRB'] > 0 
    }) & rownames(tbl) %in% df_iter$barcode[grepl('DN',df_iter$identity)]
    df_iter <- df_iter[df_iter$barcode %in% rownames(tbl)[logi | logi2], ]
    return(df_iter)
})
names(tcr_raw$data) <- sapply(strsplit(files, split = '\\.'), function(x) x[1])                          
setwd('/data1/02.private/dengyj/analysis/thymus/TCR/thymus/')                          

tcr_raw$meta$age_3 <- factor(tcr_raw$meta$age_3, levels =  c('<=12', '13-39', '>=40'))

cluster_colors <- c('pre_Treg1'='#E5D2DD','CD8+Naive'='#53A85F','CD4+Memory_like'='#808040',
                    'CD8+Memory_like'='#F3B1A0',
                             'CD4+Naive'='#D6E7A3','CD4+GZMK+T'='#57C3F3','IFN_response'='#476D87','Treg'='#E95C59',
                    'CD4+pre_T'='#83c3bf',
                            'DP_P'='#AB3282','DP_Q'='#A8B0E3','Agonist_2'='#a83117','Exhausted_T'='#8C549C',
                    'αβ_Entry'='#c28578',
                             'HSP_DP'='#9FA3A8','IFN_response_2'='#CCE0F5','CD8+pre_T'='#625D9E', 
                    'Agonist_1' = '#a6c050', 'HSP_SP' = '#d98ed6', #'αβ_Entry_2b' = '#cb89d1', 
                    'pre_Treg2' = '#9566ac' , 'None' = 'lightgrey', 
                    'ETP_1' = '#297A6B','ETP_2'='#C05EA5','ETP'='#a7aad5',
    'ETP_3'='#e6ba77','DN_P'='#D8CEEE','DN_Q1' = '#4279C0','DN_Q2'='#42BACA','DN'='#42BACA',
    'DN4'='#ff927c','γδT_P'='#ee587d'
                            )

tcr_df <- do.call(rbind, tcr_raw$data)    
tcr_df <- tcr_df[tcr_df$productive == 'True',]
rownames(tcr_df) <- NULL


tcr_df$identity_2 <- tcr_df$identity
# tcr_df$identity_2[tcr_df$identity %in% c('DN_Q1', 'DN_Q2', 'DN_P', 'DN4')] <- 'DN'
unique_tcr_df <- tcr_df[!duplicated(tcr_df$barcode), ]
tcr_cells <- table(unique_tcr_df$identity_2)
tcr_cells[order(tcr_cells, decreasing = T)]
thymus_identity_2 <- thymus_identity
#thymus_identity_2[thymus_identity %in% c('DN_Q1', 'DN_Q2', 'DN_P', 'DN4')] <- 'DN'
tcr_proprotion <- tcr_cells/table(thymus_identity_2)[names(tcr_cells)] * 100
tcr_proprotion[order(tcr_proprotion, decreasing = T)]


tra_df <- tcr_df[tcr_df$productive == 'True' & tcr_df$chain == 'TRA' & !grepl('DN', tcr_df$identity), ]
trb_df <- tcr_df[tcr_df$productive == 'True' & tcr_df$chain == 'TRB', ]

classic_T_idents <- c('DN','DP_P', 'DP_Q',  
              'αβ_Entry','CD4+pre_T', 'CD8+pre_T', 
                    'CD4+Naive', 'CD8+Naive')

trav_clu_df <- get_usage_cluster_df(tra_df, 'v_gene', clusters = setdiff(classic_T_idents, 'DN'))
traj_clu_df <- get_usage_cluster_df(tra_df, 'j_gene', clusters = setdiff(classic_T_idents, 'DN'))
trbv_clu_df <- get_usage_cluster_df(trb_df, 'v_gene', clusters = classic_T_idents)
trbj_clu_df <- get_usage_cluster_df(trb_df, 'j_gene', clusters = classic_T_idents)

tcr_chain_name <- c('trav', 'traj', 'trbv', 'trbj')

tcr_clu_test <- lapply(tcr_chain_name, function(x){
    obj <- get(paste0(x, '_clu_df'))
    result <- single_geneuseage_test(obj[, !colnames(obj) %in% c('Sample', 'age_3')], 
                                 group_var = 'cluster')
    return(result)
})
names(tcr_clu_test) <- paste0(tcr_chain_name, '_clu_test')



options(warn=-1)
p_val_cutoff <- 0.05
max_cutoff <- 2
min_cutoff <- -2
tcr_clu_plot_df <- lapply(tcr_clu_test, function(x){
    df_iter <- x[, grep('p_val_adj', colnames(x))]
    x <- x[apply(df_iter, 1, function(y) {any(y[!is.na(y)] < p_val_cutoff)}), ]
    if(nrow(x) > 0){
        x[, grep('useage', colnames(x))] <- my_scale(x[, grep('useage', colnames(x))], 'row')
        x[, grep('useage', colnames(x))][x[, grep('useage', colnames(x))] > max_cutoff] <- max_cutoff
        x[, grep('useage', colnames(x))][x[, grep('useage', colnames(x))] < min_cutoff] <- min_cutoff
        #x <- x[, grep('p_val_adj|useage', colnames(x))]
        x$gene <- rownames(x)
        x$type <- substring(x$gene,1,4)
        useage_df <- melt(x[, grep('useage|gene|type', colnames(x))], id.vars=c('gene', 'type'), 
                         variable.name = 'cluster', value.name = 'useage')
        useage_df$cluster <- gsub('_useage', '',useage_df$cluster)
        p_val_adj_df <- melt(x[, grep('p_val_adj|gene|type', colnames(x))], id.vars=c('gene', 'type'), 
                         variable.name = 'cluster', value.name = 'p_val_adj')
        p_val_adj_df$cluster <- gsub('_p_val_adj', '',p_val_adj_df$cluster) 
        df <- merge(useage_df, p_val_adj_df, by = c('gene', 'type', 'cluster'))
        return(df)        
    }
})
tcr_clu_plot_df <- do.call(rbind, tcr_clu_plot_df)
tcr_clu_plot_df$type <- factor(tcr_clu_plot_df$type, levels = c('TRAV', 'TRAJ', 'TRBV', 'TRBJ'))

tcr_clu_plot_df_gene_order <- lapply(tcr_clu_test, function(x){
    df_iter <- x[, grep('p_val_adj', colnames(x))]
    x <- x[apply(df_iter, 1, function(y) {any(y[!is.na(y)] < p_val_cutoff)}), ]
    mat <- x[, grep('useage', colnames(x))]
    colnames(mat) <- gsub('_useage', '',colnames(mat))
    mat <- order_heatmap(mat);return(rownames(mat))
    #clust <- hclust(dist(mat));return(rownames(mat)[clust$order])
})

tcr_clu_plot_df$gene <- factor(tcr_clu_plot_df$gene, 
                              levels = unlist(tcr_clu_plot_df_gene_order, use.names = F))
tcr_clu_plot_df$cluster <- factor(tcr_clu_plot_df$cluster, levels = classic_T_idents)

options(warn=0)

max_val <- max(-log10(tcr_clu_plot_df$p_val_adj))
plot_list=list()
#label_list <- list()
for(i in unique(tcr_clu_plot_df$type)){
    p <- ggplot(tcr_clu_plot_df[tcr_clu_plot_df$type == i, ], 
                aes(cluster, gene, fill = useage ,size = -log10(p_val_adj)))+
    geom_point(shape = 21, color = 'black')+
    scale_fill_gradientn(colours = viridis::viridis(20), name = 'Gene Useage', limits=c(min_cutoff, max_cutoff))+
    scale_x_discrete(label = lbls_func)+
    scale_size_continuous(range = c(2,10), name = '-log10(FDR)', 
                          limits = c(0,max_val),
                          breaks = c(0:floor(max_val)))+
    theme(panel.background = element_rect(fill = 'white', color = NA),
          panel.border = element_rect(fill = NA, color = 'black'),
          #plot.background=element_rect(fill = NA, color = 'black'),
          axis.text.x = element_text(size = 15, angle=45, hjust=1),
          axis.text.y = element_text(size = 15,  hjust=1),
          plot.title = element_textbox(size=20,hjust=0.5, vjust=0),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=20), 
         legend.title = element_text(size = 15),legend.text = element_text(size = 15))+
    guides(size=guide_legend(override.aes=list(fill='black', alpha=1)))
    p <- p+labs(title = i)
    plot_list[[i]] <- p
}


options(repr.plot.width=12, repr.plot.height=15)
tcr_clu_plot <- wrap_plots(plot_list)+plot_layout(ncol = 2, nrow=2, guides ='collect')
tcr_clu_plot



options(warn=-1)
trav_df <- as.data.frame(t(apply(table(tra_df$v_gene, tra_df$Sample), 2, FUN = function(x) x/sum(x))))
traj_df <- as.data.frame(t(apply(table(tra_df$j_gene, tra_df$Sample), 2, FUN = function(x) x/sum(x))))
trbv_df <- as.data.frame(t(apply(table(trb_df$v_gene, trb_df$Sample), 2, FUN = function(x) x/sum(x))))               
trbj_df <- as.data.frame(t(apply(table(trb_df$j_gene, trb_df$Sample), 2, FUN = function(x) x/sum(x))))
trav_df <- cbind(trav_df,'age_3' = tcr_raw$meta[match(rownames(trav_df), tcr_raw$meta$Sample), 'age_3'])  
traj_df <- cbind(traj_df,'age_3' = tcr_raw$meta[match(rownames(traj_df), tcr_raw$meta$Sample), 'age_3']) 
trbv_df <- cbind(trbv_df,'age_3' = tcr_raw$meta[match(rownames(trbv_df), tcr_raw$meta$Sample), 'age_3']) 
trbj_df <- cbind(trbj_df,'age_3' = tcr_raw$meta[match(rownames(trbj_df), tcr_raw$meta$Sample), 'age_3'])
                                 
tcr_chain_name <- c('trav', 'traj', 'trbv', 'trbj')

tcr_age_test <- lapply(tcr_chain_name, function(x){
    obj <- get(paste0(x, '_df'))
    result <- single_geneuseage_test(obj, group_var = 'age_3')
    return(result)
})
names(tcr_age_test) <- paste0(tcr_chain_name, '_age_test')
options(warn=0)                                 

p_val_cutoff <- 0.05
max_cutoff <- 1
min_cutoff <- -1
options(warn=-1)
tcr_age_plot_df <- lapply(tcr_age_test, function(x){
    df_iter <- x[, grep('p_val_adj', colnames(x))]
    x <- x[apply(df_iter, 1, function(y) {any(y[!is.na(y)] < p_val_cutoff)}), ]
    if(nrow(x) > 0){
        x[, grep('useage', colnames(x))] <- my_scale(x[, grep('useage', colnames(x))], 'row')
        x[, grep('useage', colnames(x))][x[, grep('useage', colnames(x))] > max_cutoff] <- max_cutoff
        x[, grep('useage', colnames(x))][x[, grep('useage', colnames(x))] < min_cutoff] <- min_cutoff
        #x <- x[, grep('p_val_adj|useage', colnames(x))]
        x$gene <- rownames(x)
        x$type <- substring(x$gene,1,4)
        useage_df <- melt(x[, grep('useage|gene|type', colnames(x))], id.vars=c('gene', 'type'), 
                         variable.name = 'age', value.name = 'useage')
        useage_df$age <- gsub('_useage', '',useage_df$age)
        p_val_adj_df <- melt(x[, grep('p_val_adj|gene|type', colnames(x))], id.vars=c('gene', 'type'), 
                         variable.name = 'age', value.name = 'p_val_adj')
        p_val_adj_df$age <- gsub('_p_val_adj', '',p_val_adj_df$age) 
        df <- merge(useage_df, p_val_adj_df, by = c('gene', 'type', 'age'))
        return(df)        
    }
})
tcr_age_plot_df <- do.call(rbind, tcr_age_plot_df)
tcr_age_plot_df$type <- factor(tcr_age_plot_df$type, levels = c('TRAV', 'TRAJ', 'TRBV', 'TRBJ'))


tcr_age_plot_df_gene_order <- lapply(levels(droplevels(tcr_age_plot_df$type)), function(i){
    x <- tcr_age_plot_df[tcr_age_plot_df$type == i, ,drop = F]
    mat <- data.frame(dcast(x, gene~age, value.var = 'useage'), check.names = F)
    rownames(mat) <- mat$gene
    mat <- mat[, c('<=12','13-39','>=40')]
    mat <- order_heatmap(mat);return(rownames(mat))
    #clust <- hclust(dist(mat));return(rownames(mat)[clust$order])
})


tcr_age_plot_df$gene <- factor(tcr_age_plot_df$gene, 
                              levels = unlist(tcr_age_plot_df_gene_order, use.names = F))

max_val <- max(-log10(tcr_age_plot_df$p_val_adj))
plot_list=list()
#label_list <- list()
for(i in unique(tcr_age_plot_df$type)){
    p <- ggplot(tcr_age_plot_df[tcr_age_plot_df$type == i, ], 
                aes(age, gene, fill = useage ,size = -log10(p_val_adj)))+
    geom_point(shape = 21, color = 'black')+
    scale_fill_gradientn(colours = viridis::viridis(20), name = 'Gene Useage', limits=c(min_cutoff, max_cutoff))+
    scale_x_discrete(limits = c('<=12', '13-39', '>=40'), 
                    labels = c('Prepuberal', 'Adult', 'Aged'))+
    scale_size_continuous(range = c(2,10), name = '-log10(FDR)', 
                          limits = c(0,max_val),
                          breaks = c(0:floor(max_val)))+
    theme(panel.background = element_rect(fill = 'white', color = NA),
          panel.border = element_rect(fill = NA, color = 'black'),
          #plot.background=element_rect(fill = NA, color = 'black'),
          axis.text.x = element_text(size = 15, angle=45, hjust=1),
          axis.text.y = element_text(size = 15,  hjust=1),
          plot.title = element_textbox(size=20,hjust=0.5, vjust=0),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=20), 
         legend.title = element_text(size = 15),legend.text = element_text(size = 15))+
    guides(size=guide_legend(override.aes=list(fill='black', alpha=1)))
    p <- p+labs(title = i)
    plot_list[[i]] <- p
}
options(warn=0)

options(repr.plot.width=10, repr.plot.height=8)
tcr_age_plot <- wrap_plots(plot_list)+plot_layout(ncol = 2, nrow=2, guides ='collect')
tcr_age_plot

options(warn=-1)
tcr_chain_name <- c('trav', 'traj', 'trbv', 'trbj')
tcr_clu_age_test <- lapply(tcr_chain_name, function(x){
    obj <- get(paste0(x, '_clu_df'))
    if(grepl('tra', x)){
        clusters_used <- setdiff(classic_T_idents, 'DN')
    }else{
        clusters_used <- classic_T_idents
    }
    total_result <- lapply(1:length(clusters_used), function(i){
        result <- 
        single_geneuseage_test(obj[obj$cluster == clusters_used[i], !colnames(obj) %in% c('cluster', 'Sample')], 
                                 group_var = 'age_3')   
        #result$gene <- rownames(result)
        result$cluster <- clusters_used[i]
        #rownames(result) <- NULL
        return(result)        
    })
        
})
names(tcr_clu_age_test) <- paste0(tcr_chain_name, '_clu_age_test')
options(warn=0)

options(warn=-1)
p_val_cutoff <- 0.1
max_cutoff <- 1
min_cutoff <- -1
all_col_names <- c('>=40_p_val','>=40_p_val_adj','>=40_useage',
'13-39_p_val','13-39_p_val_adj','13-39_useage','<=12_p_val','<=12_p_val_adj','<=12_useage')
tcr_clu_age_plot_df <- lapply(tcr_clu_age_test, function(x){
    li_iter <- lapply(x, function(y){
        if(!all(all_col_names %in% colnames(y))){
            lack_col_names <- setdiff(all_col_names, colnames(y))       
            y[, lack_col_names[grep('_p_val', lack_col_names)]] <- 1
            y[, lack_col_names[grep('_useage', lack_col_names)]] <- NA
        }
        y <- y[, c(all_col_names, 'cluster')]
        df_iter <- y[, grep('p_val_adj', colnames(y))]
        y <- y[apply(df_iter, 1, function(y_) {any(y_[!is.na(y_)] < p_val_cutoff)}), ]
        if(nrow(y) > 0){
            y[, grep('useage', colnames(y))] <- my_scale(y[, grep('useage', colnames(y))], 'row')
            y[, grep('useage', colnames(y))][y[, grep('useage', colnames(y))] > max_cutoff] <- max_cutoff
            y[, grep('useage', colnames(y))][y[, grep('useage', colnames(y))] < min_cutoff] <- min_cutoff 
            y$gene <- rownames(y)
            y$type <- substring(y$gene,1,4)
            useage_df <- melt(y[, grep('useage|gene|type|cluster', colnames(y))], 
                              id.vars=c('gene', 'type', 'cluster'), 
                             variable.name = 'age', value.name = 'useage')
            useage_df$age <- gsub('_useage', '',useage_df$age)
            p_val_adj_df <- melt(y[, grep('p_val_adj|gene|type|cluster', colnames(y))], 
                                 id.vars=c('gene', 'type', 'cluster'), 
                             variable.name = 'age', value.name = 'p_val_adj')
            p_val_adj_df$age <- gsub('_p_val_adj', '',p_val_adj_df$age) 
            df <- merge(useage_df, p_val_adj_df, by = c('gene', 'type', 'cluster', 'age'))
            return(df)               
        }
    })
    li_iter <- do.call(rbind, li_iter)
})
tcr_clu_age_plot_df <- do.call(rbind, tcr_clu_age_plot_df)
tcr_clu_age_plot_df$type <- factor(tcr_clu_age_plot_df$type, levels = c('TRAV', 'TRAJ', 'TRBV', 'TRBJ'))
options(warn=0)


tcr_clu_age_plot_df$age <- factor(tcr_clu_age_plot_df$age, levels = c('<=12','13-39','>=40'))
tcr_clu_age_plot_df$cluster <- factor(tcr_clu_age_plot_df$cluster, 
                                      levels = clusters_used[clusters_used %in% unique(tcr_clu_age_plot_df$cluster)])

tcr_clu_age_plot_df$cluster_gene <- paste(tcr_clu_age_plot_df$cluster, tcr_clu_age_plot_df$gene, sep='_')
options(warn=-1)
tcr_clu_age_plot_df_gene_order <- lapply(levels(tcr_clu_age_plot_df$type), function(i){
    x <- tcr_clu_age_plot_df[tcr_clu_age_plot_df$type == i, ]
    x$cluster_gene <- paste0(x$cluster, '_', x$gene)
    mat <- data.frame(dcast(x, cluster_gene~age, value.var = 'useage'), row.names = 'cluster_gene', 
                     check.names = F)
    mat <- mat[, c('<=12','13-39', '>=40')]
    mat <- order_heatmap(mat);return(rownames(mat))
    #clust <- hclust(dist(mat));return(rownames(mat)[clust$order])
})
options(warn=0)

tcr_clu_age_plot_df$cluster_gene <- factor(tcr_clu_age_plot_df$cluster_gene, 
                              levels = unlist(tcr_clu_age_plot_df_gene_order, use.names = F))


tcr_clu_age_plot_df$type <- factor(tcr_clu_age_plot_df$type, 
                                   levels = c('TRAV', 'TRAJ', 'TRBV', 'TRBJ'))

max_val <- max(-log10(tcr_clu_age_plot_df$p_val_adj))

options(repr.plot.width=15, repr.plot.height=15)
plot_list <- list()
for(i in levels(tcr_clu_age_plot_df$type)){

    color_used <- cluster_colors[gsub('_TR.*$', '', 
                        levels(droplevels(tcr_clu_age_plot_df$cluster_gene[tcr_clu_age_plot_df$type == i])))]
    p <- ggplot(tcr_clu_age_plot_df[tcr_clu_age_plot_df$type == i, ], 
                aes(age, cluster_gene, fill = useage ,size = -log10(p_val_adj)))+
    geom_point(shape = 21, color = 'black')+
    scale_fill_gradientn(colours = viridis::viridis(20), name = 'Gene Useage', limits = c(min_cutoff, max_cutoff))+
    scale_x_discrete(labels = c('Prepuberal', 'Adult', 'Aged'))+
    scale_y_discrete(labels = function(y){gsub('^.*_', '', y)})+
    scale_size_continuous(range = c(2,10), name = '-log10(FDR)', 
                          limits = c(0,max_val),
                          breaks = c(0:floor(max_val)))+
    theme(panel.background = element_rect(fill = 'white', color = NA), 
          panel.border = element_rect(fill = NA, color = 'black'),
        #panel.background = element_blank(), 
        axis.line = element_line(color='black'),
          axis.text.x = element_text(size = 15, angle=45, hjust=1),
          #axis.text.y = element_custom(size = 15,   fill = 1:2, alpha = 0.1, hjust = 0),
          axis.text.y = element_text(size = 15, hjust = 1, colour = color_used, 
                                    face = 'bold'),
          plot.title = element_textbox(size=20,hjust=0.5, vjust=0),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=20),
        panel.spacing.y = unit(0,"line"),
          legend.key=element_blank(),
          #legend.box.margin = unit(rep(1.5, 4), 'lines'), 
         legend.title = element_text(size = 15),legend.text = element_text(size = 15))+
    labs(title = i, color = 'Cluster')+
    geom_text(aes(color =cluster, label = cluster_gene, y=cluster_gene), alpha=0,x=1,size=1)+
    scale_color_manual(values = cluster_colors, 
                       limits = clusters_used[clusters_used %in% tcr_clu_age_plot_df$cluster], 
                      breaks = clusters_used[clusters_used %in% tcr_clu_age_plot_df$cluster])+
    guides(color = guide_legend(override.aes=list(labels= clusters_used[clusters_used %in% 
                                                                        tcr_clu_age_plot_df$cluster],
        size=4.5,alpha=1, fontface = 'bold'),label.position='left', label.theme=element_blank(), 
                               order = 1),
           size=guide_legend(override.aes=list(fill='black', alpha=1)))
    #p <- p+facet_grid(cluster_gene~., scale='free_y')
    plot_list[[i]] <- p
}


options(repr.plot.width=7, repr.plot.height=7)
wrap_plots(plot_list)+plot_layout(ncol = 2, nrow=2, guides ='collect')&
theme(legend.box.margin = unit(c(0,0,0,1.5), 'cm'))











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

saveRDS(tcr_clu_age_test, file = 'tcr_clu_age_test.rds')














