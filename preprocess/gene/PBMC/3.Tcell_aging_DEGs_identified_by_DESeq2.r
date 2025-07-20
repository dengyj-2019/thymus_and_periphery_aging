library(Seurat)
library(ggplot2)

library(reshape2)
library(Matrix)
library(clusterProfiler)
library(openxlsx)

library(reshape2)

library(aplot)

library(circlize)

library(patchwork)

library(ComplexHeatmap)

library(BiocParallel)


library(DESeq2, lib.loc = '/home/dengyj/R/x86_64-pc-linux-gnu-library/3.6')

library(glmGamPoi)

library(plyr)

library(ggpubr)

library(future)

module_colors <- c('1' = '#be6b63', '2' = '#86a866', 
                            '3' = '#bda25b', '4' = '#6a82aa', '5' ='#6e578c','6' ='#e5cca9', '7' ='#7cc3cd')

cluster_colors <- c('CD8+TEM-B'='#cb50a0',
                    'CD8+TEM-K'='#db72b4',#
                    'CD8+RTE'='#276d6d','CD4+TN'='#7bacd9',
                        'CD4+CTL'='#7c75ce',#
                    'CD56_Dim_NK'='#c8af1e',#
                    'CD56_Bright_NK'='#d99333',#
                    'CD27lo_Vδ1T'='#93c88e','CD8+TN'='#DDC4DA',
                    'IFN_response'='#00E8F7','CD4+TEM'='#868ae3',#'#909ee7',##'
                    'Treg'='#6e7cad','CD4+TCM'='#80a2e6',#'#9cb9e5',
                    'MAIT'='#b1b772',
                    'CD4+RTE'='#66c5b9',
                        'Vδ2Vγ9T'='#A0D7C9',
                    'NKT'='#e3cd82',#
                    'CD27+Vδ1T'='#9cb39b',#
                    'IFN_T'='#ec807c',
                    'T_un'='#747474','Tex'='#A1B4BB', 'CD4+TEFF' = '#55b3db',##
                    'CD8+TCM'='#E59CC4')

age_colors <- c('<=12'="#b7d9ee",'Prepuberal' = '#b7d9ee','Prepubertal' = '#b7d9ee',
                '13-39'="#3ca88e",'Adult' = '#3ca88e', 
                'Group1' = '#b7d9ee', 'Group2' = '#3ca88e', 'Group3' = '#ec8f46',
                '>=40'="#ec8f46", '40-99'="#ec8f46",'Aged' = '#ec8f46', '>=100' = '#c04e5c', 'Centenarian' = '#c04e5c')

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')




T_NK_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_2.rds')


#Tcell <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell.rds')




T_NK_2$age_group <- factor(T_NK_2$age_3)
levels(T_NK_2$age_group)[levels(T_NK_2$age_group) == '<=12'] <- 'Prepuberal'
levels(T_NK_2$age_group)[levels(T_NK_2$age_group) == '13-39'] <- 'Adult'
levels(T_NK_2$age_group)[levels(T_NK_2$age_group) == '40-99'] <- 'Aged'
levels(T_NK_2$age_group)[levels(T_NK_2$age_group) == '>=100'] <-  'Centenarian'


T_NK_2$age_group <- factor(T_NK_2$age_group, levels = c('Prepuberal', 'Adult', 'Aged', 'Centenarian'))

T_sub <- T_NK_2[,Idents(T_NK_2) %in% c('CD4+RTE',
'CD8+RTE',                
'CD4+TN',
'CD4+TCM',
'CD8+TN',
'CD8+TCM',
'CD4+TEFF',
'Treg',
'CD4+TEM',
'CD8+TEM-K',
'CD4+CTL',
'CD8+TEM-B',
'Tex', 'IFN_T')]


T_sub <- RenameIdents(T_sub, 'CD4+RTE'='CD4+TN', 'CD8+RTE'='CD8+TN')

levels(T_sub) <- c('CD4+RTE',               
                    'CD4+TN',
                    'CD4+TCM',
                    'CD4+TEFF',
                    'CD4+TEM',
                    'CD4+CTL',
                    'Treg',                      
                    'CD8+RTE',                    
                    'CD8+TN',
                    'CD8+TCM',                   
                    'CD8+TEM-K',
                    'CD8+TEM-B',
                    'Tex', 'IFN_T')

group_var <- factor(paste0(T_sub$age_group, '~',Idents(T_sub)))
levels(group_var) <- unlist(lapply(levels(T_sub$age_group), function(x){
    paste0(x, '~', levels(Idents(T_sub)))
}))

T_sub$idents <- Idents(T_sub)

SM <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/SM/humanSM.rds')
TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')

T_sub$age_group <- factor(T_sub$age_3)
levels(T_sub$age_group)[levels(T_sub$age_group) == '<=12'] <- 'Prepuberal'
levels(T_sub$age_group)[levels(T_sub$age_group) == '13-39'] <- 'Adult'
levels(T_sub$age_group)[levels(T_sub$age_group) == '40-99'] <- 'Aged'
levels(T_sub$age_group)[levels(T_sub$age_group) == '>=100'] <-  'Centenarian'

mat <- T_sub$RNA@counts
mat <- mat[rowSums(mat) > 0, ]
mat <- as.matrix(mat)

dds <- DESeqDataSetFromMatrix(countData = mat, 
                              colData = T_sub@meta.data, 
                              design = ~ idents * age_group)

sizeFactors(dds) <-  scran::computeSumFactors(counts(dds))

dds <- DESeq(dds, test = 'LRT',   reduced = ~idents,  fitType = "glmGamPoi", parallel = F,
             #BPPARAM = MulticoreParam(workers = 20), 
             useT = T, minmu=1e-6, minReplicatesForReplace=Inf)

saveRDS(dds, file = 'dds_result.rds')



DESeq2_features <- lapply(c('age_groupPrepuberal','age_groupAdult','age_groupAged'), function(group){
    df <- results(dds, name = group)
    df <- df[!is.na(df$padj),]
    features <- rownames(df)[df$padj < 0.05 & abs(df$log2FoldChange) < log2(1.5)]
})

DESeq2_features <- unique(unlist(DESeq2_features))



DS_res_idents <- c('identsCD4.TCM','identsCD4.TEFF','identsCD4.TEM',
'identsCD4.CTL','identsTreg','identsCD8.TN','identsCD8.TCM',
'identsCD8.TEM.K','identsCD8.TEM.B','identsTex','identsIFN_T')
DS_idents_features <- lapply(DS_res_idents, function(DS_res_ident){
    df <- results(dds, name = DS_res_ident)
    df <- df[!is.na(df$padj),]
    features <- rownames(df)[df$padj < 0.05 & abs(df$log2FoldChange) > log2(2)]    
})

DS_idents_features <- unique(unlist(DS_idents_features))

DS_res_age <- c('age_groupPrepuberal','age_groupAdult','age_groupAged')
DS_age_features <- lapply(DS_res_age, function(DS_res_ident){
    df <- results(dds, name = DS_res_ident)
    df <- df[!is.na(df$padj),]
    features <- rownames(df)[df$padj < 0.05 & abs(df$log2FoldChange) > log2(1.5)]    
})

DS_age_features <- unique(unlist(DS_age_features))



DS_age_specific_features <- setdiff(DS_age_features, DS_idents_features)

length(DS_age_specific_features)

source('fast_ssgsea.R')

go <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/c5.all.v7.0.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(go), function(x) {go[[x]][1]})
go <- sapply(1:length(go), function(x) {go[[x]][-c(1,2)]})
names(go) <- geneset_names

reactome <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/ReactomePathways.gmt"), split = '\t')
geneset_names <- sapply(1:length(reactome), function(x) {reactome[[x]][1]})
reactome <- sapply(1:length(reactome), function(x) {reactome[[x]][-c(1,2)]})
names(reactome) <- geneset_names

kegg <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/c2.cp.kegg.v7.0.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(kegg), function(x) {kegg[[x]][1]})
kegg <- sapply(1:length(kegg), function(x) {kegg[[x]][-c(1,2)]})
names(kegg) <- geneset_names

genesets_db <- c(go, kegg, reactome)



exp <- T_sub$RNA@data[rowSums(T_sub$RNA@data) > 0, ]

exp <- exp[DS_age_specific_features, ]


options(future.globals.maxSize = 64 * 1024^3)
plan('multicore', workers = 16)
T_sub_es <- ssgsea(exp = exp , genesets = genesets_db, 
                      min_size = 10, max_size = 500, exp_filter = F, normalization = F)
plan('sequential')

saveRDS(T_sub_es, file = 'T_sub_es.rds')




















get_Age_GS_test <- function(
    data,
    idents_var, 
    age_var,
    clusters, 
    file = 'Age_GS_test_list.rds', 
    xlsx_file = 'Age_GS_test_list.xlsx',
    database, 
    genes,
    ...
){
    if(file.exists(file)){
        Age_GS_test_list <- readRDS(file)
    }
    else{
        Age_GS_test_list <- list()
    }
    exists_clu <- intersect(names(Age_GS_test_list), clusters)
    clusters <- unique(setdiff(clusters, names(Age_GS_test_list)))
    if(length(exists_clu) > 0){
        print(paste0(paste(exists_clu, collapse = ', '), ' had been calculated, ignored these clusters'))
        print(paste0('Caculated the rest clusters: ', paste(clusters, collapse = ', ')))        
    }

    new_Age_GS_test_list <- lapply(clusters, function(cluster){
        data_cluster <- data[, idents_var %in%  cluster]
        idents_var_cluster <- idents_var[idents_var %in% cluster]
        age_var_cluster <- age_var[idents_var %in% cluster]
        result_cluster <- my_nonparametric_test(data = data_cluster, group_var = age_var_cluster, ...)
        result_cluster$genes <- get_pathway_gene(
            pathways = result_cluster[, grep('variable$', colnames(result_cluster))[1]], 
            database = database,
            genes = genes
        )
        return(result_cluster)
    })
    names(new_Age_GS_test_list) <- clusters
    Age_GS_test_list <- c(Age_GS_test_list, new_Age_GS_test_list)
    
    
    saveRDS(Age_GS_test_list, file = file)
    
    wb <- createWorkbook()
    for(i in names(Age_GS_test_list)){
        addWorksheet(wb, sheetName = i)
        writeData(wb, sheet = i, rowNames = F, x = Age_GS_test_list[[i]])
    }
    saveWorkbook(wb, file = xlsx_file, overwrite = TRUE) 
    
    return(Age_GS_test_list)
}

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

source('/data1/02.private/dengyj/analysis/mycode/Cal_enrichment.R')


T_sub$new_idents <- as.character(Idents(T_sub))
T_sub$new_idents[grepl('CD8.*TN|CD8.*RTE',T_sub$new_idents)] <- 'CD8+TN'
T_sub$new_idents[grepl('CD4.*TN|CD4.*RTE',T_sub$new_idents)] <- 'CD4+TN'
T_sub_identity <- as.character(T_sub$new_idents)
T_sub_age <- T_sub$age_group

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)

options(warn=-1)
Age_GS_test_list <- 
get_Age_GS_test(data = T_sub_es, clusters = unique(T_sub_identity), database = database, 
                multiprocess=T,combined_rest_group = T,
                genes = rownames(T_sub$RNA@data)[rowSums(T_sub$RNA@data) > 0],
                idents_var = T_sub_identity, age_var = T_sub_age)
options(warn=0)


plan('sequential')



Age_GS_list <- lapply(levels(T_sub) , function(idents){
    data <- Age_GS_test_list[[idents]]
    aging_group <- colnames(data)[grep('_comparison', colnames(data))]
    aging_group <- gsub('_VS_Rest_comparison', '', aging_group)
    Age_GS_list_idents <- lapply(aging_group, function(age){
        #data_iter <- data[data$cluster == age, ]
        tmp_effect_size <- paste0(age, '_VS_Rest_effect_size')
        tmp_p_val_adj <- paste0(age, '_VS_Rest_p_val_adj')
        tmp_p_val <- paste0(age, '_VS_Rest_p_val')
        tmp_pct_1 <- paste0(age, '_VS_Rest_pct_1')
        tmp_pct_2 <- paste0(age, '_VS_Rest_pct_2')
        tmp_variable <- paste0(age, '_VS_Rest_variable')
        sig_GS <- data[abs(data[, tmp_effect_size]) > 0.2 & #data[, tmp_p_val_adj] < 0.05 &
                        data[, tmp_p_val_adj] < 0.05, tmp_variable]
#         down_DEGs <- data[data[, tmp_effect_size] < -log2(1.5) & #data[, tmp_p_val_adj] < 0.05 &
#                           data[, tmp_p_val]* nrow(T_NK_2) < 0.05 & ###手动矫正P值
#                         (data[, tmp_pct_1] > 0.1 | data[, tmp_pct_2] > 0.1), tmp_variable]
        return(list(sig_GS=sig_GS))#, down_DEGs=down_DEGs))
    })
    names(Age_GS_list_idents) <- aging_group
    Age_GS_list_idents
})
names(Age_GS_list) <- levels(T_sub) 

Age_GS_list_T <- lapply(names(Age_GS_list[[1]]), function(x){
    lapply(Age_GS_list, function(y){
        y[[x]]
    })
})
names(Age_GS_list_T) <- names(Age_GS_list[[1]])

Age_GS_shared_number <- lapply(names(Age_GS_list_T), function(pair_name){
    pair <- Age_GS_list_T[[pair_name]]
    sig_GS <- unlist(lapply(pair, function(x){
        x$sig_GS
    }))
    sig_GS <- table(table(sig_GS))
#     down_GS <- unlist(lapply(pair, function(x){
#         x$down_GS
#     })) 
#     down_GS <- table(table(down_GS))
#     common <- intersect(names(up_GS), names(down_GS))
#     total <- up_GS[common]+down_GS[common]
#     total <- c(total, up_GS[setdiff(names(up_GS), common)], 
#               down_GS[setdiff(names(down_GS), common)])
    df <- data.frame(sig_GS)
    colnames(df) <- c('count', 'number')
    df$count <- as.numeric(df$count)
    df <- df[order(df$count), ]
    df$pair_name <- pair_name
    df
})
Age_GS_shared_number <- do.call(rbind, Age_GS_shared_number)

min_cutoff <- 5
Age_GS_shared <- lapply(names(Age_GS_list_T), function(pair_name){
    pair <- Age_GS_list_T[[pair_name]]
    sig_GS <- unlist(lapply(pair, function(x){
        x$sig_GS
    }))
    sig_GS <- table(sig_GS)
#     down_GS <- unlist(lapply(pair, function(x){
#         x$down_GS
#     })) 
#     down_GS <- table(down_GS)
    sig_GS_filtered <- names(sig_GS)[sig_GS >= min_cutoff]
    #down_GS_filtered <- names(down_GS)[down_GS >= min_cutoff]
    return(sig_GS_filtered)#, down_GS_filtered))
})
names(Age_GS_shared) <- names(Age_GS_list_T)
Age_GS_selected <- unique(unlist(Age_GS_shared))


length(Age_GS_selected)

options(repr.plot.width=12, repr.plot.height=8)
p <- ggplot(Age_GS_shared_number, aes(count, number))+
geom_bar(stat = 'identity')+
geom_text(mapping = aes(label = number))+
facet_wrap(.~pair_name, ncol =2)
p

GS_exp <- lapply(levels(T_sub), function(clu){
    age_tag <- c('Prepuberal', 'Adult', 'Aged', 'Centenarian')
    result <- lapply(age_tag, function(age){
        gs_used <- c(Age_GS_selected)
        exp1 <- T_sub_es[gs_used, Idents(T_sub) %in% 
                                              clu & T_sub$age_group %in% age]
        exp2 <- T_sub_es[gs_used, Idents(T_sub) %in% 
                                              clu & !T_sub$age_group %in% age]        
        val <- sapply(1:nrow(exp1), function(x){
           rcompanion::cliffDelta(x=exp1[x, ], y=exp2[x,]) 
        })
        names(val) <- gs_used
        val
    })
    result <- do.call(cbind, result)
    colnames(result) <- paste0(age_tag, '~', clu)
    result
})

GS_exp <- do.call(cbind, GS_exp)

GS_exp <- GS_exp[, levels(droplevels(group_var))]
saveRDS(GS_exp, file = 'GS_exp.rds')

GS_exp_raw <- GS_exp <- readRDS('GS_exp.rds')

GS_exp <- my_scale(GS_exp, 'row')

#GS_exp <- my_scale(GS_exp, 'row')


options(repr.plot.width=12, repr.plot.height=25)
hr <- hclust(dist(GS_exp, method = "euclidean"), method="complete") 
GS_exp_rl = cutree(hr, 12)
GS_exp_cl <- factor(sapply(strsplit(colnames(GS_exp), split = '_'), function(x) x[1]))
GS_exp_cl <- factor(GS_exp_cl, levels = unique(GS_exp_cl))                
                    
p <- 
Heatmap(GS_exp, cluster_columns = F, show_row_names = T,use_raster = F,column_split=GS_exp_cl,
        heatmap_legend_param=list(title='Expression'),
        row_split = GS_exp_rl,row_dend_reorder=T,
        row_names_gp = gpar(fontsize = 4, fontfamily = 'sans', fontface = 'italic'), 
       top_annotation = HeatmapAnnotation(
               df=data.frame(Age = sapply(strsplit(colnames(GS_exp), split = '~'), function(x) x[1]), 
                             Identity = sapply(strsplit(colnames(GS_exp), split = '~'), function(x) x[2])), 
                 col = list(Age = age_colors, Identity = cluster_colors),
               annotation_name_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = 'bold'),
               annotation_legend_param = list(
                 Age = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(GS_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ), 
               Identity = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(GS_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ))
             ))
                        
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(
    draw(p, merge_legends = T, heatmap_legend_side='right', padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p
#ggsave('PB_aging_DEGs.png', width = 12, height = 25, plot = gg_p, dpi=300) 


##combined similar module
orig_GS_exp_rl <- GS_exp_rl
GS_exp_rl[orig_GS_exp_rl %in% c(1,2)] <- 1
GS_exp_rl[orig_GS_exp_rl %in% c(5)] <- 2
GS_exp_rl[orig_GS_exp_rl %in% c(9)] <- 3
GS_exp_rl[orig_GS_exp_rl %in% c(8)] <- 3
GS_exp_rl[orig_GS_exp_rl %in% c(3,4,6)] <- 4
GS_exp_rl[orig_GS_exp_rl %in% c(7,10,11,12)] <- 5






#GS_exp <- my_scale(GS_exp, 'row')


options(repr.plot.width=12, repr.plot.height=25)
hr <- hclust(dist(GS_exp, method = "euclidean"), method="complete") 
# GS_exp_rl = cutree(hr, 8)
GS_exp_cl <- factor(sapply(strsplit(colnames(GS_exp), split = '_'), function(x) x[1]))
GS_exp_cl <- factor(GS_exp_cl, levels = unique(GS_exp_cl))                
                    
GS_exp_p <- 
Heatmap(GS_exp, cluster_columns = F, show_row_names = T,use_raster = F,column_split=GS_exp_cl,
        heatmap_legend_param=list(title='Expression'),
        row_split = GS_exp_rl,row_dend_reorder=T,
        row_names_gp = gpar(fontsize = 4, fontfamily = 'sans', fontface = 'italic'), 
       top_annotation = HeatmapAnnotation(
               df=data.frame(Age = sapply(strsplit(colnames(GS_exp), split = '~'), function(x) x[1]), 
                             Identity = sapply(strsplit(colnames(GS_exp), split = '~'), function(x) x[2])), 
                 col = list(Age = age_colors, Identity = cluster_colors),
               annotation_name_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = 'bold'),
               annotation_legend_param = list(
                 Age = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(GS_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ), 
               Identity = list(
                   legend_direction = "vertical",
                   #title = "Baaaaaaar",
                   #at = colnames(GS_exp),
                   title_gp = gpar(fontsize = 10, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 10, fontfamily = 'sans'),
                   ncol = 1
                 ))
             ))
                        
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(
    draw(GS_exp_p, merge_legends = T, heatmap_legend_side='right', padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p
#ggsave('PB_aging_DEGs.png', width = 12, height = 25, plot = gg_p, dpi=300) 




GS_exp_var <- apply(GS_exp_raw, 1, var)

plot_features <- c('Chromatin organization',
'GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
'GO_ELECTRON_TRANSPORT_CHAIN',
'GO_HISTONE_BINDING',
'GO_MITOCHONDRIAL_MEMBRANE_PART',
'GO_MITOCHONDRION_ORGANIZATION',
'GO_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS',
'GO_TYPE_I_INTERFERON_PRODUCTION',
'Innate Immune System',
'KEGG_OXIDATIVE_PHOSPHORYLATION',
'rRNA processing in the mitochondrion',
'tRNA processing in the mitochondrion',
'GO_LEUKOCYTE_CELL_CELL_ADHESION',
'GO_PROTEIN_PROCESSING',
'GO_REGULATION_OF_HISTONE_MODIFICATION',
'GO_REGULATION_OF_PROTEIN_MODIFICATION_PROCESS',
'GO_REGULATION_OF_PROTEOLYSIS',
'KEGG_CHEMOKINE_SIGNALING_PATHWAY',
'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
'KEGG_WNT_SIGNALING_PATHWAY',
'PIP3 activates AKT signaling',
'Signaling by TGF-beta family members',
'Transcriptional regulation by RUNX1',
'GO_OXIDOREDUCTASE_ACTIVITY')


####这里是按照module区分，但是只有2个module对应着Prepuberal和Adult组，所以基本上和按照年龄区分的结果差不多
GS_exp_rl_filtered <- lapply(unique(GS_exp_rl), function(iter){
    GS_exp_rl_iter <- GS_exp_rl[GS_exp_rl==iter]
    GS_exp_var_filtered <- GS_exp_var[names(GS_exp_rl_iter)]
    GS_exp_var_filtered <- GS_exp_var_filtered[order(GS_exp_var_filtered, decreasing = T)]
    result <- names(GS_exp_var_filtered)[1:100]
    result <- union(result, intersect(plot_features, 
                                      names(GS_exp_var_filtered)))
    result <- result[!is.na(result)]
    result
})

features_lvls <- intersect(rownames(GS_exp)[unlist(row_order(GS_exp_p))], unlist(GS_exp_rl_filtered))



new_GS_mean <- my_group_mean(T_sub_es[features_lvls, ], factor(T_sub$age_group))

new_GS_mean <- my_scale(new_GS_mean, do_scale = 'row')

new_GS_mean <- new_GS_mean[, c('Prepuberal','Adult','Aged','Centenarian')]

####Extended Data Figure 6e
new_GS_mean <- new_GS_mean[features_lvls, ]
plot_features <- intersect(rownames(new_GS_mean), plot_features)
text_colors <- module_colors[as.character(GS_exp_rl[intersect(features_lvls, plot_features)])]
row_split <- factor(GS_exp_rl[features_lvls], 
                    levels = as.character(sort(unique(GS_exp_rl[features_lvls]))))
p <- Heatmap(new_GS_mean,col = colorRamp2(c(-3,-1.5,0,1.5,3), 
                                                c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b')),

             #row_title = NULL,
             cluster_columns = F, cluster_rows = F, show_row_names = F,show_column_names = F,show_heatmap_legend=F,
         #row_split = GS_exp_cl, 
        row_split = row_split,
        heatmap_legend_param=list(title='Cliff’s delta', direction = "vertical", 
                                  labels_gp = gpar(fontsize = 13.75),title_gp = gpar(fontsize = 13.75)),
        row_names_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = 'italic'), 
        column_names_centered=T,
       top_annotation = HeatmapAnnotation(show_legend = F,show_annotation_name = T,
                                       
               df=data.frame(Age = colnames(new_GS_mean)
                             ), 
                 col = list(Age = age_colors[colnames(new_GS_mean)]),
                 annotation_name_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = 'bold'),                              
               annotation_legend_param = list(
                 Age = list(
                   title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
                   nrow = 4))
             ),               
            right_annotation = rowAnnotation(' ' = anno_mark(link_width = unit(2.5,'mm'),side = 'right',
                                                                          #at_adj = -0.5,
                                                                          #extra_text_args = list(vjust = 0.5),
            at = which(rownames(new_GS_mean) %in% plot_features), 
            labels = str_to_sentence_func(gsub('_', ' ',gsub('GO_|KEGG_', '', 
                                                                   plot_features))), 
            labels_gp = gpar(fontsize=13.75,col = text_colors,#'black',#text_colors, 
                             fontface='plain')), 
           'foo' = anno_empty(border = FALSE, 
                              width=
                              (max_text_width(plot_features) - unit(12,'cm')) )))
options(repr.plot.width=8.7, repr.plot.height=9.5)
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p, heatmap_legend_side="bottom",  merge_legends = T,
                                                     # annotation_legend_side="right",
                                                      legend_grouping = "original",#ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p

####
ggsave('/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig5_FigS5/PB_Tcell_total_age_DEGs.svg', 
       width = 8.7, height = 9.5) 

####Extended Data Figure 6e
Identity <- unique(sapply(strsplit(colnames(new_GS_mean), split = '~'), function(x) x[2]))
     
Age <- age_group_lbls_func(unique(sapply(strsplit(colnames(new_GS_mean), split = '~'), function(x) x[1])))
lgd2 = Legend(labels =  Age, title = "Age", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 2,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(age_colors[Age])))

                     
Module <- sort(unique(GS_exp_rl[features_lvls]))
lgd3 = Legend(labels =  Module, title = "Module", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 2,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(module_colors[as.character(Module)])))     

col_fun =  colorRamp2(c(-3,-1.5,0,1.5,3), c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))
lgd4 = Legend(col_fun = col_fun, title = "Cliff’s delta", direction = 'horizontal',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
             )                                          

pd = packLegend(lgd2,lgd3,
                lgd4, direction = "horizontal")

options(repr.plot.width=5, repr.plot.height=0.8)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(0, "cm"), 
                                              y = unit(0.2, "cm"), just = c("left", "bottom"))))
ggsave('/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig5_FigS5/PB_Tcell_total_age_DEGs_legend.svg', 
       width = 12, height = 1) 


