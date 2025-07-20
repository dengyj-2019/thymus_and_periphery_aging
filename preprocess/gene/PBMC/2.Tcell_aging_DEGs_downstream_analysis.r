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

library(future)

go <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c5.all.v7.0.symbols.gmt')
kegg <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c2.cp.kegg.v7.0.symbols.gmt')
reactome <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/ReactomePathways.gmt')

database <- rbind(go, kegg, reactome)

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


T_NK_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/T_NK_2.rds')




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

group_var <- factor(paste0(T_sub$age_group, '~',Idents(T_sub)),
levels= unlist(lapply(levels(T_sub$age_group), function(x){
    paste0(x, '~', levels(Idents(T_sub)))
})))

Age_Gene_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/PB/Age_Gene_test_list.rds')

names(Age_Gene_test_list)[names(Age_Gene_test_list) %in% c('CD8+TEM')]  <- 'CD8+TEM-K'
names(Age_Gene_test_list)[names(Age_Gene_test_list) %in% c('CD8+TEFF')]  <- 'CD8+TEM-B'

Age_DEGs_list <- lapply(levels(T_sub) , function(idents){
    data <- Age_Gene_test_list[[idents]]
    aging_group <- colnames(data)[grep('_comparison', colnames(data))]
    aging_group <- gsub('_VS_Rest_comparison', '', aging_group)
    age_DEGs_list_idents <- lapply(aging_group, function(age){
        #data_iter <- data[data$cluster == age, ]
        tmp_log2FC <- paste0(age, '_VS_Rest_log2FC')
        tmp_p_val_adj <- paste0(age, '_VS_Rest_p_val_adj')
        tmp_p_val <- paste0(age, '_VS_Rest_p_val')
        tmp_pct_1 <- paste0(age, '_VS_Rest_pct_1')
        tmp_pct_2 <- paste0(age, '_VS_Rest_pct_2')
        tmp_variable <- paste0(age, '_VS_Rest_variable')
        sig_DEGs <- data[abs(data[, tmp_log2FC]) > log2(1.5) & #data[, tmp_p_val_adj] < 0.05 &
                        data[, tmp_p_val] * nrow(T_NK_2) < 0.05 & ###手动矫正P值
                        (data[, tmp_pct_1] > 0.1 | data[, tmp_pct_2] > 0.1), tmp_variable]
#         down_DEGs <- data[data[, tmp_log2FC] < -log2(1.5) & #data[, tmp_p_val_adj] < 0.05 &
#                           data[, tmp_p_val]* nrow(T_NK_2) < 0.05 & ###手动矫正P值
#                         (data[, tmp_pct_1] > 0.1 | data[, tmp_pct_2] > 0.1), tmp_variable]
        return(list(sig_DEGs=sig_DEGs))#, down_DEGs=down_DEGs))
    })
    names(age_DEGs_list_idents) <- aging_group
    age_DEGs_list_idents
})
names(Age_DEGs_list) <- levels(T_sub) 

Age_DEGs_list_T <- lapply(names(Age_DEGs_list[[1]]), function(x){
    lapply(Age_DEGs_list, function(y){
        y[[x]]
    })
})
names(Age_DEGs_list_T) <- names(Age_DEGs_list[[1]])

Age_DEGs_shared_number <- lapply(names(Age_DEGs_list_T), function(pair_name){
    pair <- Age_DEGs_list_T[[pair_name]]
    sig_DEGs <- unlist(lapply(pair, function(x){
        x$sig_DEGs
    }))
    sig_DEGs <- table(table(sig_DEGs))
    df <- data.frame(sig_DEGs)
    colnames(df) <- c('count', 'number')
    df$count <- as.numeric(df$count)
    df <- df[order(df$count), ]
    df$pair_name <- pair_name
    df
})
Age_DEGs_shared_number <- do.call(rbind, Age_DEGs_shared_number)

options(repr.plot.width=12, repr.plot.height=8)
p <- ggplot(Age_DEGs_shared_number, aes(count, number))+
geom_bar(stat = 'identity')+
geom_text(mapping = aes(label = number))+
facet_wrap(.~pair_name, ncol =2)
p

min_cutoff <- 5
Age_DEGs_shared <- lapply(names(Age_DEGs_list_T), function(pair_name){
    pair <- Age_DEGs_list_T[[pair_name]]
    sig_DEGs <- unlist(lapply(pair, function(x){
        x$sig_DEGs
    }))
    sig_DEGs <- table(sig_DEGs)
#     down_DEGs <- unlist(lapply(pair, function(x){
#         x$down_DEGs
#     })) 
#     down_DEGs <- table(down_DEGs)
    sig_DEGs_filtered <- names(sig_DEGs)[sig_DEGs >= min_cutoff]
    #down_DEGs_filtered <- names(down_DEGs)[down_DEGs >= min_cutoff]
    return(sig_DEGs_filtered)#, down_DEGs_filtered))
})
names(Age_DEGs_shared) <- names(Age_DEGs_list_T)
Age_genes_selected <- unique(unlist(Age_DEGs_shared))


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

age_colors <- c('Prepuberal'="#b7d9ee", 'Adult'="#3ca88e", 'Aged'="#ec8f46", 'Centenarian' = '#c04e5c')

module_colors <- c('1' = '#be6b63', '2' = '#86a866', 
                            '3' = '#bda25b', '4' = '#6a82aa', '5' ='#6e578c','6' ='#e5cca9', '7' ='#7cc3cd')

pseudocount.use <- 1
gene_exp <- lapply(levels(T_sub), function(clu){
    age_tag <- c('Prepuberal', 'Adult', 'Aged', 'Centenarian')
    result <- lapply(age_tag, function(age){
        exp1 <- rowMeans(expm1(T_sub$RNA@data[Age_genes_selected, Idents(T_sub) %in% 
                                              clu & T_sub$age_group %in% age]))+pseudocount.use
        exp2 <- rowMeans(expm1(T_sub$RNA@data[Age_genes_selected, Idents(T_sub) %in% 
                                              clu & !T_sub$age_group %in% age]))+pseudocount.use
        val <- exp1/exp2
        val <- log2(val)
        val
    })
    result <- do.call(cbind, result)
    colnames(result) <- paste0(age_tag, '~', clu)
    result
})

gene_exp <- do.call(cbind, gene_exp)

gene_exp <- gene_exp[, levels(droplevels(group_var))]

saveRDS(gene_exp, file = 'gene_exp_raw.rds')

gene_exp <- gene_exp_raw <- readRDS('gene_exp_raw.rds')

gene_exp <- my_scale(gene_exp, 'row')

#gene_exp <- my_scale(gene_exp, 'row')


options(repr.plot.width=12, repr.plot.height=25)
hr <- hclust(dist(gene_exp, method = "euclidean"), method="complete") 
gene_exp_rl = cutree(hr, 7)
gene_exp_cl <- factor(sapply(strsplit(colnames(gene_exp), split = '_'), function(x) x[1]))
gene_exp_cl <- factor(gene_exp_cl, levels = unique(gene_exp_cl))                
                    
p <- 
Heatmap(gene_exp, cluster_columns = F, show_row_names = T,use_raster = F,column_split=gene_exp_cl,
        heatmap_legend_param=list(title='Expression'),
        row_split = gene_exp_rl,row_dend_reorder=T,
        row_names_gp = gpar(fontsize = 4, fontfamily = 'sans', fontface = 'italic'), 
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
             ))
                        
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(
    draw(p, merge_legends = T, heatmap_legend_side='right', padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p
#ggsave('PB_aging_DEGs.png', width = 12, height = 25, plot = gg_p, dpi=300) 


####combined module1 and module2
gene_exp_rl <- gene_exp_rl-1
gene_exp_rl[gene_exp_rl == 0] <- gene_exp_rl[gene_exp_rl == 0]+1

enrich_result_list <- lapply(unique(gene_exp_rl), function(i){
    df <- enricher(names(gene_exp_rl)[gene_exp_rl==i], 
                   universe = rownames(T_sub)[rowSums(T_sub$RNA@data) > 0], 
                   TERM2GENE = database)@result
    #pathways <- df[df$p.adjust < 0.05, 'ID']
    #pathways
    df
})

names(enrich_result_list) <- unique(gene_exp_rl)

saveRDS(enrich_result_list, file = 'enrich_result_list.rds')

pathway_list <- lapply(enrich_result_list, function(df){
    df <- df[order(df$p.adjust), ]
    df <- df[df$p.adjust < 0.05, ]
    pathways <- df$ID[1:20]
    pathways <- pathways[!is.na(pathways)]
    pathways <- gsub('^GO_|^KEGG_', '', pathways)
    pathways <- toupper(pathways)
    pathways <- gsub('_', ' ', pathways)
    #pathways <- df[df$p.adjust < 0.05, 'ID']
    #pathways
    pathways
})
names(pathway_list) <- names(enrich_result_list)


pathway_list <- lapply(pathway_list, function(x){
    x <- gsub(',', ',\n', x)
})


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
    width = max_text_width(unlist(pathway_list)) - unit(15, "cm")
))  )
                           
                        
gg_gene_exp_p <- ggplotify::as.ggplot(grid::grid.grabExpr(
    draw(gene_exp_p, merge_legends = T, heatmap_legend_side='right', ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_gene_exp_p
#ggsave('PB_aging_DEGs.png', width = 12, height = 20, plot = gg_gene_exp_p, dpi=200) 
                             

gene_exp_var <- apply(gene_exp_raw, 1, var)

plot_features <- c(
'CCR7',#'CD81',
    'ID3','GLS','LFNG','GIMAP1','GIMAP5',#'TAF10',
    'MAP2K2','IL2RG', 'ETS1', 'BCL11B', 'JUND',
'KLF2','LEF1',#'NCK2','PDE7A',
    'SATB1','SELL',
'TCF7','TGFBR2',#'SUN2',
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
    result <- names(gene_exp_var_filtered)[1:100]
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
'GO_RESPONSE_TO_INTERFERON_GAMMA',
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

gene_exp_filtered <- t(gene_exp[rownames(gene_exp) %in% unlist(gene_exp_rl_filtered), ])
gene_exp_filtered <- gene_exp_filtered[, features_lvls]
new_order <- lapply(row_order(gene_exp_p), function(x){
    intersect(rownames(gene_exp)[x], unlist(gene_exp_rl_filtered))
})
new_order <- new_order[order(names(new_order))]
features_lvls <- unlist(new_order)

plot_features <- intersect(colnames(gene_exp_filtered), plot_features)
text_colors <- module_colors[as.character(gene_exp_rl[intersect(features_lvls, plot_features)])]
row_split <- factor(gene_exp_rl[features_lvls], 
                    levels = as.character(sort(unique(gene_exp_rl[features_lvls]))))






p <- Heatmap(t(gene_exp_filtered),col = colorRamp2(c(-3,-1.5,0,1.5,3), 
                                                c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b')),
             row_title = NULL,
             cluster_columns = F, cluster_rows = F, 
             show_row_names = F,show_column_names = F,show_heatmap_legend=F,
         #row_split = gene_exp_cl, 
        row_split = factor(gene_exp_rl[features_lvls], levels = as.character(unique(gene_exp_rl[features_lvls]))),
        heatmap_legend_param=list(title='FC', direction = "vertical", 
                                  labels_gp = gpar(fontsize = 13.75),title_gp = gpar(fontsize = 13.75)), 
                   top_annotation = HeatmapAnnotation(show_legend = F,show_annotation_name = T,
                                       
               df=data.frame(Age = sapply(strsplit(rownames(gene_exp_filtered), split = '~'), function(x) x[1]), 
                             Identity = sapply(strsplit(rownames(gene_exp_filtered), split = '~'), function(x) x[2])), 
                 col = list(Age = age_colors, Identity = cluster_colors),
                 annotation_name_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = 'bold'),                              
               annotation_legend_param = list(
                 Age = list(

                   title_gp = gpar(fontsize = 13.75, fontfamily = 'sans', fontface = "bold"),
                   labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
                   ncol = 1
                 ), 
               Identity = list(

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
                    #text <- gsub(',', ',\n', text)
                    txt = paste(text, collapse = "\n")
                    #txt = paste0(txt, "\n", length(index), " rows")
                    #txt = paste(index, collapse = ",")
                    grid.text(txt, 0.03, 0.5, rot = 0,hjust =0,
                        gp = gpar(col = module_colors[levels], fontsize = 13))
                },
                width = max_text_width(unlist(pathway_list)) - unit(0.1, "cm")
            )),
left_annotation = rowAnnotation(' ' = anno_mark(link_width = unit(7.5,'mm'),side = 'left',
                                                                          #at_adj = -0.5,
                                                                          #extra_text_args = list(vjust = 0.5),
            at = which(colnames(gene_exp_filtered) %in% plot_features), 
            labels = plot_features, 
            labels_gp = gpar(fontsize=13.75,col = text_colors,#'black',#text_colors, 
                             fontface='italic')), 
           'foo' = anno_empty(border = FALSE, 
                              width=
                              (max_text_width(plot_features) - unit(2.2,'cm')) ))
                                         )
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p,merge_legend = T,  heatmap_legend_side="right",ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
options(repr.plot.width=8.7, repr.plot.height=11)
gg_p





Identity <- unique(sapply(strsplit(rownames(gene_exp_filtered), split = '~'), function(x) x[2]))
lgd1 = Legend(labels =  Identity, title = "Identity", by_row = T,
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 4,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[Identity])))
                          
Age <- unique(sapply(strsplit(rownames(gene_exp_filtered), split = '~'), function(x) x[1]))
lgd2 = Legend(labels =  age_group_lbls_func(Age), title = "Age", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              ncol = 1,title_position = 'topleft',
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
              legend_gp = gpar(col = c(module_colors[Module])))     

col_fun =  colorRamp2(c(-3,-1.5,0,1.5,3), c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))
lgd4 = Legend(col_fun = col_fun, title = "Fold change", direction = 'horizontal',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
             )                     
pd = packLegend(lgd1,lgd2,lgd3, lgd4, direction = "horizontal")

options(repr.plot.width=8, repr.plot.height=1.1)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(0, "cm"), 
                                              y = unit(0, "cm"), just = c("left", "bottom"))))
# ggsave('PB_Tcell_age_DEGs_legend.svg', width = 11.75, height = 1.2) 



save(gene_exp_cl, gene_exp_rl, file = 'gene_exp_rl_cl.RData')






