library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(Matrix)
library(clusterProfiler)
library(openxlsx)

library(aplot)

library(patchwork)

library(circlize)

go <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c5.all.v7.0.symbols.gmt')
kegg <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c2.cp.kegg.v7.0.symbols.gmt')
reactome <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/ReactomePathways.gmt')

database <- rbind(go, kegg, reactome)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/thymus/plot/support//plot_helper.r')


filtered_combined_TEC_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_TEC_2.rds')


levels(filtered_combined_TEC_2) <- c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC',
                            'Immature_TEC', 'mTEC_lo', 'mTEC_hi',
       'post_AIRE_mTEC', 'MKI67+mTEC', 'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated')

filtered_combined_TEC_2$age_group <- droplevels(factor(filtered_combined_TEC_2$Age_group, 
                                                       levels = c('Prepuberal', 'Adult', 'Aged')))


tbl <- table(Idents(filtered_combined_TEC_2),filtered_combined_TEC_2$age_group)

tbl

clusters_filtered <- rownames(tbl)[apply(tbl,1, function(x){any(x<50)})]

####
TEC_sub <- filtered_combined_TEC_2[, !Idents(filtered_combined_TEC_2) %in% c(clusters_filtered)]

group_var <- factor(paste0(TEC_sub$age_group, '~',Idents(TEC_sub)),
levels=unlist(lapply(levels(TEC_sub$age_group), function(x){
    paste0(x, '~', levels(Idents(TEC_sub)))
})))

Age_Gene_test_list <- readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/Integrated_Niche//Age_Gene_test_list_combined.rds')

Age_DEGs_list <- lapply(levels(TEC_sub) , function(idents){
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
                        data[, tmp_p_val] * nrow(filtered_combined_TEC_2) < 0.05 & ###手动矫正P值
                        (data[, tmp_pct_1] > 0.1 | data[, tmp_pct_2] > 0.1), tmp_variable]
        return(list(sig_DEGs=sig_DEGs))
    })
    names(age_DEGs_list_idents) <- aging_group
    age_DEGs_list_idents
})
names(Age_DEGs_list) <- levels(TEC_sub) 

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
#     down_DEGs <- unlist(lapply(pair, function(x){
#         x$down_DEGs
#     })) 
#     down_DEGs <- table(table(down_DEGs))
#     common <- intersect(names(up_DEGs), names(down_DEGs))
#     total <- up_DEGs[common]+down_DEGs[common]
#     total <- c(total, up_DEGs[setdiff(names(up_DEGs), common)], 
#               down_DEGs[setdiff(names(down_DEGs), common)])
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
facet_wrap(.~pair_name, ncol =3)
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


length(Age_genes_selected)

cluster_colors <- c('Fb_1' ='#E5D2DD','cTEC_hi' ='#70a887', 'MKI67+Fb'='#E59CC4', 
                   'Edo'='#F3B1A0', 'Fb_2'='#D6E7A3', 'TEC(neuron)'='#57C3F3', 'mTEC_lo'='#5a829d',#'#476D87',
'VSMCs'='#E95C59', 'mTEC_hi'='#F1BB72','Ionocyte'='#cdb979',
                    'MKI67+mTEC'='#AB3282', 'TEC(myo)'='#7db0bb',#'#91D0BE', 
                    'post_AIRE_mTEC'='#92a5a5', 'Tuft'='#c370ac', #'#3c97c4',
                    'Ionocyte'='#cdb979', 'Tuft/Ionocyte' = '#ff001c',
                    'Ciliated'='#cc5f91', 'cTEC_lo'='#c4c04f', 'MKI67+cTEC'='#ff7e00', 'Immature_TEC'='#9dabd5', 
                    'Mesothelial'='#d73e4b', 'Myelin'='#ff8ae0',#'#639791'
                        'MKI67+VSMCs'='#7ed670', 'MKI67+Endo'='#698c68'


                   )

module_colors <- c('1' = '#be6b63', '2' = '#86a866', 
                            '3' = '#bda25b', '4' = '#6a82aa', '5' ='#6e578c','6' ='#e5cca9', '7' ='#7cc3cd')

pseudocount.use <- 1
gene_exp <- lapply(levels(TEC_sub), function(clu){
    age_tag <- c('Prepuberal', 'Adult', 'Aged')
    result <- lapply(age_tag, function(age){
        exp1 <- rowMeans(expm1(TEC_sub$symbol@data[Age_genes_selected, Idents(TEC_sub) %in% 
                                              clu & TEC_sub$age_group %in% age,drop=F]))+pseudocount.use
        exp2 <- rowMeans(expm1(TEC_sub$symbol@data[Age_genes_selected, Idents(TEC_sub) %in% 
                                              clu & !TEC_sub$age_group %in% age,drop=F]))+pseudocount.use
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

options(repr.plot.width=12, repr.plot.height=25)
hr <- hclust(dist(gene_exp, method = "euclidean"), method="complete") 
gene_exp_rl = cutree(hr, 9)
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


###combined similar module
orig_gene_exp_rl <- gene_exp_rl
gene_exp_rl[orig_gene_exp_rl %in% c(4)] <- 1
gene_exp_rl[orig_gene_exp_rl %in% c(6)] <- 2
gene_exp_rl[orig_gene_exp_rl %in% c(5)] <- 3
gene_exp_rl[orig_gene_exp_rl %in% c(1,2,8)] <- 4
gene_exp_rl[orig_gene_exp_rl %in% c(3,7,9)] <- 5
# gene_exp_rl[orig_gene_exp_rl %in% c(12)] <- 3





enrich_result_list <- lapply(unique(gene_exp_rl), function(i){
    df <- enricher(names(gene_exp_rl)[gene_exp_rl==i], 
                   universe = rownames(TEC_sub)[rowSums(TEC_sub$RNA@data) > 0], 
                   TERM2GENE = database)@result
    #pathways <- df[df$p.adjust < 0.05, 'ID']
    #pathways
    df
})

names(enrich_result_list) <- unique(gene_exp_rl)

pathway_list <- lapply(enrich_result_list, function(df){
    df <- df[order(df$p.adjust), ]
    df <- df[df$p.adjust < 0.1, ]
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


saveRDS(enrich_result_list, file = 'enrich_result_list.rds')

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
    width = max_text_width(unlist(pathway_list)) - unit(5, "cm")
))  )
                           
                        
gg_gene_exp_p <- ggplotify::as.ggplot(grid::grid.grabExpr(
    draw(gene_exp_p, merge_legends = T, heatmap_legend_side='right', ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_gene_exp_p
#ggsave('TEC_aging_DEGs.png', width = 12, height = 20, plot = gg_gene_exp_p, dpi=200) 
                             

gene_exp_var <- apply(gene_exp_raw, 1, var)


gene_exp_rl_filtered <- lapply(unique(gene_exp_rl), function(iter){
    gene_exp_rl_iter <- gene_exp_rl[gene_exp_rl==iter]
    gene_exp_var_filtered <- gene_exp_var[names(gene_exp_rl_iter)]
    gene_exp_var_filtered <- gene_exp_var_filtered[order(gene_exp_var_filtered, decreasing = T)]
    result <- names(gene_exp_var_filtered)[1:100]
    result <- union(result, intersect(c('HLA-C','HLA-DPB1','HLA-DQA1','PSMA4',#'HLA-DRA',
                   'CCL25','MT-ND1',#'MT-ND6',
'MT-CO2',#'MT-CO3',
                   'ZFP36L1','MYC','CDKN1A','FOSB',#'COPE',
                   'PSMB10','PSMB9','KRT1','TXNIP',
#'HIGD1A',
                   'NENF',
                   'FABP4','CD81','EGR1','CEBPB','CEBPD',#'ARID1B',
                   'TAF10','PTEN','TGFBR3','FKBP1A','RHOA',
'UBC','ABCD4','PEX13','PEX2','WNT5B', 'KRT17', #'TUBB2A',
                   'NFKB2', 'ALDOA', 'UBB', 'ZNF703'), 
                                      names(gene_exp_var_filtered)))
    result <- result[!is.na(result)]
    result
})

pathway_list <- list('1' = c('Interferon alpha/beta signaling',
'GO_MHC_CLASS_II_PROTEIN_COMPLEX',
'GO_ANTIGEN_PROCESSING_AND_PRESENTATION',
'GO_TRANSPORT_VESICLE',
'GO_PROTEASOME_CORE_COMPLEX'),
                     '2' =c('GO_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING',
                            'GO_ANTIGEN_PROCESSING_AND_PRESENTATION',
'Mitotic Anaphase','MAPK6/MAPK4 signaling','Glycolysis'#,'DNA Replication'
                           ),
                            
                              '3' = c('rRNA processing in the mitochondrion',
'KEGG_OXIDATIVE_PHOSPHORYLATION',
'GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
'GO_RESPONSE_TO_OXYGEN_LEVELS'
#'GO_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX'
                                     ),
                             '4' = c('KEGG_PEROXISOME',
'GO_FATTY_ACID_CATABOLIC_PROCESS',
'GO_LIPID_OXIDATION') , '5' =c('Signaling by WNT',
'PTEN Regulation',
# 'GO_UBIQUITIN_LIKE_PROTEIN_TRANSFERASE_ACTIVITY',
# 'GO_COVALENT_CHROMATIN_MODIFICATION',
'GO_MESENCHYMAL_CELL_DIFFERENTIATION',
'GO_EPITHELIAL_TO_MESENCHYMAL_TRANSITION',
'GO_SMAD_BINDING',
'GO_REGULATION_OF_BMP_SIGNALING_PATHWAY',
'GO_FAT_CELL_DIFFERENTIATION',
'GO_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA')

                   )



features_lvls <- intersect(rownames(gene_exp)[unlist(row_order(gene_exp_p))], unlist(gene_exp_rl_filtered))

plot_features <- c('HLA-C','HLA-DPB1','HLA-DQA1','PSMA4',#'HLA-DRA',
                   'CCL25','MT-ND1','MT-ND6',#'MT-ND6',
'MT-CO2',#'MT-CO3',
                   'ZFP36L1','MYC',#'CDKN1A',
                   'FOSB',#'COPE',
                   'PSMB10','PSMB9','KRT1',#'TXNIP',
#'HIGD1A',
                   'NENF',
                   'FABP4','CD81','EGR1','CEBPB','CEBPD',#'ARID1B',
                   'TAF10','PTEN','TGFBR3','FKBP1A','RHOA',
'UBC','ABCD4','PEX13','PEX2','WNT5B', 'KRT17', #'TUBB2A',
                   'NFKB2', 'ALDOA', 'UBB', 'ZNF703')



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
                    #text <- gsub(',', ',\n', text)
                    txt = paste(text, collapse = "\n")
                    #txt = paste0(txt, "\n", length(index), " rows")
                    #txt = paste(index, collapse = ",")
                    grid.text(txt, 0.03, 0.5, rot = 0,hjust =0,
                        gp = gpar(col = module_colors[levels], fontsize = 13.75))
                },
                width = max_text_width(unlist(pathway_list)) - unit(3, "cm")
            )), 
            left_annotation = rowAnnotation(' ' = anno_mark(link_width = unit(2.5,'mm'),side = 'left',
                                                                          #at_adj = -0.5,
                                                                          #extra_text_args = list(vjust = 0.5),
            at = which(rownames(gene_exp_filtered) %in% plot_features), 
            labels = plot_features, 
            labels_gp = gpar(fontsize=13.75,col = text_colors,#'black',#text_colors, 
                             fontface='italic')), 
           'foo' = anno_empty(border = FALSE, 
                              width=
                              (max_text_width(plot_features) - unit(2.2,'cm')) )))
options(repr.plot.width=9, repr.plot.height=12.5)
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p, heatmap_legend_side="top", annotation_legend_side="right",legend_grouping = "original",#ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p


Identity <- unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[2]))
lgd1 = Legend(labels =  Identity, title = "Identity", by_row = T,
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 3,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[Identity])))
                          
Age <- age_group_lbls_func(unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[1])))
lgd2 = Legend(labels =  Age, title = "Age", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 3,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(age_colors[Age])))

                     
Module <- sort(unique(gene_exp_rl[features_lvls]))
lgd3 = Legend(labels =  Module, title = "Module", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 3,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(module_colors[as.character(Module)])))     

col_fun =  colorRamp2(c(-3,-1.5,0,1.5,3), c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))
lgd4 = Legend(col_fun = col_fun, title = "Fold change", direction = 'horizontal',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
             )                                          

pd = packLegend(lgd1,lgd2,lgd3,
                lgd4, direction = "horizontal")

options(repr.plot.width=8, repr.plot.height=1)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(0, "cm"), 
                                              y = unit(0, "cm"), just = c("left", "bottom"))))
#ggsave('TEC_age_DEGs_legend.svg', width = 12, height = 1) 

save(gene_exp_rl, gene_exp_cl, file = 'gene_exp_rl_cl.RData')


