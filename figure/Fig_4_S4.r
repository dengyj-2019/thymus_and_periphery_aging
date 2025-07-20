library(ggplot2)
library(Seurat)
library(MySeuratWrappers)
library(ComplexHeatmap)

library(pheatmap)

library(ggfun)

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

library(rstatix)


plot_dir <- '/data1/02.private/dengyj/analysis/thymus/Hu_plot/Fig4_FigS4'

setwd(plot_dir)



lbls_func <- function(x){
    x <- gsub('post_AIRE_mTEC', 'post-AIRE mTEC', x)
    x <- gsub('Edo', 'Endo', x)
    x <- gsub('\\(', '-', x)
    x <- gsub('\\)', '', x)    
    x <- gsub('_', '-', x)
    x <- gsub('\\+', ' ', x)
    x <- gsub('Fb', 'FB', x)
    x
}

cluster_colors <- c('Fb_1' ='#E5D2DD','cTEC_hi' ='#70a887', 'MKI67+Fb'='#E59CC4', 
                   'Endo'='#F3B1A0', 'Fb_2'='#D6E7A3', 'TEC(neuron)'='#57C3F3', 'mTEC_lo'='#5a829d',#'#476D87',
'VSMCs'='#E95C59', 'mTEC_hi'='#F1BB72','Ionocyte'='#cdb979',
                    'MKI67+mTEC'='#AB3282', 'TEC(myo)'='#7db0bb',#'#91D0BE', 
                    'post_AIRE_mTEC'='#92a5a5', 'Tuft'='#c370ac', #'#3c97c4',
                    'Ionocyte'='#cdb979', 'Tuft/Ionocyte' = '#ff001c',
                    'Ciliated'='#cc5f91', 'cTEC_lo'='#c4c04f', 'MKI67+cTEC'='#ff7e00', 'Immature_TEC'='#9dabd5', 
                    'Mesothelial'='#d73e4b', 'TEC'='#ae80dc',#'#639791'
                        'MKI67+VSMCs'='#7ed670', 'MKI67+Endo'='#698c68', 'Myelin'='#ff8ae0',
                    'ETP' = '#297A6B','TP1'='#C05EA5',
    'TP2'='#e6ba77'


                   )

module_colors <- c('1' = '#be6b63', '2' = '#86a866', 
                            '3' = '#bda25b', '4' = '#6a82aa', '5' ='#6e578c','6' ='#e5cca9', '7' ='#7cc3cd')

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

source('/data1/02.private/dengyj/analysis/mycode/crosstalk/V1.04/calculation.R')


###avoid conflict with add_label function
detach('package:psych')



filtered_combined_Niche_2  <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_Niche_2.rds')




filtered_combined_Niche_2 <- filtered_combined_Niche_2[, !Idents(filtered_combined_Niche_2) %in% c('Mesothelial', 'Myelin')]


levels(filtered_combined_Niche_2) <- 
c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC', 'Immature_TEC', 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC', 
    'MKI67+mTEC',
'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated', 'Fb_1', 'Fb_2', 'MKI67+Fb', 'VSMCs', 'MKI67+VSMCs',
                                            'Endo', 'MKI67+Endo', 'Mesothelial', 'Myelin')

###Fig.4a
plot_df <- cbind( Embeddings(filtered_combined_Niche_2, reduction = 'umap'), 
                 data.frame(Identity = Idents(filtered_combined_Niche_2)))


new_label <- as.character(1:nlevels(filtered_combined_Niche_2))
names(new_label) <- levels(filtered_combined_Niche_2)


options(repr.plot.width=6, repr.plot.height=4)
p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = Identity))+
geom_point(size= 0.5)+
#guides(color = guide_legend(override.aes = list(size = 7), ncol = 1))+
theme(legend.background = element_blank(),
      panel.background = element_rect(color = NA, fill = 'white'), 
      panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.ticks=element_blank(),legend.title=element_text(size = 13.75,margin = margin(b = -0.2,unit = 'cm')),
      axis.text=element_blank(),axis.line=element_blank(), axis.title=element_text(size = 13.75, color = 'black'),
      legend.text = element_text(size=13.75, family= 'sans', color = 'black',
                                 margin = margin(t = -0.03,unit = 'cm')))+
scale_color_manual(values = cluster_colors, label = lbls_func)
p <- add_label(p, use_new_label = T, new_label = c(new_label),add_circle_to_label = F,add_bar_clu = c(),

         label_arg = list(color='black', size=13.75/.pt, fontface = 'plain'))+
  guides(color=guide_legend(override.aes = list(size=1.4, label = c(new_label), alpha=1), ncol = 1))+
theme(legend.key.height=unit(0,"mm"),legend.key=element_rect(fill= 'white'))
p
my_saveplot_fun(p, 'filtered_combined_Niche_2_idents', width = 6, height = 4)

###Extended Data Fig.4a
p <- DotPlot(filtered_combined_Niche_2, 
             features = rev(unique(c(#'SIX1',
                                     'KRT5','PAX9','KRT8','CCL25','PSMB11','HLA-DQB1','PRSS16',
                       'FOXN1',#'CXCL12','ZBED2',
                       'PAX1','DLL4','MKI67',#'CCNA2',
                 'CDK1',
'MAOA','DPYS',         'EPCAM','KRT14','CCL19','KRT15','LYPD1',
'AIRE','FEZF2','SPIB','FXYD3','IVL','KRT1','POU2F3','OVOL3',
                                 'FOXI1','ASCL3','MYOG','DES',#'MYLPF',
                 'BEX1','NEUROD1',
                                 'ATOH1','GFI1',
 'C7',
                                 'PDGFRA',  #'MFAP5',
              'FBN1', 'RGS5', 'COX4I2', 
                                 'CD36',
              'CD300LG',
              'CLDN5','CDH5'))),
       cols = c('lightgrey', 'red'))

new_p <- ggplot(p$data, aes(features.plot,id , fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21)+
scale_size_continuous(range = c(0,5))+
scale_fill_gradientn(colours = c('white', 'red')) + 
d_theme_w(size = 13)+
theme(#axis.text.y = element_text(size = 13.75, family = 'sans'), 
      plot.margin = margin(l=0.5, unit = 'cm'),
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'horizontal',
      legend.position = 'top',
        legend.margin = margin(l = 0, r = 0, unit = 'cm'),      
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.x = element_text( family = 'sans', face = 'italic',vjust = 1, angle =60, hjust = 1),
      #legend.text = element_text(size = 13.75, family = 'sans'),
      #legend.title = element_text(size = 13.75, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 13.75, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank())+
guides(
       fill = guide_colorbar(title = "Expression", order = 1, label.position = "bottom", 
            title.position = "left", title.vjust=0.9),
size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'), order = 2, 
                    label.position = "bottom", title.vjust=0.9 ))+
scale_y_discrete(limits = rev, position = 'right', 
                 label = lbls_func)
options(repr.plot.width=11, repr.plot.height=5)
new_p
ggsave('niche_DEG.svg', width = 11, height = 5)

###Extended Data Fig.4b
filtered_combined_Niche_2$idents <- as.character(Idents(filtered_combined_Niche_2))
filtered_combined_Niche_2$idents[filtered_combined_Niche_2$idents %in% 
                                c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC', 'Immature_TEC', 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC', 'MKI67+mTEC','Tuft', 'Ionocyte',  'TEC(myo)', 'TEC(neuron)', 'Ciliated')] <- 'TEC'
filtered_combined_Niche_2$idents <- factor(filtered_combined_Niche_2$idents, 
                                           levels= c('TEC','Fb_1', 'Fb_2', 'MKI67+Fb', 'VSMCs', 'MKI67+VSMCs',
                                            'Endo', 'MKI67+Endo', 'Mesothelial', 'Myelin'))

logi <- !filtered_combined_Niche_2$sort %in% c('EPCAM') & !filtered_combined_Niche_2$batch %in% 
unique(filtered_combined_Niche_2$batch[grepl('Aged|Adult|CD45', filtered_combined_Niche_2$batch)]) & 
!filtered_combined_Niche_2$idents %in% c('Mesothelial', 'Myelin')
plot_df = plyr::ddply(data.frame(
    table(droplevels(filtered_combined_Niche_2$idents[logi]), 
          filtered_combined_Niche_2$Age_group[logi])), 
                      'Var2', transform,percentage=Freq/sum(Freq) *100)

colnames(plot_df) <- c('Cluster', 'Age', 'Freq', 'Percentage')

plot_df$Age <- factor(plot_df$Age, levels = c('Prepuberal', 'Adult', 'Aged'))



p <- ggplot(plot_df,aes(x=Age,y=Percentage,fill=Cluster,group=Cluster))+
geom_area(colour="black",size=0.1)+
scale_fill_manual(values = cluster_colors, label = lbls_func)+
#ggtitle("Proportion of thymocyte during aging ")+
#theme_bw()+
theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),
     plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 18,face = "italic"),
    axis.title.x=element_text(size=13.75),axis.title.y=element_text(size=13.75),
    axis.text.x=element_text(size=13.75, angle = 20, hjust =1),axis.text.y=element_text(size=13.75),
    legend.text=element_text(size=13.75),legend.title=element_text(size=13.75),
)+
labs(x = NULL, fill = 'Identity', y = 'Percentage %')+
guides(fill = guide_legend(ncol =1))+
scale_x_discrete( expand=c(0,0), label = age_group_lbls_func)+
scale_y_continuous( expand = c(0, 0))#+
#labs(title = bquote(paste('excluding EpCAM'^'+', ' cells')))

options(repr.plot.width=5, repr.plot.height=5)
p
#my_saveplot_fun(p,'filtered_combined_Niche_2_proportion', width = 6, height = 6)
ggsave('filtered_combined_Niche_2_proportion.svg', width = 5, height = 5)





filtered_combined_TEC_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_TEC_2.rds')




###Extended Data Fig.4c
tbl <- table(Idents(filtered_combined_TEC_2),filtered_combined_TEC_2$Age_group)####Statistic of cell number in each age-group


tmp_tbl <- tbl[rownames(tbl)[apply(tbl,1, function(x){all(x>50)})],]
colnames(tmp_tbl) <- paste0('Group ', as.roman(1:3))


logi <- Idents(filtered_combined_TEC_2) %in% rownames(tmp_tbl)

options(repr.plot.width=5.5, repr.plot.height=6)
plot_df = data.frame(table(Idents(filtered_combined_TEC_2)[logi], 
                                       filtered_combined_TEC_2$batch[logi]))
colnames(plot_df) <- c('Identity', 'batch', 'Freq')
plot_df$Age_group <- 
filtered_combined_TEC_2$Age_group[match(plot_df$batch, filtered_combined_TEC_2$batch)]
plot_df <- plot_df[plot_df$Freq > 0, ]
plot_df <- plyr::ddply(plot_df, c('Identity', 'Age_group'), transform,percentage=Freq/sum(Freq) *100)

color_used <- c('#E5D2DD', '#53A85F', '#F1BB72',  '#D6E7A3', '#57C3F3', '#476D87',
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

plot_df$Identity <- factor(plot_df$Identity, levels= intersect(levels(filtered_combined_TEC_2), 
                                                               unique(plot_df$Identity)))
plot_df$Age_group[plot_df$Age_group=='Prepuberal'] <- 'Group I'
plot_df$Age_group[plot_df$Age_group=='Adult'] <- 'Group II'
plot_df$Age_group[plot_df$Age_group=='Aged'] <- 'Group III'
plot_df$batch <- as.character(plot_df$batch)
logi <- grepl('^A43|^T03|^T07|^A16|^T06', plot_df$batch)
plot_df$batch[logi] <- gsub('_.*','', plot_df$batch[logi])
logi <- grepl('GSM4466784|GSM4466785', plot_df$batch)
plot_df$batch[logi] <- 'Parent_10M'
logi <- grepl('Aged_.*\\d', plot_df$batch)
plot_df$batch[logi] <- 'donor40'
logi <- grepl('Aged.*-03$', plot_df$batch)
plot_df$batch[logi] <- 'donor41'
p <- ggplot(plot_df,aes(x=Identity,y=percentage,
                   fill=batch,
                   group=batch))+
geom_bar(stat='identity')+
#scale_fill_manual(values = cluster_colors, label = lbls_func)+
#ggtitle("Proportion of thymocyte during aging ")+
theme(strip.background = element_blank(),
      panel.border = element_blank(),#element_rect(color = 'black', fill = NA), legend.position = 'right',
      panel.background = element_blank(),#element_rect(color = NA, fill = 'white'), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_blank(),
     plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 18,face = "italic"),
    axis.title.x=element_text(size=13.75, family='sans'),
    axis.title.y=element_text(size=13.75, angle =90, family='sans'),
    axis.text.x=element_text(size=13.75, family='sans', angle=45, hjust=1),
    axis.text.y=element_text(size=13.75, family='sans'),
    legend.text = element_text(size = 13.75, family='sans'),legend.spacing.y = unit(1e-3, 'npc'),
    legend.title = element_text(size = 13.75, family='sans')
)+guides(fill = guide_legend(ncol = 1, byrow = TRUE))+
scale_x_discrete(expand=c(0,0), label =lbls_func)+
scale_y_continuous( expand = c(0, 0))+
facet_wrap(.~Age_group)+
labs(x = NULL, y = 'Percentage %')+
guides(fill =F)+
scale_fill_manual(values = color_used)+
d_theme_w(size =13.75, italic_text = '')+
theme(plot.title = element_text(size = 20, face = 'bold'))+
labs(title = 'Filled by donors')
options(repr.plot.width=6.5, repr.plot.height=5)
p
ggsave('TEC_donor_filled.svg', width = 6.5, height = 5)

###Extended Data Fig.4d

levels(filtered_combined_TEC_2) <- 
c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC', 'Immature_TEC', 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC', 
    'MKI67+mTEC',
'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated')


num_of_sample <- table(filtered_combined_TEC_2$Age_group)
plot_df <- table(Idents(filtered_combined_TEC_2),filtered_combined_TEC_2$Age_group)


plot_list <- lapply(c('Prepuberal', 'Adult', 'Aged'), function(i){
    df <- data.frame(number = plot_df[, i], idents = factor(rownames(plot_df), 
                                                            levels =levels(filtered_combined_TEC_2))) 
    p <- mypie(data = df, area_var = 'number', lbl_size = 13.75/.pt,
         fill_var = 'idents', #labels = lbls,
         plot_label = T, lbl_adjust = 0.8,digits = 1,
         pie_colors = 'black', lbl_face = 'bold', pie_size = 0.1, min_cutoff=0.03) + 
    theme(legend.position = 'right', legend.direction = "horizontal",
          #plot.title = element_text(hjust = 0.5 ,size = 10), 
          axis.text.x = element_blank(), legend.text = element_text(size = 13.75), 
          legend.title = element_text(size = 13.75), 
         plot.title = element_text(size = 13.75,hjust=0.5))+
    scale_fill_manual(values =cluster_colors, label = lbls_func)+
    guides(fill = guide_legend(nrow = 3))+
    labs(fill = NULL, title = age_group_lbls_func(i))
    return(p)
})


options(repr.plot.width=12, repr.plot.height=5)
wrap_plots(plot_list)+plot_layout(ncol = 3, guides = 'collect')&
theme(legend.position = 'bottom')
ggsave('filtered_combined_TEC_2_pct.svg', width = 12, height = 5)
#my_saveplot_fun(p,'filtered_combined_TEC_3_pct', width = 5.5, height = 6)







gene_exp <- gene_exp_raw <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Gene/Integrated_Niche/gene_exp_raw.rds')

gene_exp <- my_scale(gene_exp, 'row')

gene_exp_var <- apply(gene_exp_raw, 1, var)

load('/data1/02.private/dengyj/analysis/thymus/Gene/Integrated_Niche/gene_exp_rl_cl.RData')

####这里是按照module区分，但是只有2个module对应着Prepuberal和Adult组，所以基本上和按照年龄区分的结果差不多
gene_exp_rl_filtered <- lapply(unique(gene_exp_rl), function(iter){
    gene_exp_rl_iter <- gene_exp_rl[gene_exp_rl==iter]
    gene_exp_var_filtered <- gene_exp_var[names(gene_exp_rl_iter)]
    gene_exp_var_filtered <- gene_exp_var_filtered[order(gene_exp_var_filtered, decreasing = T)]
    result <- names(gene_exp_var_filtered)[1:50]
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
                            
                              '3' = c('tRNA processing',#'rRNA processing in the mitochondrion',
                                      'Respiratory electron transport',
'KEGG_OXIDATIVE_PHOSPHORYLATION',
#'GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
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

###Extended Data Fig.4e
features_lvls <- intersect(rownames(gene_exp)[unlist(row_order(gene_exp_p))], unlist(gene_exp_rl_filtered))

plot_features <- c('HLA-C','HLA-DPB1','HLA-DQA1','PSMA4',
                   'CCL25','MT-ND1','MT-ND6',
'MT-CO2',
                   'ZFP36L1','MYC',
                   'FOSB',
                   'PSMB10','PSMB9','KRT1',

                   'NENF',
                   'FABP4','CD81','EGR1','CEBPB','CEBPD',
                   'TAF10','PTEN','TGFBR3','FKBP1A','RHOA',
'UBC','ABCD4','PEX13','PEX2','WNT5B', 'KRT17', 
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
                    
                    text <- gsub('Pten', 'PTEN', text)  
                    text <- gsub('NF KAPPAB', 'NF-κB', text, ignore.case = T)
                    text <- gsub('KAPPAB', 'κB', text, ignore.case = T)
                    #text <- gsub(',', ',\n', text)
                    txt = paste(text, collapse = "\n")
                    #txt = paste0(txt, "\n", length(index), " rows")
                    #txt = paste(index, collapse = ",")
                    grid.text(txt, 0.03, 0.5, rot = 0,hjust =0,
                        gp = gpar(col = module_colors[levels], fontsize = 13.75))
                },
                width = max_text_width(unlist(pathway_list)) - unit(4, "cm")
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
options(repr.plot.width=8, repr.plot.height=8)
gg_p <- ggplotify::as.ggplot(grid::grid.grabExpr(draw(p, heatmap_legend_side="top", annotation_legend_side="right",legend_grouping = "original",#ht_gap = unit(c(5), 'cm'),
padding = unit(c(0, 0, 0, 0), "mm"))))
gg_p
ggsave('TEC_aging_DEGs.svg', width = 8, height = 8)

###Extended Data Fig.4e legend
Identity <- unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[2]))
lgd1 = Legend(labels =  Identity, title = "Identity", by_row = T,
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 5,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[Identity])))
                          
Age <- age_group_lbls_func(unique(sapply(strsplit(colnames(gene_exp_filtered), split = '~'), function(x) x[1])))
lgd2 = Legend(labels =  Age, title = "Age", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 5,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(age_colors[Age])))

                     
Module <- sort(unique(gene_exp_rl[features_lvls]))
lgd3 = Legend(labels =  Module, title = "Module", 
              title_gp = gpar(fontsize = 13.75, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(7, "mm"),
              nrow = 5,title_position = 'topleft',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(module_colors[as.character(Module)])))     

col_fun =  colorRamp2(c(-3,-1.5,0,1.5,3), c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b'))
lgd4 = Legend(col_fun = col_fun, title = "Fold change", direction = 'vertical',
              labels_gp = gpar(fontsize = 13.75, fontfamily = 'sans'), 
              title_gp = gpar(fontsize = 13.75, fontfamily = 'sans'),
             )                                          

pd = packLegend(lgd1,lgd2,lgd3,
                lgd4, direction = "horizontal")

options(repr.plot.width=6.5, repr.plot.height=1.3)


ggplotify::as.ggplot(grid::grid.grabExpr(draw(pd,x = unit(0, "cm"), 
                                              y = unit(0, "cm"), just = c("left", "bottom"))))
ggsave('TEC_age_DEGs_legend.svg', width = 6.5, height = 1.3      ) 







filtered_TRA_TEC <- readRDS('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/filtered_TRA_TEC_proccess.rds')

TRAs_symbol=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_symbol.rds')

TRAs_df=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_df.rds')

TRAs_tau=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_tau.rds')

TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')



plot_df <- cbind(Embeddings(filtered_TRA_TEC, reduction = 'umap'), idents = Idents(filtered_TRA_TEC),
                 FetchData(filtered_TRA_TEC, vars = c('TRAs_consensus', 'Age_group')))

max_value <- quantile(plot_df$TRAs_consensus, 0.9)
min_value <- quantile(plot_df$TRAs_consensus, 0.1)
plot_df$TRAs_consensus[plot_df$TRAs_consensus > max_value] <- max_value
plot_df$TRAs_consensus[plot_df$TRAs_consensus < min_value] <- min_value

plot_df$Age_group <- factor(plot_df$Age_group, levels = c('Prepuberal','Adult','Aged'))
plot_list <- lapply(levels(plot_df$Age_group), function(x){
    p <- ggplot(plot_df[plot_df$Age_group==x, ], aes(UMAP_1, UMAP_2, color = TRAs_consensus))+
    geom_point(size = 0.5)+
    scale_color_gradientn(colors = c('lightgrey', 'red'))+
    #facet_wrap(.~Age_group,ncol =3)+
    theme(
          strip.background.x = element_blank(),#element_rect(fill = 'white', color = NA),
          panel.background = element_rect(color = NA, fill = 'white'), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
          axis.line = element_blank(), axis.ticks = element_blank(), 
          #axis.text.y = element_text(size = 30),
          legend.text = element_text(size = 13.75),
          legend.title =element_text(size = 13.75),
          axis.title = element_blank(),#element_text(size = 13.75),
          strip.text.x = element_blank(),#element_text(size = 13.75),
          axis.text = element_blank())+
    labs(color = 'TRAs Score')+
    scale_color_gradientn(colours = c('#2948a4', '#5eb4d2', '#e4f1cf', '#ff7443', '#bc002b')) +
    annotation_custom(grid::textGrob(label = age_group_lbls_func(x),hjust = 0.5,
                               x=grid::unit(0.75,"npc") ,                                
                               y=grid::unit(0.92,"npc"), 
                            gp=grid::gpar(col = 'black', fontface = 'plain', fontsize = 13.75)))
    p
})

###Fig.4b
p1 <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = idents))+
    geom_point(size = 0.5)+
    scale_color_gradientn(colors = c('lightgrey', 'red'))+
    #facet_wrap(.~Age_group,ncol =3)+
scale_color_manual(values =cluster_colors)+d_theme_w(size = 13.75)+
labs(x = NULL, y = NULL, color = 'Identity')+
theme( axis.ticks = element_blank(), 
          
          axis.text = element_blank())+
guides(color = guide_legend(override.aes = list(size = 5)))
options(repr.plot.width=6, repr.plot.height=3.75)
p1
#plot_list <- c(list(p1), plot_list)
my_saveplot_fun(p1,'TRA_TEC_idents', width = 6, height = 3.75)

###Fig.4c
options(repr.plot.width=4*3, repr.plot.height=4)
p <- wrap_plots(plot_list) + plot_layout(ncol = 3, guides = 'collect')&
theme(legend.position = 'right')
p
my_saveplot_fun(p,'TRA_score', width = 4*3, height = 4)








loom <- open_loom('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/filtered_TRA_TEC_auc.loom', mode="r")
regulons_incidMat <- get_regulons(loom, attrName = 'Regulons')
filtered_TRA_TEC_regulons <- regulonsToGeneLists(regulons_incidMat)
filtered_TRA_TEC_regulonsAUC <- get_regulonsAuc(loom,attrName = 'RegulonsAUC')
close_loom(loom)
filtered_TRA_TEC_adj <- read.csv('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/filtered_TRA_TEC_adj.csv')

skeletal_muscle_TRAs <- TRAs_consensus_df$symbol[TRAs_consensus_df$tau >=0.8 & 
                                                 grepl('skeletal.*muscle', TRAs_consensus_df$tissue)]

DefaultAssay(filtered_TRA_TEC) <- 'symbol'
markers_TEC_myo <- readRDS('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/markers_TEC_myo.rds')

TEC_myo_TFs <- intersect(rownames(markers_TEC_myo)[markers_TEC_myo$p_val_adj < 0.01 & 
                                                   markers_TEC_myo$avg_logFC > log(1.5)], TF)







###Extended Data Fig.4f
TRAs_consensus_df$tissue <- gsub(' ', '_', TRAs_consensus_df$tissue)
TRAs_consensus_df$tissue <- gsub('\\-', '_', TRAs_consensus_df$tissue)
TRAs_consensus_df$tissue <- gsub(',', '', TRAs_consensus_df$tissue)

tissue <- unique(TRAs_consensus_df$tissue)



options(repr.plot.width=10, repr.plot.height=5)
levels(filtered_TRA_TEC) <-  c( 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC','MKI67+mTEC', 'Tuft', 'Ionocyte','TEC(myo)', 
                            'TEC(neuron)', 'Ciliated')
p <- DotPlot(filtered_TRA_TEC, features = tissue, cols = c('lightgrey', 'red'))+
#coord_flip()+
theme(axis.text.x =element_text(angle=45, hjust=1))+
guides(color =guide_colourbar(title='TRAs'), size = guide_legend(title='Cell %'))
#p
plot_mat <- reshape2::dcast(p$data[, c('id', 'features.plot', 'avg.exp.scaled')], 
                            as.formula(paste0('`features.plot`','~','id')),  
                            value.var = 'avg.exp.scaled') 
rownames(plot_mat) <- plot_mat[,'features.plot']
plot_mat <- plot_mat[, !colnames(plot_mat) %in% 'features.plot']
plot_mat <- order_heatmap(plot_mat)

new_p <- ggplot(p$data, aes(id, features.plot, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21)+
scale_size_continuous(range = c(0,4.5))+
scale_fill_gradientn(colours = c('#4feb11', '#ea9833','#ec261e')) + 
theme(plot.margin = margin(l = 1, unit = 'cm'), 
      axis.text.y = element_text(size = 13.75, family = 'sans', angle =0, hjust = 1,vjust=0.5), 
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      #legend.box = 'horizontal',
      legend.position = 'right',
      #axis.ticks = element_blank(), 
      axis.title = element_text(size = 13.75), 
     axis.text.x = element_text(size = 13.75, family = 'sans',angle =45, hjust=1),
      legend.text = element_text(size = 13.75, family = 'sans'),
      legend.title = element_text(size = 13.75, family = 'sans'),

      axis.line =element_blank())+
guides(
       fill = guide_colorbar(title = "TRAs Score", order = 1),
size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'), order = 2))+
scale_x_discrete(limits = c( 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC','MKI67+mTEC', 'Tuft', 'Ionocyte','TEC(myo)', 
                            'TEC(neuron)', 'Ciliated'), 
                label = function(x){
    x <- gsub('post_AIRE_mTEC', 'post-AIRE mTEC', x)
    x <- gsub('Edo', 'Endo', x)                    
    x <- gsub('\\(', '-', x)
    x <- gsub('\\)', '', x)    
    x <- gsub('_', '-', x)
    x
})+
scale_y_discrete(limits=rev(rownames(plot_mat)),position = 'right',label = function(x){

    x <- gsub('_', ' ', x)
    x
})+
labs(x = NULL, y = NULL)




options(repr.plot.width=7.5, repr.plot.height=9.3)
new_p
ggsave('TRA_exp.svg', width = 7.5, height = 9.3)



###Extended Data Fig.4g
plot_df <- lapply(TEC_myo_TFs[paste0(TEC_myo_TFs, '(+)') %in% names(filtered_TRA_TEC_regulons)], function(TF){
    df <- data.frame(Freq = 
                     sum(skeletal_muscle_TRAs %in% filtered_TRA_TEC_regulons[[paste0(TF,'(+)')]])/length(skeletal_muscle_TRAs))
    df$TF <- TF
    df
})
plot_df <- do.call(rbind, plot_df)
plot_df <- plot_df[order(plot_df$Freq, decreasing = T),]
plot_df$TF <- factor(plot_df$TF, levels = unique(plot_df$TF))

p <- ggplot(plot_df, aes(TF, Freq*100, fill = TF))+
geom_bar(stat='identity')+
labs(y='Proprotion of regulated skeletal_muscle TRAs %', x =NULL, fill ='TFs')+
d_theme_w(size=13.75)+
theme(axis.text.x = element_text(face='italic', angle=60,hjust=1,vjust=1),
      legend.text = element_text(face='italic'),
               panel.border = element_rect(color ='black', fill =NA),
          panel.background= element_rect(fill ='white', color =NA),
      plot.margin = margin(t = 0.8,unit = 'mm'),#axis.ticks.x = element_blank(),
    text=element_text(size=13.75), legend.spacing.y = unit(0, 'cm'),
      
      legend.key.height = unit(2,'mm'),
      legend.key.width = unit(4,'mm'),
      axis.text.y=element_text(size=13.75), 

     legend.title=element_text(size=13.75), axis.title=element_text(size=13.75))+
guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2)))
options(repr.plot.width=6, repr.plot.height=5)
p
#my_saveplot_fun(p, 'TRAs_regulated_by_TFs', width = 5.5, height = 6.5)
ggsave('TRAs_regulated_by_TFs.svg', width = 6, height = 5)


###Extended Data Fig.4h
DefaultAssay(filtered_TRA_TEC) <- 'symbol'
cor_plot_df <- FetchData(filtered_TRA_TEC[, Idents(filtered_TRA_TEC) %in% c('TEC(myo)')], 
                     vars = c(intersect(plot_df$TF[plot_df$Freq > 0.2], gsub('\\(.*\\)', '', names(filtered_TRA_TEC_regulons))), 'skeletal_muscle'))
cor_plot_df <- cbind(cor_plot_df, 
                     FetchData(filtered_TRA_TEC[, Idents(filtered_TRA_TEC) %in% c('TEC(myo)')],
                               vars = c('batch', 'Age_group')))

saveRDS(cor_plot_df, file = 'TRAs_regulated_by_TFs_regression.rds')
cor_plot_df$batch <- cor_plot_df$Age_group <- NULL
cor_plot_df <- reshape2::melt(cor_plot_df, id.vars = c('skeletal_muscle'), variable.name = 'features', value.name = 'exp')



p <- ggplot(cor_plot_df, aes(skeletal_muscle, exp, color =features))+
geom_point(alpha = 0.5)+
geom_smooth(method = 'lm', se = T)+
stat_cor(method = "pearson",size=13.75/.pt, label.y.npc = 0.995, show.legend = F)+
labs(y='Expression', x = 'TRAs % (skeletal_muscle)', color = 'TFs')+
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
         axis.title = element_text(size = 13.75),
          plot.title = element_text(size = 13.75, hjust=0.5),
          axis.text = element_text(size = 13.75),
          legend.justification = 0.5,
        legend.title = element_text(size = 13.75),
          legend.text = element_text(size = 13.75, face = 'italic'),
          #legend.position = c(0.9, 0.2),
          panel.border = element_rect(color ='black', fill =NA),
          panel.background= element_rect(fill ='white', color =NA),
          
         )+
guides(color = guide_legend(override.aes = list(size=5)))+
labs(title = 'Pearson correlation')
options(repr.plot.width=5.5, repr.plot.height=5)
p
ggsave('TRAs_regulated_by_TFs_regression.svg', width = 5.5, height = 5)
#my_saveplot_fun(p, 'TRAs_regulated_by_TFs_regression', width = 5, height = 6.5)



filtered_TRA_TEC$Age_group <- factor(filtered_TRA_TEC$Age_group , levels = c('Prepuberal', 'Adult', 'Aged'))

plot_df <- cbind(FetchData(filtered_TRA_TEC, vars = c('AIRE', 'FEZF2','Age_group')), 
                    identity = Idents(filtered_TRA_TEC))
plot_df <- plot_df[plot_df$identity == 'mTEC_hi', ]

###Fig.4d
plot_df <- cbind(FetchData(filtered_TRA_TEC, vars = c('AIRE', 'FEZF2','Age_group', 'batch')), 
                    identity = Idents(filtered_TRA_TEC))
plot_df <- plot_df[plot_df$identity == 'mTEC_hi', ]
plot_df <- melt(plot_df, id.vars = c('identity', 'Age_group', 'batch'),
                variable.name = 'features', value.name = 'Exp')
# plot_df <- plot_df[plot_df$Age_group %in% c('Prepuberal', 'Adult'),]
plot_df$Age_group <- factor(plot_df$Age_group , levels = c('Prepuberal', 'Adult', 'Aged'))
stat_test <- lapply(unique(plot_df$features), function(feature){
    tmp_data <- plot_df[plot_df$features %in% feature, ]
    ###only one mTEC_hi cell in Aged group
    #tmp_data <- tmp_data[!tmp_data$Age_group %in% c('Aged'), ]
    df <- rstatix::wilcox_test(
         Exp ~ Age_group, data = tmp_data
    )
    df$features <- feature
    
    max_val <- max(tmp_data$Exp)
    by_val <- 0.1* max_val
    df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = nrow(df))
    df
})

stat_test <- do.call(rbind, stat_test)
stat_test$p.adj <- p.adjust(stat_test$p, method = 'fdr')
stat_test$p.adj.signif <- 'ns'
stat_test$p.adj.signif[stat_test$p.adj < 0.05] <- '*'
stat_test$p.adj.signif[stat_test$p.adj < 0.01] <- '**'
stat_test$p.adj.signif[stat_test$p.adj < 0.001] <- '***'
stat_test$p.adj.signif[stat_test$p.adj < 0.0001] <- '****'

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
stat_test$p.adj[logi]

options(repr.plot.width=3.5, repr.plot.height=3.5)
p <- ggplot(plot_df, aes(Age_group, Exp))+
geom_violin(mapping = aes(fill = Age_group))+
# geom_boxplot(fill = 'white', width = 0.1, outlier.colour = NA)+
stat_summary(aes(y = Exp, group = Age_group),
             fun.y = mean, 
               fun.ymin = mean, 
               fun.ymax = mean,  geom="crossbar",color = 'red',
             size = 0.1, alpha = 1, width = 0.25*2)+
  stat_summary(fun.data = 'get_sem', geom = "errorbar", colour = "red",size = 0.25,
               width = 0.15*2,position = position_dodge( .9))+
geom_point(position = 'jitter', size = 0.1, alpha = 0.2)+

stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75/.pt)+

theme(strip.text.x = element_text(size = 13.75, face = 'italic'),
      strip.background.x = element_blank(),
        panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 13.75),
        axis.title.y = element_text(size = 13.75),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45,hjust = 1,size = 13.75),
        legend.text = element_text(size = 13.75),
        plot.title = element_text(size = 13.75, hjust=0.5, face = 'italic'),
        legend.title =element_text(size = 13.75)
        )+
scale_fill_manual(values = age_colors)+
guides(fill =F)+
scale_y_continuous(expand = c(0,0.29))+
labs(y = 'Expression Level')+
scale_x_discrete(limits = c('Prepuberal', 'Adult', 'Aged'), label = age_group_lbls_func)+
facet_wrap(.~features, scales = 'free_y')
p
ggsave('mTEC_hi_TRA.svg', width = 3.5, height = 3.5)

saveRDS(plot_df, file = 'mTEC_hi_TRA.rds')

###Fig.4e
plot_df <- cbind(FetchData(filtered_TRA_TEC, vars = c('AIRE', 'FEZF2','Age_group', 'batch')), 
                    identity = Idents(filtered_TRA_TEC))
plot_df <- plot_df[plot_df$identity == 'post_AIRE_mTEC', ]

plot_df <- melt(plot_df, id.vars = c('identity', 'Age_group', 'batch'),variable.name = 'features', value.name = 'Exp')
# plot_df <- plot_df[plot_df$Age_group %in% c('Prepuberal', 'Adult'),]
plot_df$Age_group <- factor(plot_df$Age_group , levels = c('Prepuberal', 'Adult', 'Aged'))
stat_test <- lapply(unique(plot_df$features), function(feature){
    tmp_data <- plot_df[plot_df$features %in% feature, ]
    ###only one mTEC_hi cell in Aged group
    #tmp_data <- tmp_data[!tmp_data$Age_group %in% c('Aged'), ]
    df <- rstatix::wilcox_test(
         Exp ~ Age_group, data = tmp_data
    )
    df$features <- feature
    
    max_val <- max(tmp_data$Exp)
    by_val <- 0.1* max_val
    df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = nrow(df))
    df
})

stat_test <- do.call(rbind, stat_test)
stat_test$p.adj <- p.adjust(stat_test$p, method = 'fdr')
stat_test$p.adj.signif <- 'ns'
stat_test$p.adj.signif[stat_test$p.adj < 0.05] <- '*'
stat_test$p.adj.signif[stat_test$p.adj < 0.01] <- '**'
stat_test$p.adj.signif[stat_test$p.adj < 0.001] <- '***'
stat_test$p.adj.signif[stat_test$p.adj < 0.0001] <- '****'

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
stat_test$p.adj[logi]

options(repr.plot.width=3.5, repr.plot.height=3.5)
p <- ggplot(plot_df, aes(Age_group, Exp))+
geom_violin(mapping = aes(fill = Age_group))+
# geom_boxplot(fill = 'white', width = 0.1, outlier.colour = NA)+
stat_summary(aes(y = Exp, group = Age_group),
             fun.y = mean, 
               fun.ymin = mean, 
               fun.ymax = mean,  geom="crossbar",color = 'red',
             size = 0.1, alpha = 1, width = 0.25*2)+
  stat_summary(fun.data = 'get_sem', geom = "errorbar", colour = "red",size = 0.25,
               width = 0.15*2,position = position_dodge( .9))+
geom_point(position = 'jitter', size = 0.1, alpha = 0.2)+

stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75/.pt)+

theme(strip.text.x = element_text(size = 13.75, face = 'italic'),
      strip.background.x = element_blank(),
        panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 13.75),
        axis.title.y = element_text(size = 13.75),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45,hjust = 1,size = 13.75),
        legend.text = element_text(size = 13.75),
        plot.title = element_text(size = 13.75, hjust=0.5, face = 'italic'),
        legend.title =element_text(size = 13.75)
        )+
scale_fill_manual(values = age_colors)+
guides(fill =F)+
scale_y_continuous(expand = c(0,0.29))+
labs(y = 'Expression Level')+
scale_x_discrete(limits = c('Prepuberal', 'Adult', 'Aged'), label = age_group_lbls_func)+
facet_wrap(.~features, scales = 'free_y')
p
ggsave('post_AIRE_mTEC_TRA.svg', width = 3.5, height = 3.5)
saveRDS(plot_df, file = 'post_AIRE_mTEC_TRA.rds')



###Extended Data Fig.4i
plot_df <- cbind(FetchData(filtered_TRA_TEC, vars = c('MYF6', 'MEF2C', 'MYOD1','TEAD4','Age_group', 'batch')), 
                    identity = Idents(filtered_TRA_TEC))
plot_df <- plot_df[plot_df$identity == 'TEC(myo)', ]

plot_df <- melt(plot_df, id.vars = c('identity', 'Age_group', 'batch'),variable.name = 'features', value.name = 'Exp')
# plot_df <- plot_df[plot_df$Age_group %in% c('Prepuberal', 'Adult'),]
plot_df$Age_group <- factor(plot_df$Age_group , levels = c('Prepuberal', 'Adult', 'Aged'))
stat_test <- lapply(unique(plot_df$features), function(feature){
    tmp_data <- plot_df[plot_df$features %in% feature, ]
    #tmp_data <- tmp_data[!tmp_data$Age_group %in% c('Aged'), ]
    df <- rstatix::wilcox_test(
         Exp ~ Age_group, data = tmp_data
    )
    df$features <- feature
    
    max_val <- max(tmp_data$Exp)
    by_val <- 0.08* max_val
    df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = nrow(df))
    df
})

stat_test <- do.call(rbind, stat_test)
stat_test$p.adj <- p.adjust(stat_test$p, method = 'fdr')
stat_test$p.adj.signif <- 'ns'
stat_test$p.adj.signif[stat_test$p.adj < 0.05] <- '*'
stat_test$p.adj.signif[stat_test$p.adj < 0.01] <- '**'
stat_test$p.adj.signif[stat_test$p.adj < 0.001] <- '***'
stat_test$p.adj.signif[stat_test$p.adj < 0.0001] <- '****'

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 2)


options(repr.plot.width=6, repr.plot.height=5)
p <- ggplot(plot_df, aes(Age_group, Exp))+
geom_violin(mapping = aes(fill = Age_group))+
# geom_boxplot(fill = 'white', width = 0.1, outlier.colour = NA)+
stat_summary(aes(y = Exp, group = Age_group),
             fun.y = mean, 
               fun.ymin = mean, 
               fun.ymax = mean,  geom="crossbar",color = 'red',
             size = 0.1, alpha = 1, width = 0.25*2)+
  stat_summary(fun.data = 'get_sem', geom = "errorbar", colour = "red",size = 0.25,
               width = 0.15*2,position = position_dodge( .9))+
geom_point(position = 'jitter', size = 0.1, alpha = 0.2)+

stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75/.pt, tip.length = 0.01)+

theme(strip.text.x = element_text(size = 13.75, face = 'italic'),
      strip.background.x = element_blank(),
        panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 13.75),
        axis.title.y = element_text(size = 13.75),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45,hjust = 1,size = 13.75),
        legend.text = element_text(size = 13.75),
        plot.title = element_text(size = 13.75, hjust=0.5, face = 'italic'),
        legend.title =element_text(size = 13.75)
        )+
scale_fill_manual(values = age_colors)+
guides(fill =F)+
scale_y_continuous(expand = c(0,0.29))+
labs(y = 'Expression Level')+
scale_x_discrete(limits = c('Prepuberal', 'Adult', 'Aged'), label = age_group_lbls_func)+
facet_wrap(.~features, scales = 'free_y', ncol = 4)
p
ggsave('TEC_myo_TRA.svg', width =11.5, height = 4)
saveRDS(plot_df, file = 'TEC_myo_TRA.rds')

###Fig.4f
plot_df <- cbind(FetchData(filtered_TRA_TEC, vars = c('TTN', 'RYR1','Age_group', 'batch')), 
                    identity = Idents(filtered_TRA_TEC))
plot_df <- plot_df[plot_df$identity == 'TEC(myo)', ]

plot_df <- melt(plot_df, id.vars = c('identity', 'Age_group', 'batch'),variable.name = 'features', value.name = 'Exp')
# plot_df <- plot_df[plot_df$Age_group %in% c('Prepuberal', 'Adult'),]
plot_df$Age_group <- factor(plot_df$Age_group , levels = c('Prepuberal', 'Adult', 'Aged'))
stat_test <- lapply(unique(plot_df$features), function(feature){
    tmp_data <- plot_df[plot_df$features %in% feature, ]
    #tmp_data <- tmp_data[!tmp_data$Age_group %in% c('Aged'), ]
    df <- rstatix::wilcox_test(
         Exp ~ Age_group, data = tmp_data
    )
    df$features <- feature
    
    max_val <- max(tmp_data$Exp)
    by_val <- 0.08* max_val
    df$y.position <- seq(from = max_val + by_val, by = by_val , length.out = nrow(df))
    df
})

stat_test <- do.call(rbind, stat_test)
stat_test$p.adj <- p.adjust(stat_test$p, method = 'fdr')
stat_test$p.adj.signif <- 'ns'
stat_test$p.adj.signif[stat_test$p.adj < 0.05] <- '*'
stat_test$p.adj.signif[stat_test$p.adj < 0.01] <- '**'
stat_test$p.adj.signif[stat_test$p.adj < 0.001] <- '***'
stat_test$p.adj.signif[stat_test$p.adj < 0.0001] <- '****'

logi <- stat_test$p.adj > 0.05 & stat_test$p.adj < 0.1
stat_test$p.adj.signif[logi] <- 
round(stat_test$p.adj[logi], 3)


options(repr.plot.width=3.5, repr.plot.height=3.5)
p <- ggplot(plot_df, aes(Age_group, Exp))+
geom_violin(mapping = aes(fill = Age_group))+
# geom_boxplot(fill = 'white', width = 0.1, outlier.colour = NA)+
stat_summary(aes(y = Exp, group = Age_group),
             fun.y = mean, 
               fun.ymin = mean, 
               fun.ymax = mean,  geom="crossbar",color = 'red',
             size = 0.1, alpha = 1, width = 0.25*2)+
  stat_summary(fun.data = 'get_sem', geom = "errorbar", colour = "red",size = 0.25,
               width = 0.15*2,position = position_dodge( .9))+
geom_point(position = 'jitter', size = 0.1, alpha = 0.2)+

stat_pvalue_manual(data = stat_test, label = 'p.adj.signif',  size = 13.75/.pt, tip.length = 0.01)+

theme(strip.text.x = element_text(size = 13.75, face = 'italic'),
      strip.background.x = element_blank(),
        panel.background = element_rect(color = NA, fill = 'white'), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 13.75),
        axis.title.y = element_text(size = 13.75),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45,hjust = 1,size = 13.75),
        legend.text = element_text(size = 13.75),
        plot.title = element_text(size = 13.75, hjust=0.5, face = 'italic'),
        legend.title =element_text(size = 13.75)
        )+
scale_fill_manual(values = age_colors)+
guides(fill =F)+
scale_y_continuous(expand = c(0,0.29))+
labs(y = 'Expression Level')+
scale_x_discrete(limits = c('Prepuberal', 'Adult', 'Aged'), label = age_group_lbls_func)+
facet_wrap(.~features, scales = 'free_y')
p
ggsave('TEC_myo_MG.svg', width = 3.5, height = 3.5)
saveRDS(plot_df, file = 'TEC_myo_MG.rds')







crosstalk <- readRDS('/data1/02.private/dengyj/analysis/thymus/crosstalk/250222//crosstalk.rds')

CD45_neg_idents_lvls <-  CD45_neg_idents <- c(
    'cTEC_hi', 'cTEC_lo', 'MKI67+cTEC', 'Immature_TEC', 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC', 
    'MKI67+mTEC',
'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated', 'Fb_1', 'Fb_2', 'MKI67+Fb', 'VSMCs', 'MKI67+VSMCs',
                                            'Endo', 'MKI67+Endo','Myelin', 'Mesothelial')

TEC_idents <- c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC', 'Immature_TEC', 'mTEC_lo', 'mTEC_hi', 
       'post_AIRE_mTEC', 
    'MKI67+mTEC',
'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated')

Fb_idents <- c('Fb_1', 'Fb_2', 'MKI67+Fb')

ETP_idents_lvls <- ETP_idents <- c('ETP', 'TP1', 'TP2')

con_result <- get_sig_LR_result(crosstalk, slot_used = 'condition', signal = 'intercellular',inter_p_threshold = 0.01)


con_result <- con_result[!con_result$ligand %in% names(crosstalk@complex_info) & 
                         !con_result$receptor %in% names(crosstalk@complex_info), ]

tbl <- table(crosstalk@identity,crosstalk@condition)#


clusters_filtered <- rownames(tbl)[apply(tbl,1, function(x){any(x<50)})]


ETP_DN_Rec_con_result <- con_result[con_result$ligand_cluster %in% 
                                    setdiff(CD45_neg_idents, clusters_filtered) & 
                                    con_result$receptor_cluster %in% setdiff(ETP_idents, clusters_filtered), ]

ETP_DN_Rec_con_result_intercellular <- get_plot_df(object = crosstalk, slot_used = 'condition',clusters_pairs = 
                                             unique(ETP_DN_Rec_con_result$clu_pairs),
                                do_scale = T,do_log = F,
                              LR_pairs = unique(ETP_DN_Rec_con_result$interaction_name), 
                             data_to_plot  = 'intercellular', value_max_cutoff =0.8, value_min_cutoff = 0.2)

ETP_DN_Rec_con_result_pnull <- get_plot_df(object = crosstalk, slot_used = 'condition',clusters_pairs = 
                                             unique(ETP_DN_Rec_con_result$clu_pairs),
                                do_scale = F,do_log = T,
                              LR_pairs = unique(ETP_DN_Rec_con_result$interaction_name), 
                             data_to_plot  = 'pnull')

ETP_DN_Rec_con_result_pnull_filtered <- ETP_DN_Rec_con_result_pnull[ETP_DN_Rec_con_result_pnull$greater_pnull >=2, ]

ETP_DN_Rec_con_result_plot_df <- merge(ETP_DN_Rec_con_result_intercellular, ETP_DN_Rec_con_result_pnull)

ETP_DN_Rec_con_result_plot_df$clusters_pairs <- 
factor(ETP_DN_Rec_con_result_plot_df$clusters_pairs, 
      levels = intersect(unlist(lapply(CD45_neg_idents_lvls, function(x){paste0(x, '~', ETP_idents_lvls)})), 
                         unique(ETP_DN_Rec_con_result_plot_df$clusters_pairs)))

order_mat <- function(mat){
    max_col <- apply(mat, 1, which.max)
    mat_scale <- t(apply(mat, 1, scale))
    order_mat <- lapply(1:ncol(mat), function(col) {
        mat_iter <- mat[max_col == col, , drop = F]
        mat_scale_iter <- mat_scale[max_col == col, , drop = F]
        if (nrow(mat_iter) > 0) {
            if (col == 1) {
                mat_iter <- mat_iter[order(mat_scale_iter[, col + 
                  1], decreasing = F), , drop = F]
            }
            else if (col == ncol(mat)) {
                mat_iter <- mat_iter[order(mat_scale_iter[, col - 
                  1], decreasing = T), , drop = F]
            }
            else {
                mat_iter_diff <- mat_iter[, col - 
                  1] - mat_iter[, col + 1]
                mat_iter <- mat_iter[order(mat_iter_diff, 
                  decreasing = T), , drop = F]
            }
            return(mat_iter)
        }
    })
    order_mat <- do.call(rbind, order_mat) 
    return(order_mat)
}

###Fig.4g
plot_df <- ETP_DN_Rec_con_result_plot_df[
    ETP_DN_Rec_con_result_plot_df$ligand_cluster %in% ETP_DN_Rec_con_result_pnull_filtered$ligand_cluster&
    ETP_DN_Rec_con_result_plot_df$ligand_cluster %in% TEC_idents & 
    ETP_DN_Rec_con_result_plot_df$interaction_name %in% ETP_DN_Rec_con_result_pnull_filtered$interaction_name&
    ETP_DN_Rec_con_result_plot_df$greater_pnull >=2, ]

plot_df <- plot_df[plot_df$interaction_name %in% c('CCL25~CCR9','CCL19~CCR10',
                                                   'DLL4~NOTCH1', 'IL7~IL7R', 'JAG1~NOTCH1',
                                                   'DLL4~NOTCH3',
                                                   # 'IL15~IL2RG',
                          'CXCL12~CXCR4', 'JAG1~NOTCH1',  'CCL19~CCR10', 'IL15~IL2RG', 
                          'FN1~CD44', 'ICAM1~SPN', 'APOE~LRP5', 'VCAN~ITGB1'), ]
stat_plot_df <- table(plot_df$interaction_name)
stat_plot_df <- stat_plot_df[order(stat_plot_df, decreasing = T)]
#stat_plot_df <- stat_plot_df[1:100]
plot_df <- plot_df[plot_df$interaction_name %in% names(stat_plot_df), ]
# plot_df$condition <- as.character(plot_df$condition)
plot_df$condition[plot_df$condition == 'Prepuberal'] <- 'Group I'
plot_df$condition[plot_df$condition == 'Adult'] <- 'Group II'
plot_df$condition[plot_df$condition == 'Aged'] <- 'Group III'
plot_df$condition <- factor(plot_df$condition, levels = c('Group I', 'Group II', 'Group III'))

p <- ggplot(plot_df, 
            aes(clusters_pairs,interaction_name, fill = intercellular, size = greater_pnull))+
geom_point(shape=21, stroke=0.1)+

scale_fill_gradientn(colours =c('lightgrey', 'red'))+

scale_y_discrete(limits=
                rev(
                    rownames(order_mat(apply(
                        table(plot_df$interaction_name, plot_df$condition), 2, function(x){x}))
                            )), position = 'right')+
  theme(panel.border = element_rect(color ='black', fill =NA),strip.background = element_blank(),
        panel.grid.major=element_line(color ='lightgrey'),
        panel.background = element_rect(color =NA, fill ='white'),
        strip.text.x = element_text(size = 13.75, family='sans'),
         legend.text = element_text(size = 13.75, family='sans'),
        legend.title = element_text(size = 13.75, family='sans'),
      #axis.text.y = element_text(size = 6, face = 'bold'), 
        axis.text.y = element_text(size=13.5),#axis.ticks = element_blank(),
        #axis.text.x = element_text(size = 3, angle =90, hjust = 1, face ='bold'),
        axis.text.x = element_blank(),legend.direction = 'vertical',
       )+
  labs(fill = 'Interaction', x = NULL, y= NULL)+
scale_size_continuous(range = c(1,3))+
guides(size = guide_legend(title = "-log10 P", override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Intensity"))+
facet_wrap(.~condition)

label_plot_df <- data.frame(clusters_pairs = levels(droplevels(plot_df$clusters_pairs)))
label_plot_df$ligand_cluster <- sapply(strsplit(label_plot_df$clusters_pairs, split = '~'), function(x){x[1]})
label_plot_df$receptor_cluster <- sapply(strsplit(label_plot_df$clusters_pairs, split = '~'), function(x){x[2]})
label_plot_df$clusters_pairs <- factor(label_plot_df$clusters_pairs, levels = 
                                       levels(droplevels(plot_df$clusters_pairs)))
label_plot_df <- lapply(intersect(c('Group I', 'Group II', 'Group III'), plot_df$condition), function(age){
    new_df <- label_plot_df
    new_df$age <- age
    new_df
})
label_plot_df <- do.call(rbind, label_plot_df)
label_plot_df$ligand_cluster<- factor(label_plot_df$ligand_cluster, levels = unique(label_plot_df$ligand_cluster))

label_plot_df$receptor_cluster<- factor(label_plot_df$receptor_cluster, levels = unique(label_plot_df$receptor_cluster))

options(repr.plot.width=15, repr.plot.height=0.5)
label_plot_ligand <- ggplot(label_plot_df)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = ligand_cluster))+
scale_fill_manual(values =cluster_colors, label = lbls_func)+
guides(fill = guide_legend(ncol  =2))+
#scale_y_discrete(limits = rev(c('Hema', 'Epi', 'Fibro', 'Endo')))+
theme(legend.text = element_text(size = 13.75, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())+
facet_wrap(.~age)



label_plot_receptor <- ggplot(label_plot_df)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = receptor_cluster))+
scale_fill_manual(values =cluster_colors)+
#guides(fill = F)+
#scale_y_discrete(limits = rev(c('Hema', 'Epi', 'Fibro', 'Endo')))+
theme(
    legend.text = element_text(size = 13.75, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())+
facet_wrap(.~age)


options(repr.plot.width=9, repr.plot.height=4.5)
p <- wrap_plots(list(p, label_plot_ligand, label_plot_receptor))+plot_layout(ncol = 1,guides = 'collect',
                                                                        heights = c(19,1,1))&
theme(plot.margin = margin(0,0,0,1,'lines'), legend.position = 'bottom')
p


#my_saveplot_fun(p, 'TEC_crosstalk', width = 11,height = 5.5)
ggsave('TEC_crosstalk.svg', width = 9,height = 4.5)


###Fig.4h
plot_df <- ETP_DN_Rec_con_result_plot_df[
    ETP_DN_Rec_con_result_plot_df$ligand_cluster %in% ETP_DN_Rec_con_result_pnull_filtered$ligand_cluster&
    ETP_DN_Rec_con_result_plot_df$ligand_cluster %in% Fb_idents & 
    ETP_DN_Rec_con_result_plot_df$interaction_name %in% ETP_DN_Rec_con_result_pnull_filtered$interaction_name&
    ETP_DN_Rec_con_result_plot_df$greater_pnull >=2, ]

plot_df <- plot_df[plot_df$interaction_name %in% c('FN1~ITGB1', #'CD58~CD2', 
                                                   'VCAM1~ITGB2', 'COL1A1~CD93', 
                                                   'VCAM1~ITGA4',
                          'JAG1~NOTCH1', 'JAG1~NOTCH2', 'TGFB2~TGFBR2', 'COL1A2~ITGA2', 
                          'IGF2~IGF1R', 'ICAM1~ITGAL', 'CXCL12~CXCR4'), ]
stat_plot_df <- table(plot_df$interaction_name)
stat_plot_df <- stat_plot_df[order(stat_plot_df, decreasing = T)]
#stat_plot_df <- stat_plot_df[1:100]
plot_df <- plot_df[plot_df$interaction_name %in% names(stat_plot_df), ]

plot_df$condition[plot_df$condition == 'Prepuberal'] <- 'Group I'
plot_df$condition[plot_df$condition == 'Adult'] <- 'Group II'
plot_df$condition[plot_df$condition == 'Aged'] <- 'Group III'
plot_df$condition <- factor(plot_df$condition, levels = c('Group I', 'Group II', 'Group III'))

p <- ggplot(plot_df, 
            aes(clusters_pairs,interaction_name, fill = intercellular, size = greater_pnull))+
geom_point(shape=21, stroke=0.1)+

scale_fill_gradientn(colours =c('lightgrey', 'red'))+
#scale_y_discrete(limits = rev(names(stat_plot_df)))+
scale_y_discrete(limits=
                rev(
                    rownames(order_mat(apply(
                        table(plot_df$interaction_name, plot_df$condition), 2, function(x){x}))
                            )), position = 'right')+
  theme(panel.border = element_rect(color ='black', fill =NA),strip.background = element_blank(),
        panel.grid.major=element_line(color ='lightgrey'),
        panel.background = element_rect(color =NA, fill ='white'),
        strip.text.x = element_text(size = 13.75, family='sans'),
         legend.text = element_text(size = 13.75, family='sans'),
        legend.title = element_text(size = 13.75, family='sans'),
      #axis.text.y = element_text(size = 6, face = 'bold'), 
        axis.text.y = element_text(size=13.75),#axis.ticks = element_blank(),
        #axis.text.x = element_text(size = 3, angle =90, hjust = 1, face ='bold'),
        axis.text.x = element_blank(),legend.direction = 'vertical',
       )+
  labs(fill = 'Interaction', x = NULL, y= NULL)+
scale_size_continuous(range = c(2,3))+
guides(size = guide_legend(title = "-log10 P", order = 1,override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Intensity", order = 2))+
facet_wrap(.~condition)

label_plot_df <- data.frame(clusters_pairs = levels(droplevels(plot_df$clusters_pairs)))
label_plot_df$ligand_cluster <- sapply(strsplit(label_plot_df$clusters_pairs, split = '~'), function(x){x[1]})
label_plot_df$receptor_cluster <- sapply(strsplit(label_plot_df$clusters_pairs, split = '~'), function(x){x[2]})
label_plot_df$clusters_pairs <- factor(label_plot_df$clusters_pairs, levels = 
                                       levels(droplevels(plot_df$clusters_pairs)))
label_plot_df <- lapply(intersect(c('Group I', 'Group II', 'Group III'), plot_df$condition), function(age){
    new_df <- label_plot_df
    new_df$age <- age
    new_df
})
label_plot_df <- do.call(rbind, label_plot_df)
label_plot_df$ligand_cluster<- factor(label_plot_df$ligand_cluster, levels = unique(label_plot_df$ligand_cluster))

label_plot_df$receptor_cluster<- factor(label_plot_df$receptor_cluster, levels = unique(label_plot_df$receptor_cluster))

options(repr.plot.width=15, repr.plot.height=0.5)
label_plot_ligand <- ggplot(label_plot_df)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = ligand_cluster))+
scale_fill_manual(values =cluster_colors, label = lbls_func)+
guides(fill = guide_legend(ncol  =1))+
#scale_y_discrete(limits = rev(c('Hema', 'Epi', 'Fibro', 'Endo')))+
theme(legend.text = element_text(size = 13.75, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())+
facet_wrap(.~age)



label_plot_receptor <- ggplot(label_plot_df)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = receptor_cluster))+
scale_fill_manual(values =cluster_colors)+

theme(legend.text = element_text(size = 13.75, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 13.75, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())+
facet_wrap(.~age)


options(repr.plot.width=9, repr.plot.height=4.5)
p <- wrap_plots(list(p, label_plot_ligand, label_plot_receptor))+plot_layout(ncol = 1,guides = 'collect',
                                                                        heights = c(19,1,1))&
theme(plot.margin = margin(0,0,0,1,'lines'), legend.position = 'bottom')
p


ggsave('Fb_crosstalk.svg', width = 9,height = 4.5)












