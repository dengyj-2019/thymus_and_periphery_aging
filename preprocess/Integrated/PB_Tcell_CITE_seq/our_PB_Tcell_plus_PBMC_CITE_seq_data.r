library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)

library(harmony)


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')

###public CITE-seq data of PBMC(GSE164378)
PBMC <- readRDS('/data1/02.private/dengyj/analysis/public_data/GSE164378/PBMC_0d.rds')

Idents(PBMC) <- PBMC$celltype.l2

####our PB Tcell data
Tcell <-readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell.rds')

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(PBMC, reduction = 'tsne',group.by = 'celltype.l2',label = T)

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(Tcell, label = T)

PBMC$project <- 'a'
Tcell$project <- 'b'

PBMC$PBMC_idents <- as.character(Idents(PBMC))
Tcell$Tcell_Idents <- as.character(Idents(Tcell))


cluster_colors <- c('CD8+TEM-B'='#cb50a0',
                    'CD8+TEM-K'='#db72b4',#
                    'CD8+RTE'='#276d6d','CD4+TN'='#c1d4eb','CD4+CTL'='#7c75ce',#
                    'CD27lo_Vδ1T'='#93c88e','CD8+TN'='#DDC4DA',
                    'IFN_response'='#00E8F7','CD4+TEM'='#909ee7',##'
                    'Treg'='#6e7cad','CD4+TCM'='#9cb9e5',
                    'MAIT'='#b1b772',
                    'CD4+RTE'='#52a298','Vδ2Vγ9T'='#A0D7C9',
                    'CD27+Vδ1T'='#9cb39b',#
                    'IFN_T'='#ec807c',
                    'T_un'='#747474','Tex'='#A1B4BB', 'CD4+TEFF' = '#55b3db',##
                    'CD8+TCM'='#E59CC4')

combined_PB <- merge(PBMC[, PBMC$celltype.l1 %in% c('CD4 T', 'CD8 T', 'other T')], ###Tcell
                    Tcell)

combined_PB <- NormalizeData(combined_PB)
combined_PB <- FindVariableFeatures(combined_PB, selection.method = "vst", nfeatures = 2000)
combined_PB <- ScaleData(combined_PB)
combined_PB <- RunPCA(combined_PB, npcs = 50)

combined_PB <- RunHarmony(combined_PB, 'batch')

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combined_PB, reduction = 'harmony',ndims = 50)

combined_PB <- RunUMAP(combined_PB, reduction = "harmony", dims = 1:30)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(combined_PB, reduction = 'umap',label = T, group.by = 'project')
DimPlot(combined_PB, reduction = 'umap',label = T, group.by = 'batch')

sub_T1 <- combined_PB[, combined_PB$project == 'a']

options(repr.plot.width=11, repr.plot.height=9)
DimPlot(sub_T1, reduction = 'umap', pt.size = 1, label = T)+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title =element_text(size = 15))+
guides(color = guide_legend(override.aes = list(size = 5)))

sub_T2 <- combined_PB[, combined_PB$project == 'b']

levels(sub_T2) <- c('CD4+RTE', 'CD4+TN', 'CD4+TCM', 'CD8+RTE','CD8+TN', 'CD8+TCM',
           'CD4+TEFF', 'Treg', 'CD4+TEM', 'CD8+TEM-K', 
           'CD4+CTL', 'CD8+TEM-B',
           'MAIT', 'Vδ2Vγ9T', 'CD27+Vδ1T',
           'CD27lo_Vδ1T', 'Tex', 'IFN_T','T_un')

p <- DimPlot(sub_T2,  label = F,pt.size = 0.5)+
scale_color_manual(values = cluster_colors)
new_label <- as.character(1:nlevels(sub_T2))
names(new_label) <- levels(sub_T2)

options(repr.plot.width=12.6, repr.plot.height=8)
p <- add_label(p, use_new_label = T, new_label = c(new_label),
          location_adj= list('4'=c(0.01, 0.035)),
         label_arg = list(color='black', size=6))+
  guides(color=guide_legend(override.aes = list(size=3.3, label = c(new_label)), ncol=2))+
theme(axis.title = element_blank(),axis.ticks = element_blank(), 
     axis.text = element_blank(), axis.line  = element_blank(), 
     panel.border = element_rect(color = 'black', fill = NA), legend.key.width = 
      unit(2, "lines"),
   legend.key.height = unit(2, "lines"),
      legend.text = element_text(size=20), legend.position = 'right', 
      legend.justification='center')
p

DefaultAssay(sub_T1) <- 'ADT'
sub_T1 <- NormalizeData(sub_T1,normalization.method = 'CLR')
DefaultAssay(sub_T1) <- 'RNA'

options(repr.plot.width=16, repr.plot.height=16)#Tcell
DefaultAssay(sub_T1) <- 'ADT'
plot_list <- lapply(c('CD45RA', 'CD45RO', 'CD4-2', 'CD8a', 'CD127', 'CD279', 'TIGIT', 'CD27',
              'CD28', 'CD274', 'CD278', 'CX3CR1', 'CD69','CD31', 'CD103','CD38-1'
), function(feature){
    p <- 
    FeaturePlot(sub_T1,reduction = 'umap',
            features = feature, 
           #cols = c('lightgrey', 'red'), 
                min.cutoff= 'q1', max.cutoff='q99', pt.size=0.25)+
    theme(panel.border = element_rect(color = 'black', fill = NA), 
     panel.background = element_rect(color = NA, fill = 'white'),
     axis.title = element_blank(),axis.ticks = element_blank(), 
     axis.text = element_blank(), axis.line  = element_blank(),
     plot.title = element_blank(),
     legend.position = c(0.55,0.775),
          legend.direction = 'horizontal',
    legend.key.size = unit(0.045, 'npc') , 
     legend.text = element_text(size=15, angle = 45, hjust= 1))+
    annotation_custom(
        grid::textGrob(label = feature,hjust = 0.5,vjust=0.5,
                       x=grid::unit(0.835,"npc") , 
                       y=grid::unit(0.925,"npc"), 
                       gp=gpar(col = 'black', fontsize = 20))) +
    scale_color_gradientn(colours = c('lightgrey', 'red'), 
                          label = function(x){sprintf(paste0('%0.', 1, 'f'), x)})
    return(p)
})

DefaultAssay(sub_T1) <- 'RNA'

options(repr.plot.width=16, repr.plot.height=16)
wrap_plots(plot_list)+plot_layout(ncol =4)

DefaultAssay(sub_T1) <- 'RNA'
plot_df <- 
cbind(FetchData(sub_T1, vars = c('CD45RA', 'CD45RO', 'CD4-2', 'CD8a', 'CD31','CD38-1','CD38-2', 'SOX4'), 
                    slot = 'data'), Embeddings(sub_T1, reduction = 'umap'))
colnames(plot_df) <- gsub('adt_', '', colnames(plot_df))


plot_list <- lapply(c('CD45RA', 'CD45RO', 'CD4-2', 'CD8a', 'CD31','CD38-1','CD38-2', 'SOX4'), function(feature){
    if(feature == 'SOX4'){
        colors <- c('lightgrey', 'red')
        fontface='italic'
        
    }else{
        colors <- c('#d1cdd5', '#2100ff')
        fontface='plain'
    }
    plot_df_iter <- plot_df[, c('UMAP_1', 'UMAP_2', feature)]
    min_cutoff <- quantile(plot_df_iter[, feature], 0.05)
    max_cutoff <- quantile(plot_df_iter[, feature], 0.95)
    plot_df_iter[, feature][plot_df_iter[, feature] <= min_cutoff] <- min_cutoff
    plot_df_iter[, feature][plot_df_iter[, feature] >= max_cutoff] <- max_cutoff
    p <- ggplot(plot_df_iter, aes_string('UMAP_1', 'UMAP_2', color =paste0('`', feature, '`')))+
    geom_point(size = 0.1)+ 

    theme(panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'), 
        #element_blank(),#element_text(size = 25, face = 'italic', family = 'sans'),
        legend.text=element_text(size =14, angle=45, hjust=0.5, family = 'sans'),
      axis.ticks=element_blank(),axis.title=element_blank(),#element_text(size = 15),
      axis.text=element_blank(),axis.line=element_blank(), legend.key.size = unit(0.03, 'npc'),
       legend.position = c(0.8, 0.82), legend.title=element_blank(),legend.background =element_blank(),
       legend.direction = 'horizontal')+
    #labs(title = feature)
    annotation_custom(grid::textGrob(label = feature,hjust = 0,
                                   x=grid::unit(0.74,"npc") ,                                
                                   y=grid::unit(0.95,"npc"), 
                                gp=grid::gpar(col = 'black', fontsize = 15, fontface=fontface)))

    p <- p+ scale_color_gradientn(colours = colors, 
                                    label = function(x){sprintf(paste0('%0.', 1, 'f'), x)})
    return(p)
})

options(repr.plot.width=16, repr.plot.height=8)
wrap_plots(plot_list)+plot_layout(ncol =4)







saveRDS(combined_PB, file = 'combined_PB.rds')
