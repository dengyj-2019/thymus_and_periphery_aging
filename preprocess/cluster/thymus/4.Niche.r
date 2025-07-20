library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(openxlsx)
library(ggpubr)
library(harmony)
library(SCopeLoomR)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
# gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
#                          columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
load('/data1/02.private/dengyj/analysis/database/specie_conversion/specie_conversion.rds')
# load('/data1/02.private/dengyj/analysis/thymus/TRAs/TRAs_df.rds')
library(Rmisc)

# load('/data1/02.private/dengyj/analysis/thymus/TRAs/AIRE_genes/AIRE_genes.rds')
# load('/data1/02.private/dengyj/analysis/thymus/TRAs/AIRE_genes/APS_1.rds')

library(patchwork)

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle.R')

source('/data1/02.private/dengyj/analysis/mycode//plot_helper.r')



cluster_colors <- c('Fb_1' ='#E5D2DD','cTEC_hi' ='#70a887', 'MKI67+Fb'='#E59CC4', 
                   'Endo'='#F3B1A0', 'Fb_2'='#D6E7A3', 'TEC(neuron)'='#57C3F3', 'mTEC_lo'='#5a829d',#'#476D87',
'VSMCs'='#E95C59', 'mTEC_hi'='#F1BB72','Ionocyte'='#cdb979',
                    'MKI67+mTEC'='#AB3282', 'TEC(myo)'='#7db0bb',#'#91D0BE', 
                    'post_AIRE_mTEC'='#92a5a5', 'Tuft'='#c370ac', #'#3c97c4',
                    'Ionocyte'='#cdb979', 'Tuft/Ionocyte' = '#ff001c',
                    'Ciliated'='#cc5f91', 'cTEC_lo'='#c4c04f', 'MKI67+cTEC'='#ff7e00', 'Immature_TEC'='#9dabd5', 
                    'Mesothelial'='#d73e4b', 'TEC'='#ae80dc',#'#639791'
                        'MKI67+VSMCs'='#7ed670', 'MKI67+Endo'='#698c68', 'Myelin'='#ff8ae0'


                   )





color_used <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
                "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
                "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
                "#40E0D0","#5F9EA0","#FF1493",
                "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
                "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

Idents_proj <- function(
    obj1, 
    obj2
){
    new_identity <- rbind(data.frame(b = as.character(Idents(obj1)), row.names = names(Idents(obj1))),
                    data.frame(b = as.character(Idents(obj2)), 
                    row.names = names(Idents(obj2)))[setdiff(colnames(obj2), colnames(obj1)),,drop=F])
    obj2 <- AddMetaData(obj2, metadata = new_identity, col.name ='new_identity')            
    Idents(obj2) <- obj2$new_identity
    return(obj2)
}

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')
SM <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/SM/humanSM.rds')

combined_Niche <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/prepare_data_with_ensembl_id/combined_Niche_filter_mt.rds')

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
'#968175'
)

DefaultAssay(combined_Niche) <- 'symbol'

combined_Niche <- NormalizeData(combined_Niche)
combined_Niche <- FindVariableFeatures(combined_Niche, selection.method = "vst", nfeatures = 2000)
combined_Niche <- CellCycleScoring(combined_Niche, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
combined_Niche$CC.Difference <- combined_Niche$S.Score - combined_Niche$G2M.Score
# plan('multicore', workers = 16)
# options(future.globals.maxSize = 10 * 1024^3)
combined_Niche <- ScaleData(combined_Niche)#, vars.to.regress = c('S.Score', 'G2M.Score'))
# plan('sequential')
combined_Niche <- RunPCA(combined_Niche, npcs = 50)

combined_Niche <- RunHarmony(combined_Niche, 'batch', assay.use = 'symbol')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(combined_Niche, ndims = 50, reduction = 'harmony')

#combined_Niche <- RunTSNE(combined_Niche, reduction = "harmony", dims = 1:30)
combined_Niche <- RunUMAP(combined_Niche, reduction = "harmony", dims = 1:10, verbose = F)
combined_Niche <- FindNeighbors(combined_Niche, reduction = "harmony", dims = 1:10, verbose = F)

combined_Niche <- AddMetaData(combined_Niche, Idents(orig_combined_Niche), col.name = 'orig_idents')



options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(combined_Niche, group.by = 'Anno_level_5', label = T)+
scale_color_manual(values = c(my36colors))

options(repr.plot.width=10, repr.plot.height=7.5)
DimPlot(combined_Niche, group.by = 'orig_idents', label = T)+
scale_color_manual(values = c(my36colors))

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combined_Niche, group.by = 'Phase', label = T)



options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_Niche, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche, features = c('CD3E', 'CD3D', 'CD3G'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche, features = c('GATA1', 'KLF1', 'GYPA'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche, features = c('CD36', 'CLDN5', 'CDH5'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_Niche, features = c('COL1A1', 'COL1A2', 'PDGFRA', 'PDGFRB'))

options(repr.plot.width=16, repr.plot.height=16)#Endo
FeaturePlot(combined_Niche, features = c('CD36', 'CLDN5', 'CDH5', 'RGS5', 'PECAM1', 'SOX7', 'CAV1', 
                                        'SPP1', 'COL15A1'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_Niche, features = c('PCNA', 'MKI67', 'CDK1'))

options(repr.plot.width=16, repr.plot.height=16)####Mes特征，c5是双细胞
FeaturePlot(combined_Niche, features = c('FRZB', 'FBN1', 'COL1A1', 'COL3A1', 'COLEC11', 'MKI67', 
                                        'PDGFRA', 'SFRP2', 'PTN'))



options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_Niche, reduction = 'umap',
            features = c('KRT5', 'KRT8','PRSS16','EPCAM'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche, features = c('CD3E', 'CD3D', 'CD3G', 'PTPRC'))

options(repr.plot.width=12, repr.plot.height=12)#########
FeaturePlot(combined_Niche, features = c('PTPRC'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche, features = c('IGHA1', 'IGKC', 'JCHAIN', 'MZB1'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(combined_Niche, features = c('MPZ'))

plan('multicore', workers = 16)
options(future.globals.maxSize = 20 * 1024^3)
combined_Niche <- FindClusters(combined_Niche, resolution = seq(0.1,1, 0.05))
plan('sequential')

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('symbol_snn_res.', i)
    p <- DimPlot(combined_Niche, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))+
    scale_color_manual(values = color_used)
    print(p)
} 

# 
# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
# combined_Niche_markers <- FindAllMarkers(combined_Niche, pseudocount.use = 0, logfc.threshold = log(1.5))
# plan('sequential')

# saveRDS(combined_Niche_markers, file = 'combined_Niche_markers.rds')
# openxlsx::write.xlsx(combined_Niche_markers, file = 'combined_Niche_markers.xlsx')



combined_TEC <- combined_Niche[, combined_Niche$symbol_snn_res.0.1 %in% c('2','5','7')]
# combined_TEC <- NormalizeData(combined_TEC)
combined_TEC <- FindVariableFeatures(combined_TEC, selection.method = "vst", nfeatures = 2000)
combined_TEC <- ScaleData(combined_TEC)#, vars.to.regress = c('S.Score', 'G2M.Score'))
combined_TEC <- RunPCA(combined_TEC, npcs = 50)

combined_TEC <- RunHarmony(combined_TEC, 'batch', assay.use = 'symbol')

options(repr.plot.width=7.5, repr.plot.height=7.5)
ElbowPlot(combined_TEC, ndims = 50, reduction = 'harmony')

combined_TEC <- RunUMAP(combined_TEC, reduction = "harmony", dims = 1:5, verbose = F)



options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(combined_TEC, group.by = 'Phase', label = T)

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(combined_TEC, group.by = 'orig_idents', label = T)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, features = c('PDGFRA', 'COL1A1', 'COL3A1', 'PDGFRB', 'CLDN5', 
                                      'CDH5', 'CD36', 'RGS5', 'FBN1'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, features = c('CD3E', 'CD3D', 'CD3G'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_TEC, features = c('CCL19'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_TEC, features = c('AIRE'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_TEC, features = c('KRT1'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(combined_TEC, features = c('MPZ'))

options(repr.plot.width=8, repr.plot.height=8)####可能是
FeaturePlot(combined_TEC, features = c('PTH')) 

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, #reduction = 'tsne',
            features = c('CCL25', 'KRT5', 'KRT8', 'PSMB11', 'PRSS16', 
                                      'FOXN1', 'EPCAM', 'PAX1', 'PAX9'))

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(combined_TEC, #reduction = 'tsne',
            features = c('CCL19', 'CCL21'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, #reduction = 'tsne',
            features = c('DES', 'MYOG', 'MYLPF', 'MYF6'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, #reduction = 'tsne',
            features = c('KRT14', 'MYF6', 'BEX1', 'POU2F3', 'FOXI1', 'MSLN', 
                        'ATOH1'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_TEC, features = c('CD3E', 'CD3D', 'CD3G'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_TEC, features = c('CD36', 'CLDN5', 'CDH5'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC, features = c('COL1A1', 'COL1A2', 'PDGFRA', 'PDGFRB'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combined_TEC, features = c('CD36', 'CLDN5', 'CDH5', 'RGS5', 'PECAM1', 'SOX7', 'CAV1', 
                                        'SPP1', 'COL15A1'))

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combined_TEC, features = c('DLK2'))



options(repr.plot.width=10, repr.plot.height=8)
DimPlot(combined_TEC, group.by = 'orig_idents', label = T)+
scale_color_manual(values = c(my36colors))

combined_TEC <- FindNeighbors(combined_TEC, reduction = "harmony", dims = 1:5)

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
combined_TEC <- FindClusters(combined_TEC, resolution = seq(0.2, 1, 0.05), verbose = F)
plan('sequential')

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('symbol_snn_res.', i)
    p <- DimPlot(combined_TEC, reduction = 'umap',label = T, label.size = 20/.pt,
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

combined_TEC_S1 <- combined_TEC[, combined_TEC$symbol_snn_res.0.35 %in% c('5')]

combined_TEC_S1 <- RunUMAP(combined_TEC_S1, reduction = 'harmony', verbose = F, dims = 1:20)

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(combined_TEC_S1, group.by = 'orig_idents', label = T)+
scale_color_manual(values = c(my36colors))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC_S1, #reduction = 'tsne',
            features = c( 'POU2F3', 'FOXI1', 'MSLN'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC_S1, #reduction = 'tsne',
            features = c( 'PDGFRA', 'CLDN5', 'COL1A2', 'COL1A1'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_TEC_S1, #reduction = 'tsne',
            features = c( 'KRT5', 'KRT8', 'EPCAM'))

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(combined_TEC_S1, #reduction = 'tsne',
            features = c( 'AIRE', 'KRT1'))

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(combined_TEC_S1, #reduction = 'tsne',
            features = c( 'CCL19', 'CCL21'))

combined_TEC_S1 <- FindNeighbors(combined_TEC_S1, reduction = 'harmony', dims = 1:20, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multisession', workers = 10)
combined_TEC_S1 <- FindClusters(combined_TEC_S1, resolution = seq(0.2, 1.5, 0.05), verbose = F)
# plan('sequential')
###输出不同分辨率的聚类结果



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('symbol_snn_res.', i)
    p <- DimPlot(combined_TEC_S1, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

# options(future.globals.maxSize = 10 * 1024^3)
# plan('multicore', workers = 10)
Idents(combined_TEC_S1) <- combined_TEC_S1$`symbol_snn_res.0.5`
combined_TEC_S1_markers <- FindAllMarkers(combined_TEC_S1, pseudocount.use = 0, logfc.threshold = log(1.5))
# plan('sequential')
saveRDS(combined_TEC_S1_markers, file = 'combined_TEC_S1_markers.rds')
openxlsx::write.xlsx(combined_TEC_S1_markers, file = 'combined_TEC_S1_markers.xlsx')


options(repr.plot.width=16, repr.plot.height=16)
for(i in unique(combined_TEC_S1_markers$cluster)){
    genes <- combined_TEC_S1_markers[combined_TEC_S1_markers$cluster == i & 
                                   combined_TEC_S1_markers$avg_logFC > log(1.5) &
                          combined_TEC_S1_markers$p_val_adj < 0.05, 'gene'][1:16]
    genes <- genes[!is.na(genes)]
    if(length(genes) > 0){
        p <- FeaturePlot(combined_TEC_S1, reduction = 'umap', cols = c('lightgrey', 'red'),
                features = genes)
        print(p)
        print(paste0('cluster: ', i)) 
    }
}  

combined_TEC_S1_S1 <- combined_TEC_S1[, combined_TEC_S1$symbol_snn_res.0.5 %in% c('4')]

combined_TEC_S1_S1 <- FindNeighbors(combined_TEC_S1_S1, reduction = 'harmony', dims = 1:20, verbose = F)


# options(future.globals.maxSize = 10 * 1024^3)
# plan('multisession', workers = 10)
combined_TEC_S1_S1 <- FindClusters(combined_TEC_S1_S1, resolution = seq(0.2, 1, 0.05), verbose = F)
# plan('sequential')
###输出不同分辨率的聚类结果



options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 0.5, 0.05)){
    res <- paste0('symbol_snn_res.', i)
    p <- DimPlot(combined_TEC_S1_S1, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

combined_TEC_S1_S1$Identity <- as.character(combined_TEC_S1_S1$`symbol_snn_res.0.5`)
combined_TEC_S1_S1$Identity[combined_TEC_S1_S1$Identity %in% c('1')]<- 'Tuft'
combined_TEC_S1_S1$Identity[combined_TEC_S1_S1$Identity %in% c('0')]<- 'Ionocyte'
Idents(combined_TEC_S1_S1) <- combined_TEC_S1_S1$Identity
##亚群投射
##combined_TEC_S1_S1 <- Idents_proj(combined_TEC_S1_S1_sub, combined_TEC_S1_S1)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combined_TEC_S1_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


combined_TEC_S1$Identity <- as.character(combined_TEC_S1$`symbol_snn_res.0.5`)
combined_TEC_S1$Identity[combined_TEC_S1$Identity %in% c('2','5')]<- 'mTEC_lo'
combined_TEC_S1$Identity[combined_TEC_S1$Identity %in% c('0','1','6','7','8','9')]<- 'doublets'
combined_TEC_S1$Identity[combined_TEC_S1$Identity %in% c('3')]<- 'Mesothelial'
Idents(combined_TEC_S1) <- combined_TEC_S1$Identity
##亚群投射
combined_TEC_S1 <- Idents_proj(combined_TEC_S1_S1, combined_TEC_S1)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combined_TEC_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))


combined_TEC$Identity <- as.character(combined_TEC$`symbol_snn_res.0.35`)
combined_TEC$Identity[combined_TEC$Identity %in% c('0')]<- 'cTEC_hi'
combined_TEC$Identity[combined_TEC$Identity %in% c('3', '1')]<- 'cTEC_lo'
combined_TEC$Identity[combined_TEC$Identity %in% c('4')]<- 'Immature_TEC'
combined_TEC$Identity[combined_TEC$Identity %in% c('2')]<- 'mTEC_lo'
combined_TEC$Identity[combined_TEC$Identity %in% c('8')]<- 'mTEC_hi'
combined_TEC$Identity[combined_TEC$Identity %in% c('9')]<- 'TEC(myo)'
combined_TEC$Identity[combined_TEC$Identity %in% c('6')]<- 'TEC(neuron)'
combined_TEC$Identity[combined_TEC$Identity %in% c('7')]<- 'post_AIRE_mTEC'
combined_TEC$Identity[combined_TEC$Identity %in% c('10')]<- 'Ciliated'
Idents(combined_TEC) <- combined_TEC$Identity
##亚群投射
combined_TEC <- Idents_proj(combined_TEC_S1, combined_TEC)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combined_TEC, reduction = 'umap',label = T, label.size = 15/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 15))+
guides(color = guide_legend(override.aes = list(size = 7)))






combined_Niche_S1 <- combined_Niche[, combined_Niche$`symbol_snn_res.0.1` %in% c('6')]

combined_Niche_S1 <- FindVariableFeatures(combined_Niche_S1, selection.method = "vst", nfeatures = 2000, verbose = F)
combined_Niche_S1 <- ScaleData(combined_Niche_S1, verbose = F)#, vars.to.regress = c('S.Score', 'G2M.Score'))
combined_Niche_S1 <- RunPCA(combined_Niche_S1, npcs = 50, verbose = FALSE)

combined_Niche_S1 <- RunHarmony(combined_Niche_S1, 'batch', assay.use = 'symbol')

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combined_Niche_S1, ndims = 50, reduction = 'harmony')

combined_Niche_S1 <- RunUMAP(combined_Niche_S1, reduction = "harmony", dims = 1:10, verbose  = F)

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche_S1, features = c('GATA1', 'KLF1', 'GYPA'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche_S1, features = c('CD3E', 'CD3D', 'CD3G', 'PTPRC'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche_S1, features = c('IGKC', 'IGHA1', 'JCHAIN', 'MZB1'))

options(repr.plot.width=16, repr.plot.height=16)#########
FeaturePlot(combined_Niche_S1, features = c('CD36', 'CLDN5', 'CDH5'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combined_Niche_S1, features = c('COL1A1', 'COL1A2', 'PDGFRA', 'PDGFRB'))

options(repr.plot.width=16, repr.plot.height=16)#Endo
FeaturePlot(combined_Niche_S1, features = c('CD36', 'CLDN5', 'CDH5', 'RGS5', 'PECAM1', 'SOX7', 'CAV1', 
                                        'SPP1', 'COL15A1'))

options(repr.plot.width=16, repr.plot.height=16)#Endo
FeaturePlot(combined_Niche_S1, features = c('CCL25', 'KRT5', 'KRT8', 'PRSS16'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_Niche_S1, features = c('AIRE'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_Niche_S1, features = c('MPZ'))

options(repr.plot.width=8, repr.plot.height=8)#Endo
FeaturePlot(combined_Niche_S1, features = c('CCL19'))

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combined_Niche_S1, reduction = 'umap',label = T, group.by = 'Phase')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combined_Niche_S1, reduction = 'umap',label = T, group.by = 'orig_idents')+
scale_color_manual(values = color_used)





combined_Niche_S1 <- FindNeighbors(combined_Niche_S1, reduction = 'harmony', dims = 1:10, verbose = F)

combined_Niche_S1 <- FindClusters(combined_Niche_S1, resolution = seq(0.2, 1, 0.05), verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('symbol_snn_res.', i)
    p <- DimPlot(combined_Niche_S1, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)+
    scale_color_manual(values = color_used)

    print(p)
    print(res)
}

Idents(combined_Niche_S1) <- combined_Niche_S1$`symbol_snn_res.0.2`
options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 10)
combined_Niche_S1_markers <- FindAllMarkers(combined_Niche_S1, pseudocount.use = 0, logfc.threshold = log(1.5))
plan('sequential')
saveRDS(combined_Niche_S1_markers, file = 'combined_Niche_S1_markers.rds')
openxlsx::write.xlsx(combined_Niche_S1_markers, file = 'combined_Niche_S1_markers.xlsx')

options(repr.plot.width=16, repr.plot.height=16)
for(i in unique(combined_Niche_S1_markers$cluster)){
    genes <- combined_Niche_S1_markers[combined_Niche_S1_markers$cluster == i & 
                                   combined_Niche_S1_markers$avg_logFC > log(1.5) &
                                       combined_Niche_S1_markers$pct.1 > 0.2 & 
                          combined_Niche_S1_markers$p_val_adj < 0.05, 'gene'][1:16]
    genes <- genes[!is.na(genes)]
    if(length(genes) > 0){
        p <- FeaturePlot(combined_Niche_S1, reduction = 'umap', cols = c('lightgrey', 'red'),#raster = F,
                features = genes)
        print(p)
        print(paste0('cluster: ', i)) 
    }
}   


combined_Niche_S1$Identity <- as.character(combined_Niche_S1$`symbol_snn_res.0.7`)
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('0','10','17',
                                                            '3','5','6','15','16','18','19','20',
                                                            '12')]<- 'doublets'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('2' ,'7')]<- 'MKI67+cTEC'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('11')]<- 'MKI67+mTEC'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('4','13','14')]<- 'MKI67+Endo'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('1')]<- 'MKI67+Fb'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('8')]<- 'MKI67+VSMCs'
combined_Niche_S1$Identity[combined_Niche_S1$Identity %in% c('9')]<- 'Myelin'
Idents(combined_Niche_S1) <- combined_Niche_S1$Identity
##combined_Niche_S1 <- Idents_proj(combined_Niche_S1_sub, combined_Niche_S1)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combined_Niche_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 7)))









combined_Niche$Age_group <- ''
combined_Niche$Age_group[grep('w', combined_Niche$Age)] <- 'Fetal'
combined_Niche$Age_group[grep('d|m', combined_Niche$Age)] <- 'Prepuberal'
combined_Niche$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche$Age )) <= 12), 
                                  grep('y', combined_Niche$Age))] <- 'Prepuberal'
combined_Niche$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche$Age )) <= 39 &
                                   as.numeric(gsub('[a-z]', '',combined_Niche$Age )) >= 13), 
                                  grep('y', combined_Niche$Age))] <- 'Adult'
combined_Niche$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche$Age )) >= 40), 
                                  grep('y', combined_Niche$Age))] <- 'Aged'

combined_TEC$Age_group <- ''
combined_TEC$Age_group[grep('w', combined_TEC$Age)] <- 'Fetal'
combined_TEC$Age_group[grep('d|m', combined_TEC$Age)] <- 'Prepuberal'
combined_TEC$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_TEC$Age )) <= 12), 
                                  grep('y', combined_TEC$Age))] <- 'Prepuberal'
combined_TEC$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_TEC$Age )) <= 39 &
                                   as.numeric(gsub('[a-z]', '',combined_TEC$Age )) >= 13), 
                                  grep('y', combined_TEC$Age))] <- 'Adult'
combined_TEC$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_TEC$Age )) >= 40), 
                                  grep('y', combined_TEC$Age))] <- 'Aged'

combined_Niche_S1$Age_group <- ''
combined_Niche_S1$Age_group[grep('w', combined_Niche_S1$Age)] <- 'Fetal'
combined_Niche_S1$Age_group[grep('d|m', combined_Niche_S1$Age)] <- 'Prepuberal'
combined_Niche_S1$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche_S1$Age )) <= 12), 
                                  grep('y', combined_Niche_S1$Age))] <- 'Prepuberal'
combined_Niche_S1$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche_S1$Age )) <= 39 &
                                   as.numeric(gsub('[a-z]', '',combined_Niche_S1$Age )) >= 13), 
                                  grep('y', combined_Niche_S1$Age))] <- 'Adult'
combined_Niche_S1$Age_group[intersect(which(as.numeric(gsub('[a-z]', '',combined_Niche_S1$Age )) >= 40), 
                                  grep('y', combined_Niche_S1$Age))] <- 'Aged'









combined_TEC_2 <- merge(combined_TEC[, !Idents(combined_TEC) %in% c('doublets', 'Mesothelial')], 
                       combined_Niche_S1[, grepl('TEC', Idents(combined_Niche_S1))])

combined_TEC_2 <- FindVariableFeatures(combined_TEC_2, selection.method = "vst", nfeatures = 2000)
combined_TEC_2 <- ScaleData(combined_TEC_2)#, vars.to.regress = c('S.Score', 'G2M.Score'))
combined_TEC_2 <- RunPCA(combined_TEC_2, npcs = 50)
combined_TEC_2 <- RunHarmony(combined_TEC_2, 'batch', assay.use = 'symbol')

combined_TEC_2 <- RunUMAP(combined_TEC_2, reduction = "harmony", dims = 1:5, verbose = F)
# combined_TEC_2 <- FindNeighbors(combined_TEC_2, reduction = "harmony", dims = 1:5)

# combined_TEC_2 <- RunTSNE(combined_TEC_2, reduction = "harmony", dims = 1:5)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combined_TEC_2, reduction = 'umap',label = T)

options(repr.plot.width=8, repr.plot.height=13.5)
p <- DotPlot(combined_TEC_2, features = c('SIX1','KRT5','PAX9','KRT8',
                     'CCL25','PSMB11','HLA-DQB1','PRSS16',#'PAX9',
                       'FOXN1',
                                          'CXCL12',#'ZBED2',
                       #'PAX1',
                                          'DLL4',
                                          'MKI67','PCNA','CDK1',
'HMGB1',#'IGFBP6','NNMT','IGFBP5','MAOA','DPYS','FKBP5','GLUL',
            'EPCAM','ASCL1','KRT14','CCL19','KRT15','LYPD1','CCL21',
'AIRE',#'CD52',
                     'FEZF2','CLDN4','SPIB','FXYD3','IVL','KRT1','LYPD2','POU2F3','GNAT3', 'PLCB2',
                     'FOXI1', 'CFTR', 'ASCL3',
                    # 'MPZ','TGFBR2','SOX10',
                'MYOG','DES','MYLPF','BEX1','NEUROD1','CHGA',
'ATOH1','GFI1',
'LHX3','FOXJ1'),#stacked=T,pt.size=0, 
       cols = c('lightgrey', 'red'))

new_p <- ggplot(p$data, aes(id, features.plot, fill = avg.exp.scaled, size = pct.exp))+
geom_point(shape = 21)+
scale_size_continuous(range = c(2,9))+
scale_fill_gradientn(colours = c('white', 'red')) + 
theme(axis.text.x = element_text(size = 17.5, family = 'sans', angle =45, hjust = 1), 
      panel.background = element_rect(color = NA, fill ='white'),
      panel.border = element_rect(color = 'black', fill =NA),
      legend.box = 'vertical',
      legend.position = 'right',
      axis.ticks = element_blank(), axis.title = element_blank(), 
     axis.text.y = element_text(size = 17.5, family = 'sans', face = 'italic'),
      legend.text = element_text(size = 17.5, family = 'sans'),
      legend.title = element_text(size = 17.5, family = 'sans'),
      #strip.background.y = element_blank(), 
     #strip.text.y = element_text(size = 17.5, hjust = 0, family = 'sans', face = 'italic'),

      axis.line =element_blank(),
     plot.margin  = unit(c(0,0,0,2), 'lines'))+
guides(
       fill = guide_colorbar(title = "Expression", order = 1),
size = guide_legend(title = "% Cells", override.aes = list(fill = 'black'), order = 2))+
scale_x_discrete(limits = c('cTEC_hi', 'cTEC_lo', 'MKI67+cTEC',
                            'Immature_TEC', 'mTEC_lo', 'mTEC_hi',
       'post_AIRE_mTEC', 'MKI67+mTEC', 'Tuft', 'Ionocyte',  'TEC(myo)', 
                            'TEC(neuron)', 'Ciliated'))

new_p







combined_Niche$Identity <- as.character(combined_Niche$`symbol_snn_res.0.5`)
combined_Niche$Identity[combined_Niche$Identity %in% c('12', '13')]<- 'doublets'
combined_Niche$Identity[combined_Niche$Identity %in% c('15')]<- 'MKI67+Fb'
combined_Niche$Identity[combined_Niche$Identity %in% c('1' ,'4')]<- 'Endo'
combined_Niche$Identity[combined_Niche$Identity %in% c('5', '11')]<- 'VSMCs'
combined_Niche$Identity[combined_Niche$Identity %in% c('0')]<- 'Fb_1'
combined_Niche$Identity[combined_Niche$Identity %in% c('2')]<- 'Fb_2'

Idents(combined_Niche) <- combined_Niche$Identity
combined_Niche <- Idents_proj(combined_TEC, combined_Niche)
combined_Niche <- Idents_proj(combined_Niche_S1, combined_Niche)
options(repr.plot.width=12, repr.plot.height=10)
DimPlot(combined_Niche, reduction = 'umap',label = T, label.size = 15/.pt)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 15))+
guides(color = guide_legend(override.aes = list(size = 7)))



filtered_combined_Niche_2 <- combined_Niche[, !Idents(combined_Niche) %in% c('doublets') & 
                                           !combined_Niche$Age_group %in% c('Fetal')]

filtered_combined_Niche_2 <- RunUMAP(filtered_combined_Niche_2, reduction = 'harmony', dims = 1:10, verbose = F)

options(repr.plot.width=12, repr.plot.height=10)
DimPlot(filtered_combined_Niche_2, reduction = 'umap',label = T, label.size = 15/.pt, repel = T)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 15))+
guides(color = guide_legend(override.aes = list(size = 7)))

options(repr.plot.width=12, repr.plot.height=10)
DimPlot(filtered_combined_Niche_2, reduction = 'umap',label = T, label.size = 15/.pt, repel = T)+
scale_color_manual(values = color_used)+
theme(legend.text = element_text(size = 15))+
guides(color = guide_legend(override.aes = list(size = 7)))+
scale_color_manual(values = cluster_colors)











age_colors <- c('Fetal'='#cfbdd8', 'Prepuberal'="#b7d9ee", 'Adult'="#3ca88e", 'Aged'="#ec8f46")

filtered_TRA_TEC <- combined_TEC_2[, !combined_TEC_2$Age_group %in% c('Fetal') & 
                           !Idents(combined_TEC_2) %in% 
                           c('cTEC_hi','cTEC_lo','Immature_TEC','MKI67+cTEC')]

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(filtered_TRA_TEC, ndims = 50, reduction = 'harmony')

filtered_TRA_TEC <- RunUMAP(filtered_TRA_TEC, reduction = "harmony", dims = 1:5, verbose = F)


options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(filtered_TRA_TEC, reduction = 'umap',label = T)+
scale_color_manual(values =color_used)

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(filtered_TRA_TEC, reduction = 'umap',label = T, split.by = 'Age_group')+
scale_color_manual(values =color_used)



filtered_combined_TEC_2 <- combined_TEC_2[, !combined_TEC_2$Age_group %in% c('Fetal')]

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(filtered_combined_TEC_2, ndims = 50, reduction = 'harmony')

filtered_combined_TEC_2 <- RunUMAP(filtered_combined_TEC_2, reduction = "harmony", dims = 1:5, verbose = F)

filtered_combined_TEC_2 <- RunTSNE(filtered_combined_TEC_2, reduction = "harmony", dims = 1:5, verbose = F)

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(filtered_combined_TEC_2, reduction = 'tsne',label = T)+
scale_color_manual(values =color_used)

options(repr.plot.width=7.5, repr.plot.height=6)
DimPlot(filtered_combined_TEC_2, reduction = 'umap',label = T, label.size = 15/.pt)+
scale_color_manual(values =color_used)+
scale_color_manual(values =cluster_colors)+d_theme_w(size = 15)











saveRDS(combined_Niche, file = 'combined_Niche.rds')

saveRDS(combined_TEC, file = 'combined_TEC.rds')

save(combined_TEC_S1, combined_TEC_S1_S1, file = 'combined_TEC_subsets.RData')

saveRDS(combined_Niche_S1, file = 'combined_Niche_S1.rds')

saveRDS(combined_TEC_2, file = 'combined_TEC_2.rds')

saveRDS(filtered_combined_Niche_2, file = 'filtered_combined_Niche_2.rds')

saveRDS(filtered_TRA_TEC, file = 'filtered_TRA_TEC.rds')

saveRDS(filtered_combined_TEC_2, file = 'filtered_combined_TEC_2.rds')


