library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(data.table)
library(openxlsx)
library(ggpubr)
library(harmony)

library(Rmisc)
library(patchwork)

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle2.R')


library(AUCell)
library(SCopeLoomR)
library(SCENIC)

source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

TRAs_symbol=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_symbol.RData')

TRAs_df=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_df.RData')

TRAs_tau=load('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/TRAs_tau.RData')

TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')



filtered_TRA_TEC <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_TRA_TEC.rds')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(filtered_TRA_TEC, label=T, reduction = 'umap')



####using gene ensembl id
DefaultAssay(filtered_TRA_TEC) <- 'RNA'
filtered_TRA_TEC$TRAs_consensus <- 
PercentageFeatureSet(filtered_TRA_TEC, features = 
                     intersect(unique(TRAs_consensus_df[TRAs_consensus_df$tau >=0.8, 'ensembl']), 
                               rownames(filtered_TRA_TEC)))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(filtered_TRA_TEC, reduction = 'umap',min.cutoff = 'q5',max.cutoff = 'q95',
            features = c('TRAs_consensus'), cols = c('lightgrey', 'red'))



TRAs_consensus_df$tissue <- gsub(' ', '_', TRAs_consensus_df$tissue)
TRAs_consensus_df$tissue <- gsub('\\-', '_', TRAs_consensus_df$tissue)
TRAs_consensus_df$tissue <- gsub(',', '', TRAs_consensus_df$tissue)

DefaultAssay(filtered_TRA_TEC) <- 'RNA'
filtered_TRA_TEC$TRAs <- 
PercentageFeatureSet(filtered_TRA_TEC, features = 
                     intersect(unique(TRAs_consensus_df[TRAs_consensus_df$tau >=0.8, 'ensembl']), 
                               rownames(filtered_TRA_TEC)))

tissue <- unique(TRAs_consensus_df$tissue)

for(i in tissue){
    assign('features_iter', intersect(unique(TRAs_consensus_df[TRAs_consensus_df$tissue == i & 
                                                         TRAs_consensus_df$tau >= 0.8, 'ensembl']), 
                                      rownames(filtered_TRA_TEC)))
    filtered_TRA_TEC[[i]] <- PercentageFeatureSet(filtered_TRA_TEC, features = features_iter)
}




markers_TEC_myo <- FindMarkers(filtered_TRA_TEC, `ident.1` = c('TEC(myo)'), pseudocount.use = 0, logfc.threshold = log(1.5))

saveRDS(markers_TEC_myo, file = 'markers_TEC_myo.rds')

saveRDS(filtered_TRA_TEC, file = 'filtered_TRA_TEC_process.rds')
