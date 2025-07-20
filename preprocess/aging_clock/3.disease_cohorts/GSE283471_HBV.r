library(SCopeLoomR)
library(Seurat)
library(ggplot2)
library(Matrix)
library(harmony)
library(plyr)
library(pheatmap)
library(future)
library(patchwork)


source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


library(glmnet)

library(data.table)

library(DoubletFinder)

library(clusterProfiler)

library(openxlsx)

library(future)

library(future.apply)



TN <- readRDS('TN_GSE283471.rds')

features_df_V3 <- 
  read.table('/data1/01.project/Thymus/thymus_first_data/RNA/result/T201910161F-10X5/filtered_feature_bc_matrix/features.tsv.gz')

features_df_V3 <- features_df_V3[, c('V1', 'V2')]
colnames(features_df_V3) <- c('ensembl', 'symbol')
features_df_V3$symbol <- make.unique(features_df_V3$symbol)

features_df_V6 <- read.table('/data2/thymus_PB_plus/ABFC20211317-040511-SCS-result/2.Cellranger/P2021113001M-10X5/filtered_feature_bc_matrix/features.tsv.gz')
features_df_V6 <- features_df_V6[, c('V1', 'V2')]
colnames(features_df_V6) <- c('ensembl', 'symbol')
features_df_V6$symbol <- make.unique(features_df_V6$symbol)

features_df_2022_NA <- read.table('/data1/02.private/dengyj/analysis/thymus/Integrated/2022_NA_aging_PBMC/data/GSM5684308_OH17_features.tsv.gz')

features_df_2022_NA <- features_df_2022_NA[, c(1,2)]
colnames(features_df_2022_NA) <- c('ensembl', 'symbol')

features_df_2022_NA$symbol <- make.unique(features_df_2022_NA$symbol)

combined_features_df <- read.csv('/data1/02.private/dengyj/analysis/thymus/doublets/combined_features_df.csv', 
                                 row.names = 1)

features_df_GSE161918 <- read.table('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE161918/features_df.tsv')

features_df_GSE161918 <- features_df_GSE161918[, c(1,3)]
colnames(features_df_GSE161918) <- c('ensembl', 'symbol')

features_df_GSE161918$symbol <- make.unique(features_df_GSE161918$symbol)

all(rownames(TN) %in% features_df_V6$symbol)

exp <- TN$RNA@counts

if(all(rownames(exp) %in% features_df_V6$symbol)){
   rownames(exp) <- features_df_V6$ensembl[match(rownames(exp), features_df_V6$symbol)] 
}

all(rownames(exp) %in% features_df_V6$ensembl)

TN <- CreateSeuratObject(counts = exp, meta.data = TN@meta.data)

TN <- NormalizeData(TN, verbose = F)

all(table(TN$batch) > 50)

load('/data1/02.private/dengyj/analysis/thymus/Integrated/aging_clock/250404_final_model.RData')

exp_mean_orig <- my_group_mean(expm1(TN$RNA@data), factor(TN$batch))

exp_mean <- log1p(t(exp_mean_orig))

metadata <- TN[[]]
metadata <- metadata[!duplicated(metadata$batch),]
rownames(metadata) <- metadata$batch
metadata <- metadata[rownames(exp_mean),]

coef_mat <- coef(final_model)[,1]
coef_val <- names(coef_mat)#[coef_mat > 0]
coef_val <- setdiff(coef_val, '(Intercept)')
setdiff_val <- setdiff(coef_val,colnames(exp_mean))

exp_mean <- cbind(exp_mean, 
                  matrix(0, dimnames = list(rownames(exp_mean), setdiff_val),
                         ncol = length(setdiff_val), 
                                   nrow = nrow(exp_mean)))

ind_plot_df <- data.frame(predict(final_model, newx = exp_mean[,coef_val]), metadata$Age)
colnames(ind_plot_df) <- c('pred', 'true')  
ind_plot_df$pred_adj <- ind_plot_df$pred - predict(loess_m, ind_plot_df$true)

options(repr.plot.width=5, repr.plot.height=5)
ind_plot_df <- ind_plot_df
ind_plot_df$batch <- rownames(ind_plot_df)
ind_plot_df$source <- NULL


ind_plot_df <- ind_plot_df
ind_plot_df <- melt(ind_plot_df, id.vars = 'batch', variable.name = 'class' , value.name = 'Age')
ind_plot_df$class <- factor(ind_plot_df$class, levels=c('true','pred', 'pred_adj'))
ind_plot_df$group <- TN$Sample_title[match(ind_plot_df$batch, TN$batch)]
#ind_plot_df <- ind_plot_df[!grepl('AHB', ind_plot_df$group),]

options(repr.plot.width=2.5, repr.plot.height=4.25)
ggplot(ind_plot_df[!ind_plot_df$class %in%c('pred') ,], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0.1, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Real', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'HBV cohort')+
theme(plot.title = element_text(hjust= 0.5))



saveRDS(ind_plot_df, file = 'HBV_ind_plot_df.rds')
