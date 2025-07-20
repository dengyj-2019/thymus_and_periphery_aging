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

library(openxlsx)

library(data.table)

library(DoubletFinder)

library(fgsea)

color_used <- c('#E5D2DD', '#53A85F', '#F1BB72',  '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#F3B1A0',
'#968175',"#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4"
)


library(future.apply)

options(future.globals.maxSize = 10 * 1024^3)
plan('multicore', workers = 4)

# library(clusterProfiler)

library(openxlsx)


TN_GSE135779 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE135779/TN_GSE135779.rds')

TN_GSE283471 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE283471/TN_GSE283471.rds')

TN_GSE227835 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE227835/TN_GSE227835.rds')

TN_GSE243905 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE243905/TN_GSE243905.rds')

head(rownames(TN_GSE135779))

head(rownames(TN_GSE283471))

head(rownames(TN_GSE243905))







HIV_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE243905/HIV_ind_plot_df.rds')
HBV_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE283471/HBV_ind_plot_df.rds')

HIV_ind_plot_df$disease <- 'HIV'

HBV_ind_plot_df$disease <- 'HBV'

MG_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE227835/MG_ind_plot_df.rds')
SLE_ind_plot_df <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE135779/SLE_ind_plot_df.rds')


MG_ind_plot_df$disease <- 'MG'

SLE_ind_plot_df$disease <- 'SLE'

chronic_infection_df <- rbind(HIV_ind_plot_df, HBV_ind_plot_df)

autoimmune_disease_df <- rbind(MG_ind_plot_df, SLE_ind_plot_df)

table(autoimmune_disease_df$disease)/3



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





obj_list <- c('TN_GSE135779','TN_GSE227835','TN_GSE243905',
              'TN_GSE283471')

if(exists(x = 'exp_mean_list')){
    exp_mean_list <- exp_mean_list
}else{
    exp_mean_list <- list()
}
for(obj_name in setdiff(obj_list, names(exp_mean_list))){
    obj <- get(obj_name)
    exp <- expm1(obj$RNA@data)
    exp_mean <- future_lapply(unique(obj$batch), function(sample){
        exp_iter <- exp[, obj$batch %in% sample, drop = F]
        val <- rowMeans(exp_iter)
    })
    exp_mean <- do.call(cbind, exp_mean)
    colnames(exp_mean) <- unique(obj$batch)
    exp_mean_list[[obj_name]] <- exp_mean  
}



ensembl_used <- c()
for(i in c('TN_GSE135779','TN_GSE227835','TN_GSE283471')){
    ensembl_used <- c(ensembl_used, rownames(exp_mean_list[[i]]))
}
ensembl_used <- unique(ensembl_used)

map_df <- data.frame(ensembl =ensembl_used)
map_df$symbol <- combined_features_df$symbol[match(map_df$ensembl, combined_features_df$ensembl)]
# logi <-is.na(map_df$symbol)
# map_df$symbol[logi] <- features_df_V3$symbol[match(map_df$ensembl[logi], features_df_V3$ensembl)]

logi <-is.na(map_df$symbol)
map_df$symbol[logi] <- features_df_2022_NA$symbol[match(map_df$ensembl[logi], features_df_2022_NA$ensembl)]


logi <-is.na(map_df$symbol)
map_df$symbol[logi] <- features_df_GSE161918$symbol[match(map_df$ensembl[logi], 
                                                          features_df_GSE161918$ensembl)]

logi <-is.na(map_df$symbol)
map_df$symbol[logi] <- map_df$ensembl[logi]  



map_df$symbol[grepl('ENSG', map_df$symbol)]

for(i in c('TN_GSE135779','TN_GSE227835',#'TN_GSE243905',
              'TN_GSE283471')){
    rownames(exp_mean_list[[i]]) <- map_df$symbol[match(rownames(exp_mean_list[[i]]), 
                                                        map_df$ensembl)]
}

if(exists(x = 'TN_exp_mean_orig')){
    TN_exp_mean_orig <- TN_exp_mean_orig
    logi <- sapply(exp_mean_list, function(x){
        all(colnames(x) %in% colnames(TN_exp_mean_orig))
    })
    obj_names <- names(exp_mean_list)[!logi]
    if(length(obj_names) > 0){
        for(obj_name in obj_names){
            df_iter <- data.frame(exp_mean_list[[obj_name]], check.rows = F, check.names = F)
            TN_exp_mean_orig <- merge(TN_exp_mean_orig, df_iter, by = 'row.names', all = T)
            rownames(TN_exp_mean_orig) <- TN_exp_mean_orig$Row.names
            TN_exp_mean_orig$Row.names <- NULL
            TN_exp_mean_orig[is.na(TN_exp_mean_orig)] <- 0
        }           
    }
 
}else{
    for(i in 1:length(exp_mean_list)){
        if(i==1){
            TN_exp_mean_orig <- data.frame(exp_mean_list[[i]], check.rows = F, check.names = F)
        }else{
            df_iter <- data.frame(exp_mean_list[[i]], check.rows = F, check.names = F)
            TN_exp_mean_orig <- merge(TN_exp_mean_orig, df_iter, by = 'row.names', all = T)
            rownames(TN_exp_mean_orig) <- TN_exp_mean_orig$Row.names
            TN_exp_mean_orig$Row.names <- NULL
            TN_exp_mean_orig[is.na(TN_exp_mean_orig)] <- 0
        }
    }
}


TN_exp_mean <- log1p(t(TN_exp_mean_orig))

chronic_infection_df_d <- dcast(chronic_infection_df, batch~class, value.var = 'Age')
chronic_infection_residuals <- chronic_infection_df_d$pred_adj-chronic_infection_df_d$true

chronic_infection_cor_residuals <- psych::corr.test(chronic_infection_residuals,
                                  TN_exp_mean[chronic_infection_df_d$batch,])

r <- t(chronic_infection_cor_residuals$r);r[is.na(r)] <- 0
p.adj <- t(chronic_infection_cor_residuals$p.adj);p.adj[is.na(p.adj)] <- 1
p <- t(chronic_infection_cor_residuals$p);p[is.na(p)] <- 1

chronic_infection_cor_residuals_df <- data.frame(r , p,  p.adj)
colnames(chronic_infection_cor_residuals_df) <- c('r', 'p', 'p.adj')



chronic_infection_cor_residuals_df <- chronic_infection_cor_residuals_df[order(chronic_infection_cor_residuals_df$r, decreasing = T), ]
chronic_infection_cor_residuals_df$genes <- factor(rownames(chronic_infection_cor_residuals_df), 
                                                   levels = rev(unique(rownames(chronic_infection_cor_residuals_df))))
chronic_infection_cor_residuals_df$class <- ifelse(chronic_infection_cor_residuals_df$r > 0, 'yes', 'no')
chronic_infection_cor_residuals_df$class <- factor(chronic_infection_cor_residuals_df$class, levels = c('yes','no'))

openxlsx::write.xlsx(chronic_infection_cor_residuals_df,
                     file ='chronic_infection_cor_residuals_df.xlsx', row.names=T)

saveRDS(chronic_infection_cor_residuals_df,file = 'chronic_infection_cor_residuals_df.rds')

autoimmune_disease_df_d <- dcast(autoimmune_disease_df, batch~class, value.var = 'Age')
autoimmune_disease_residuals <- autoimmune_disease_df_d$pred_adj-autoimmune_disease_df_d$true

autoimmune_disease_cor_residuals <- psych::corr.test(autoimmune_disease_residuals,
                                  TN_exp_mean[autoimmune_disease_df_d$batch,])

r <- t(autoimmune_disease_cor_residuals$r);r[is.na(r)] <- 0
p.adj <- t(autoimmune_disease_cor_residuals$p.adj);p.adj[is.na(p.adj)] <- 1
p <- t(autoimmune_disease_cor_residuals$p);p[is.na(p)] <- 1

autoimmune_disease_cor_residuals_df <- data.frame(r , p,  p.adj)
colnames(autoimmune_disease_cor_residuals_df) <- c('r', 'p', 'p.adj')




autoimmune_disease_cor_residuals_df <- autoimmune_disease_cor_residuals_df[order(autoimmune_disease_cor_residuals_df$r, decreasing = T), ]
autoimmune_disease_cor_residuals_df$genes <- factor(rownames(autoimmune_disease_cor_residuals_df), 
                                                   levels = rev(unique(rownames(autoimmune_disease_cor_residuals_df))))
autoimmune_disease_cor_residuals_df$class <- ifelse(autoimmune_disease_cor_residuals_df$r > 0, 'yes', 'no')
autoimmune_disease_cor_residuals_df$class <- factor(autoimmune_disease_cor_residuals_df$class, levels = c('yes','no'))

openxlsx::write.xlsx(autoimmune_disease_cor_residuals_df,
                     file ='autoimmune_disease_cor_residuals_df.xlsx', row.names=T)

saveRDS(autoimmune_disease_cor_residuals_df,file = 'autoimmune_disease_cor_residuals_df.rds')

save(exp_mean_list, map_df,TN_exp_mean_orig, 
     TN_exp_mean, 
     chronic_infection_df,
     autoimmune_disease_df, file = 'chronic_disease.RData')

####GSEA
chronic_infection_pcc <- chronic_infection_cor_residuals_df$r
names(chronic_infection_pcc) <- chronic_infection_cor_residuals_df$genes
chronic_infection_pcc <- sort(chronic_infection_pcc, decreasing = T)

write.table(data.frame(names(chronic_infection_pcc), chronic_infection_pcc), 
            file = 'chronic_infection_pcc.rnk', quote = F,sep = '\t',
            col.names = F, row.names = F)

#Linux environment
gsea-cli.sh GSEAPreranked -rnk chronic_infection_pcc.rnk \
-gmx /data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/GO_KEGG_REACTOME.gmt \
-collapse No_Collapse \
-set_max 10000 \
-rpt_label chronic_infection_pcc

####GSEA
autoimmune_disease_pcc <- autoimmune_disease_cor_residuals_df$r
names(autoimmune_disease_pcc) <- autoimmune_disease_cor_residuals_df$genes
autoimmune_disease_pcc <- sort(autoimmune_disease_pcc, decreasing = T)

write.table(data.frame(names(autoimmune_disease_pcc), autoimmune_disease_pcc), 
            file = 'autoimmune_disease_pcc.rnk', quote = F,sep = '\t',
            col.names = F, row.names = F)

#Linux environment
gsea-cli.sh GSEAPreranked -rnk autoimmune_disease_pcc.rnk \
-gmx /data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/GO_KEGG_REACTOME.gmt \
-collapse No_Collapse \
-set_max 10000 \
-rpt_label autoimmune_disease_pcc






