library(SCopeLoomR)
library(Seurat)
library(ggplot2)
library(Matrix)
library(harmony)
library(plyr)
library(future)
library(pheatmap)
library(patchwork)



source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')



library(circlize)

library(ComplexHeatmap)


source('/data1/02.private/dengyj/analysis/thymus/plot/support//plot_helper.r')


age_colors <- c('<=12'="#b7d9ee",'Prepuberal' = '#b7d9ee','Prepubertal' = '#b7d9ee',
                '13-39'="#3ca88e",'Adult' = '#3ca88e', 
                'Group1' = '#b7d9ee', 'Group2' = '#3ca88e', 'Group3' = '#ec8f46',
                '>=40'="#ec8f46", '40-99'="#ec8f46",'Aged' = '#ec8f46', '>=100' = '#c04e5c', 
                'Centenarian' = '#c04e5c')


library(openxlsx)

library(clusterProfiler)

library(reshape2)

library(ggpubr)

library(future.apply)

library(glmnet)


color_used <- c('#E5D2DD', '#53A85F', '#F1BB72',  '#57C3F3', '#476D87',
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

options(future.globals.maxSize = 100 * 1024^3)
plan('multicore', workers = 16)

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



CD8TN_2023_Imm <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_Immunity_aging_PBMC/CD8TN.rds')

CD8TN_2023_Imm$identity <- as.character(CD8TN_2023_Imm$Cluster_names)
CD8TN_2023_Imm$source <- '2023_Imm'
CD8TN_2023_Imm$batch <- paste0(CD8TN_2023_Imm$Donor_id, '_', CD8TN_2023_Imm$Age)
# CD8TN_2023_Imm$batch2 <- CD8TN_2023_Imm$Batch

CD4TN_2023_Imm <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_Immunity_aging_PBMC/CD4TN.rds')

CD4TN_2023_Imm$identity <- as.character(CD4TN_2023_Imm$Cluster_names)
CD4TN_2023_Imm$source <- '2023_Imm'
CD4TN_2023_Imm$batch <- paste0(CD4TN_2023_Imm$Donor_id, '_', CD4TN_2023_Imm$Age)
# CD4TN_2023_Imm$batch2 <- CD4TN_2023_Imm$Batch

TN_2023_Imm <- merge(CD4TN_2023_Imm, CD8TN_2023_Imm)

all(rownames(TN_2023_Imm) %in% features_df_V6$symbol)

exp <- TN_2023_Imm$RNA@counts

if(all(rownames(exp) %in% features_df_V6$symbol)){
   rownames(exp) <- features_df_V6$ensembl[match(rownames(exp), features_df_V6$symbol)] 
}

all(rownames(exp) %in% features_df_V6$ensembl)

TN_2023_Imm <- CreateSeuratObject(counts = exp, meta.data = TN_2023_Imm@meta.data)

TN_2023_Imm <- NormalizeData(TN_2023_Imm, verbose = F)



Tcell <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell.rds')

Tcell$Age <- Tcell$age


Tcell$source <- 'thymus_aging'

TN_this_study <- Tcell[, grepl('RTE|TN',Idents(Tcell))]

all(rownames(TN_this_study$RNA@counts) %in% combined_features_df$symbol)

exp <- TN_this_study$RNA@counts

if(all(rownames(exp) %in% combined_features_df$symbol)){
   rownames(exp) <- combined_features_df$ensembl[match(rownames(exp), combined_features_df$symbol)] 
}

all(rownames(exp) %in% combined_features_df$ensembl)

TN_this_study <- CreateSeuratObject(counts = exp, meta.data = TN_this_study@meta.data)

TN_this_study <- NormalizeData(TN_this_study, verbose = F)



TN_2019_PNAs <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2019_PNAs_centenarian//centenarian_TN.rds')

TN_2019_PNAs$Age <- as.numeric(gsub('s', '', TN_2019_PNAs$Age))
TN_2019_PNAs$age <- TN_2019_PNAs$Age
TN_2019_PNAs$identity <- 'Naive'
TN_2019_PNAs$source <- '2019_PNAs'
TN_2019_PNAs$batch <- TN_2019_PNAs$ID

# TN_2019_PNAs$batch2 <- TN_2019_PNAs$source

TN_2022_NA <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2022_NA_aging_PBMC/TN.rds')

TN_2022_NA$Age <- TN_2022_NA$age
TN_2022_NA$identity <- as.character(Idents(TN_2022_NA))

TN_2022_NA$source <- '2022_NA'

# TN_2022_NA$batch2 <- TN_2022_NA$source





exp <- TN_2019_PNAs$RNA@counts

logi <- grepl('--', rownames(exp))
rownames(exp)[logi] <- gsub('--', '__', rownames(exp)[logi])
logi <- grepl('-A\\.2', rownames(exp))
rownames(exp)[logi] <- gsub('-A\\.2', '_A.2', rownames(exp)[logi])

if(all(rownames(exp) %in% features_df_2022_NA$symbol)){
   rownames(exp) <- features_df_2022_NA$ensembl[
       match(rownames(exp), features_df_2022_NA$symbol)] 
}else{
    print('stop')
}

all(rownames(exp) %in% features_df_2022_NA$ensembl)

TN_2019_PNAs <- CreateSeuratObject(counts = exp, meta.data = TN_2019_PNAs@meta.data)

TN_2019_PNAs <- NormalizeData(TN_2019_PNAs, verbose = F)



exp <- TN_2022_NA$RNA@counts

logi <- grepl('--', rownames(exp))
rownames(exp)[logi] <- gsub('--', '__', rownames(exp)[logi])
logi <- grepl('-A\\.2', rownames(exp))
rownames(exp)[logi] <- gsub('-A\\.2', '_A.2', rownames(exp)[logi])

if(all(rownames(exp) %in% features_df_2022_NA$symbol)){
   rownames(exp) <- features_df_2022_NA$ensembl[
       match(rownames(exp), features_df_2022_NA$symbol)] 
}else{
    print('stop')
}

all(rownames(exp) %in% features_df_2022_NA$ensembl)

TN_2022_NA <- CreateSeuratObject(counts = exp, meta.data = TN_2022_NA@meta.data)

TN_2022_NA <- NormalizeData(TN_2022_NA, verbose = F)



CD4TN_2024_AIDA_V1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2024_AIDA_V1/CD4TN.rds')
CD4TN_2024_AIDA_V1$batch <- CD4TN_2024_AIDA_V1$donor_id
CD4TN_2024_AIDA_V1$identity <- 'CD4TN'
CD4TN_2024_AIDA_V1$source <- '2024_AIDA_V1'
CD4TN_2024_AIDA_V1$Age <- as.numeric(gsub('-year-old stage', '', CD4TN_2024_AIDA_V1$development_stage))
# CD4TN_2024_AIDA_V1$batch2 <- CD4TN_2024_AIDA_V1$source

CD8TN_2024_AIDA_V1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2024_AIDA_V1/CD8TN.rds')
CD8TN_2024_AIDA_V1$batch <- CD8TN_2024_AIDA_V1$donor_id
CD8TN_2024_AIDA_V1$identity <- 'CD8TN'
CD8TN_2024_AIDA_V1$source <- '2024_AIDA_V1'
CD8TN_2024_AIDA_V1$Age <- as.numeric(gsub('-year-old stage', '', CD8TN_2024_AIDA_V1$development_stage))
# CD8TN_2024_AIDA_V1$batch2 <- CD8TN_2024_AIDA_V1$source

TN_2024_AIDA_V1 <- merge(CD4TN_2024_AIDA_V1, CD8TN_2024_AIDA_V1)

TN_2024_AIDA_V1 <- NormalizeData(TN_2024_AIDA_V1, verbose = F)



TN_2025_NI <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2025_NI_Hs_PBMC_aging/TN.rds')

TN_2025_NI$identity <- as.character(TN_2025_NI$secondary_type)
TN_2025_NI$source <- '2025_NI'
TN_2025_NI$batch <- TN_2025_NI$sampleName
# TN_2025_NI$batch2 <- TN_2025_NI$source

all(rownames(TN_2025_NI) %in% features_df_V3$symbol)

exp <- TN_2025_NI$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

TN_2025_NI <- CreateSeuratObject(counts = exp, meta.data = TN_2025_NI@meta.data)

TN_2025_NI <- NormalizeData(TN_2025_NI, verbose = F)



TN_2023_SA <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_SA_aging_PBMC/2023_SA_aging_TN.rds')

TN_2023_SA$identity <- 'Naive'
TN_2023_SA$source <- '2023_SA'

TN_2023_SA$batch2 <- TN_2023_SA$source

all(rownames(TN_2023_SA) %in% features_df_V3$symbol)

exp <- TN_2023_SA$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

TN_2023_SA <- CreateSeuratObject(counts = exp, meta.data = TN_2023_SA@meta.data)

TN_2023_SA <- NormalizeData(TN_2023_SA, verbose = F)







TN_GSE158055 <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE158055/TN_GSE158055.rds')


all(rownames(TN_GSE158055) %in% features_df_V3$symbol)

exp <- TN_GSE158055$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

TN_GSE158055 <- CreateSeuratObject(counts = exp, 
                                               meta.data = TN_GSE158055@meta.data)

TN_GSE158055 <- NormalizeData(TN_GSE158055, verbose = F)



TN_GSE161918 <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE161918//GSE161918_TN.rds')


exp <- TN_GSE161918$RNA@counts

logi <- grepl('--', rownames(exp))
rownames(exp)[logi] <- gsub('--', '__', rownames(exp)[logi])
logi <- grepl('-A\\.2', rownames(exp))
rownames(exp)[logi] <- gsub('-A\\.2', '_A.2', rownames(exp)[logi])

logi <- grepl('Y-RNA', rownames(exp))
rownames(exp)[logi] <- gsub('Y-RNA', 'Y_RNA', rownames(exp)[logi])

if(all(rownames(exp) %in% features_df_GSE161918$symbol)){
   rownames(exp) <- features_df_GSE161918$ensembl[
       match(rownames(exp), features_df_GSE161918$symbol)] 
}else{
    print('stop')
}

all(rownames(exp) %in% features_df_GSE161918$ensembl)

TN_GSE161918 <- CreateSeuratObject(counts = exp, 
                                               meta.data = TN_GSE161918@meta.data)

TN_GSE161918 <- NormalizeData(TN_GSE161918, verbose = F)







obj_list <- c('TN_2023_Imm','TN_2019_PNAs', 'TN_2022_NA', 
'TN_2024_AIDA_V1', 'TN_2025_NI', 'TN_2023_SA',
             'TN_GSE158055','TN_this_study', 'TN_GSE161918')

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

saveRDS(exp_mean_list, file = 'exp_mean_list.rds')

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
saveRDS(TN_exp_mean_orig, file = 'TN_exp_mean_orig.rds')
TN_exp_mean <- t(TN_exp_mean_orig)

common_meta <- c('source', 'Age', 'batch')

meta_list <- lapply(obj_list, function(obj){
    obj <- get(obj)
    meta <- obj[[common_meta]]
})
meta_list <- do.call(rbind, meta_list)
saveRDS(meta_list, file = 'meta_list.rds')



TN_meta <- meta_list
TN_meta <- TN_meta[!duplicated(TN_meta$batch), c('source','Age','batch')]
rownames(TN_meta) <- TN_meta$batch

saveRDS(TN_meta, file = 'TN_meta.rds')

identical(rownames(TN_exp_mean),rownames(TN_meta))

samples_used <- unique(meta_list$batch)
samples_tbl <- table(meta_list$batch)
samples_used <- intersect(samples_used,
                          names(samples_tbl)[samples_tbl >= 50])

names(samples_tbl)[samples_tbl < 50]

modeling_samples <- samples_used


modeling_samples <- setdiff(modeling_samples, 
                            TN_meta$batch[TN_meta$source=='2019_PNAs' & TN_meta$Age < 100])



####log1p process
source_for_gene_detected <- c('2023_Imm', '2024_AIDA_V1', '2025_NI')
if(exists(x = 'lm_list')){
    lm_list <- lm_list
}else{
    lm_list <- list()
}
for(source in setdiff(source_for_gene_detected, names(lm_list))){
    logi <- TN_meta$source == source & TN_meta$batch %in% modeling_samples
    tmp_exp <- TN_exp_mean[logi, ]
    tmp_exp <- tmp_exp[, colSums(tmp_exp) > 0]
    data_used <- cbind(TN_meta[logi, ], tmp_exp)    
    lm_res <- future_lapply(colnames(tmp_exp), function(gene){
        formula <- as.formula(paste0('`', gene, '`', " ~ Age"))
        # 建立线性回归模型
        model <- lm(formula, data = data_used)
        pval <- summary(model)$coefficients[, 4][2]
        slope = coef(model)[2]
        #print(which(colnames(TN_exp_pseudobulk_orig)==gene))
        df <- data.frame(pval = pval, gene = gene, slope = slope)         
    })
    lm_res <- do.call(rbind, lm_res)
    lm_res$source <- source
    lm_res$padj <- p.adjust(lm_res$pval, method = 'fdr')
    lm_list[[source]] <-   lm_res  
}

saveRDS(lm_list, file = 'lm_list.rds')

source_for_gene_detected2 <- c('2024_AIDA_V1', '2023_Imm')

up_gene_list <- lapply(lm_list[source_for_gene_detected2], function(df){
    genes <- df$gene[df$padj < 0.01 & df$slope > 1e-3]
})
up_gene_list <- Reduce(intersect, up_gene_list)

down_gene_list <- lapply(lm_list[source_for_gene_detected2], function(df){
    genes <- df$gene[df$padj < 0.01 & df$slope < -1e-3]
})
down_gene_list <- Reduce(intersect, down_gene_list)

gene_used <- c(up_gene_list, down_gene_list)

length(gene_used)



source_for_modeling <- c('2024_AIDA_V1',  '2025_NI')

frail_cohort <- unique(TN_2022_NA$batch[TN_2022_NA$age_group=='frail'])
#covid19_cohort1 <- TN_meta$batch[TN_meta$source %in% c('2022_NC_COVID19_aging')]
centenarians_cohort <- TN_meta$batch[TN_meta$Age >=100]
covid19_cohort1 <- TN_GSE158055$batch[!TN_GSE158055$CoVID.19.severity %in% c('control')]
covid19_cohort2 <- TN_GSE161918$batch[!TN_GSE161918$Class %in% c('HC')]


ind_exclude_sample <- union(
    c(frail_cohort,covid19_cohort1,centenarians_cohort,covid19_cohort2), 
c('S-HC013','S-HC014',unique(TN_2023_Imm$batch)))

save(frail_cohort, centenarians_cohort, covid19_cohort1, covid19_cohort2,
     ind_exclude_sample, file = 'cohorts.RData')

#log1p process
logi <- TN_meta$source %in% source_for_modeling#
batch_for_modeling <- intersect(unique(TN_meta$batch[logi]), modeling_samples)
batch_for_modeling <- setdiff(batch_for_modeling, centenarians_cohort)
# batch_for_modeling <- setdiff(batch_for_modeling, covid19_cohort2)
set.seed(1234)
train_samples <- sample(batch_for_modeling, size = 0.7 * length(batch_for_modeling))
test_samples <- setdiff(batch_for_modeling, train_samples)

logi <- TN_meta$batch %in% train_samples
train_x <- log1p(TN_exp_mean[logi, gene_used])
train_y <- TN_meta$Age[logi]

logi <- TN_meta$batch %in% test_samples
test_x <- log1p(TN_exp_mean[logi, gene_used])
test_y <- TN_meta$Age[logi]



ind_samples <- setdiff(modeling_samples, batch_for_modeling)

logi <- TN_meta$batch %in% ind_samples
ind_x <- log1p(TN_exp_mean[logi, gene_used])
ind_y <- TN_meta$Age[logi]




alpha <- 0.5
seed <- 1
set.seed(seed)
cv_model <- cv.glmnet(x = train_x, y = train_y, alpha = alpha)
final_model <- glmnet(train_x, train_y, alpha = alpha, lambda = cv_model$lambda.min)

library(glmnet)


train_plot_df <- data.frame(predict(final_model, newx = train_x), train_y)
colnames(train_plot_df) <- c('pred', 'true')    
test_plot_df <- data.frame(predict(final_model, newx = test_x), test_y)
colnames(test_plot_df) <- c('pred', 'true')

 
ind_plot_df <- data.frame(predict(final_model, newx = ind_x), ind_y)
colnames(ind_plot_df) <- c('pred', 'true')  


train_val <- mean(abs(train_plot_df$pred-train_plot_df$true))
test_val <- mean(abs(test_plot_df$pred-test_plot_df$true))

ind_val <- mean(abs(ind_plot_df$pred-ind_plot_df$true))   

train_cor <- cor(train_plot_df$pred,train_plot_df$true)
test_cor <- cor(test_plot_df$pred,test_plot_df$true)

ind_cor <- cor(ind_plot_df$pred,ind_plot_df$true)
train_plot_df$source <- TN_meta$source[match(rownames(train_plot_df), TN_meta$batch)]
ind_plot_df$source <- TN_meta$source[match(rownames(ind_plot_df), TN_meta$batch)]


options(repr.plot.width=4, repr.plot.height=4)
ggplot(train_plot_df, 
       aes(true, pred))+
geom_point(mapping = aes( ))+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')+
d_theme_w(size =10)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(train_plot_df, 
       aes(true, pred-true))+
geom_point(mapping = aes( ))+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')+
d_theme_w(size =10)

###loess correction for over or under prediction
train_plot_df$age_diff <- train_plot_df$pred-train_plot_df$true
loess_m <- loess(age_diff~ true, data = train_plot_df, span = 0.75,
                control = loess.control(surface = "direct"))
train_plot_df$pred_adj <- train_plot_df$pred - predict(loess_m, train_plot_df$true)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(train_plot_df, 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5))+
labs(x = 'Actual age', y='Predicted age', title = 'Training cohort')
ggsave('TN_aging_clock_train_scatter.svg', width = 4, height = 4)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(train_plot_df, 
       aes(true, pred_adj-true))+
geom_point(mapping = aes())+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')+
d_theme_w(size =10)

healthy_ind_plot_df <- ind_plot_df[!rownames(ind_plot_df) %in% ind_exclude_sample,]

healthy_ind_plot_df$age_diff <- healthy_ind_plot_df$pred-healthy_ind_plot_df$true
healthy_ind_plot_df$pred_adj <- healthy_ind_plot_df$pred - predict(loess_m, healthy_ind_plot_df$true)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(healthy_ind_plot_df[ ], 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5),
     plot.margin = margin(r = 2.4, unit = 'mm'))+
labs(x = 'Actual age', y='Predicted age', title = 'Healthy cohort')
ggsave('TN_aging_clock_test_scatter.svg', width = 4, height = 4)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(healthy_ind_plot_df[ ], 
       aes(true, pred))+
geom_point(mapping = aes( ))+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')+
d_theme_w(size=10)

options(repr.plot.width=4, repr.plot.height=4)
ggplot(healthy_ind_plot_df, 
       aes(true, pred_adj-true))+
geom_point(mapping = aes())+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')+
d_theme_w(size=10)

options(repr.plot.width=2.5, repr.plot.height=4)
tmp_df <- healthy_ind_plot_df
tmp_df$batch <- rownames(tmp_df)
tmp_df$source <- NULL
tmp_df$age_diff <- NULL
tmp_df2 <- tmp_df
tmp_df2 <- melt(tmp_df2, 
               id.vars = 'batch', variable.name = 'class' , value.name = 'Age')
tmp_df2$class <- factor(tmp_df2$class, levels=c('true','pred', 'pred_adj'))
ggplot(tmp_df2[!tmp_df2$class %in% c('pred'),], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0.05, size = 13.75/.pt)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'Healthy cohort')+
theme(plot.title = element_text(hjust= 0.5))
ggsave('TN_healthy_boxplot.svg', width = 2.5, height = 4)

options(repr.plot.width=2.5, repr.plot.height=4)
tmp_df <- ind_plot_df[rownames(ind_plot_df) %in% c(centenarians_cohort),]
tmp_df$pred_adj <- tmp_df$pred - predict(loess_m, tmp_df$true)
tmp_df$batch <- rownames(tmp_df)
tmp_df$source <- NULL
ggplot(tmp_df, 
       aes(true, pred_adj))+
geom_point()+
ggpubr::stat_cor()+
stat_smooth(method = 'lm')
tmp_df$age_diff <- NULL
tmp_df2 <- tmp_df
tmp_df2 <- melt(tmp_df2, 
               id.vars = 'batch', variable.name = 'class' , value.name = 'Age')
tmp_df2$class <- factor(tmp_df2$class, levels=c('true','pred', 'pred_adj'))
ggplot(tmp_df2[!tmp_df2$class %in% c('pred'),], aes(class, Age,fill=class))+
geom_boxplot()+
geom_point(mapping = aes())+geom_line(mapping = aes(group = batch)) +
ggpubr::stat_compare_means(method = 't.test',label.x.npc = 0, size = 13.75/.pt, hjust = 0.25)+
d_theme_w(size=13.75, italic_text = '')+
scale_x_discrete(label = function(x){
    x <-gsub('true', 'Actual', x)
    x <-gsub('pred_adj', 'Predicted', x)
})+
guides(fill = F)+
scale_fill_manual(values = c('true'='#f6a000', 'pred_adj'='#017cf2'))+
labs(x =NULL, title = 'Centenarian cohort')+
theme(plot.title = element_text(hjust= 0.5))+
scale_y_continuous(limits = c(70,112))
ggsave('TN_centenarians_boxplot.svg', width = 2.5, height = 4)

save(final_model, loess_m, file = '250404_final_model.RData')

tmp_df <- ind_plot_df[rownames(ind_plot_df) %in% c(centenarians_cohort),]
tmp_df$pred_adj <- tmp_df$pred - predict(loess_m, tmp_df$true)
residuals <- tmp_df$pred_adj-tmp_df$true

cor_residuals <- psych::corr.test(residuals,
                                  TN_exp_mean[rownames(tmp_df),])

map_df <- data.frame(ensembl = colnames(TN_exp_mean))
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

r <- t(cor_residuals$r);r[is.na(r)] <- 0
p.adj <- t(cor_residuals$p.adj);p.adj[is.na(p.adj)] <- 1
p <- t(cor_residuals$p);p[is.na(p)] <- 1

cor_residuals_df <- data.frame(r , p,  p.adj)
colnames(cor_residuals_df) <- c('r', 'p', 'p.adj')
cor_residuals_df$ensembl <- rownames(cor_residuals_df)
cor_residuals_df <- merge(cor_residuals_df, map_df)

cor_residuals_df <- cor_residuals_df[order(cor_residuals_df$r, decreasing = T), ]
cor_residuals_df$ensembl <- factor(cor_residuals_df$ensembl, levels = rev(unique(cor_residuals_df$ensembl)))
cor_residuals_df$class <- ifelse(cor_residuals_df$r > 0, 'yes', 'no')
cor_residuals_df$class <- factor(cor_residuals_df$class, levels = c('yes','no'))
saveRDS(cor_residuals_df, file = 'cor_residuals_df.rds')
write.xlsx(cor_residuals_df,file ='cor_residuals_df.xlsx', row.names=T)

####GSEA
pcc <- cor_residuals_df$r
names(pcc) <- cor_residuals_df$symbol
pcc <- sort(pcc, decreasing = T)

write.table(data.frame(names(pcc), pcc), file = 'centenarian_pcc.rnk', quote = F,sep = '\t',
            col.names = F, row.names = F)

#Linux environment
gsea-cli.sh GSEAPreranked -rnk centenarian_pcc.rnk \
-gmx /data1/02.private/dengyj/analysis/thymus/Integrated/thymus_PB/Enrichment/GO_KEGG_REACTOME.gmt \
-collapse No_Collapse \
-set_max 10000 \
-rpt_label centenarian_pcc







save(train_plot_df, ind_plot_df, healthy_ind_plot_df, file = 'TN_aging_clock_plot_df.RData')


