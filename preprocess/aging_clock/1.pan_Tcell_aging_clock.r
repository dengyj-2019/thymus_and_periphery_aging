library(SCopeLoomR)
library(Seurat)
library(ggplot2)
library(Matrix)
library(harmony)
library(plyr)
library(future)
library(pheatmap)
library(patchwork)


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

options(future.globals.maxSize = 100 * 1024^3)
plan('multicore', workers = 4)



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



CD8T_2023_Imm <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_Immunity_aging_PBMC/CD8T.rds')

CD8T_2023_Imm$identity <- as.character(CD8T_2023_Imm$Cluster_names)
CD8T_2023_Imm$source <- '2023_Imm'
CD8T_2023_Imm$batch <- paste0(CD8T_2023_Imm$Donor_id, '_', CD8T_2023_Imm$Age)
# CD8T_2023_Imm$batch2 <- CD8T_2023_Imm$Batch

CD4T_2023_Imm <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_Immunity_aging_PBMC/CD4T.rds')

CD4T_2023_Imm$identity <- as.character(CD4T_2023_Imm$Cluster_names)
CD4T_2023_Imm$source <- '2023_Imm'
CD4T_2023_Imm$batch <- paste0(CD4T_2023_Imm$Donor_id, '_', CD4T_2023_Imm$Age)
# CD4T_2023_Imm$batch2 <- CD4T_2023_Imm$Batch

T_2023_Imm <- merge(CD4T_2023_Imm, CD8T_2023_Imm)

all(rownames(T_2023_Imm) %in% features_df_V6$symbol)

exp <- T_2023_Imm$RNA@counts

if(all(rownames(exp) %in% features_df_V6$symbol)){
   rownames(exp) <- features_df_V6$ensembl[match(rownames(exp), features_df_V6$symbol)] 
}

all(rownames(exp) %in% features_df_V6$ensembl)

T_2023_Imm <- CreateSeuratObject(counts = exp, meta.data = T_2023_Imm@meta.data)

T_2023_Imm <- NormalizeData(T_2023_Imm, verbose = F)





Tcell <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/T_NK/Tcell.rds')

Tcell$Age <- Tcell$age
Tcell$source <- 'thymus_aging'

T_this_study <- Tcell[, grepl('RTE|T',Idents(Tcell))]

all(rownames(T_this_study$RNA@counts) %in% combined_features_df$symbol)

exp <- T_this_study$RNA@counts

if(all(rownames(exp) %in% combined_features_df$symbol)){
   rownames(exp) <- combined_features_df$ensembl[match(rownames(exp), combined_features_df$symbol)] 
}

all(rownames(exp) %in% combined_features_df$ensembl)

T_this_study <- CreateSeuratObject(counts = exp, meta.data = T_this_study@meta.data)

T_this_study <- NormalizeData(T_this_study, verbose = F)





T_2022_NA <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2022_NA_aging_PBMC/Tcell.rds')

T_2022_NA$Age <- T_2022_NA$age
T_2022_NA$identity <- as.character(Idents(T_2022_NA))

T_2022_NA$source <- '2022_NA'

# T_2022_NA$batch2 <- T_2022_NA$source



exp <- T_2022_NA$RNA@counts

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

T_2022_NA <- CreateSeuratObject(counts = exp, meta.data = T_2022_NA@meta.data)

T_2022_NA <- NormalizeData(T_2022_NA, verbose = F)



CD4T_2024_AIDA_V1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2024_AIDA_V1/CD4T.rds')
CD4T_2024_AIDA_V1$batch <- CD4T_2024_AIDA_V1$donor_id
CD4T_2024_AIDA_V1$identity <- 'CD4T'
CD4T_2024_AIDA_V1$source <- '2024_AIDA_V1'
CD4T_2024_AIDA_V1$Age <- as.numeric(gsub('-year-old stage', '', CD4T_2024_AIDA_V1$development_stage))
# CD4T_2024_AIDA_V1$batch2 <- CD4T_2024_AIDA_V1$source

CD8T_2024_AIDA_V1 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2024_AIDA_V1/CD8T.rds')
CD8T_2024_AIDA_V1$batch <- CD8T_2024_AIDA_V1$donor_id
CD8T_2024_AIDA_V1$identity <- 'CD8T'
CD8T_2024_AIDA_V1$source <- '2024_AIDA_V1'
CD8T_2024_AIDA_V1$Age <- as.numeric(gsub('-year-old stage', '', CD8T_2024_AIDA_V1$development_stage))
# CD8T_2024_AIDA_V1$batch2 <- CD8T_2024_AIDA_V1$source

T_2024_AIDA_V1 <- merge(CD4T_2024_AIDA_V1, CD8T_2024_AIDA_V1)

T_2024_AIDA_V1 <- NormalizeData(T_2024_AIDA_V1, verbose = F)

T_2025_NI <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2025_NI_Hs_PBMC_aging/Tcell.rds')

T_2025_NI$identity <- as.character(T_2025_NI$secondary_type)
T_2025_NI$source <- '2025_NI'
T_2025_NI$batch <- T_2025_NI$sampleName
# T_2025_NI$batch2 <- T_2025_NI$source

all(rownames(T_2025_NI) %in% features_df_V3$symbol)

exp <- T_2025_NI$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

T_2025_NI <- CreateSeuratObject(counts = exp, meta.data = T_2025_NI@meta.data)

T_2025_NI <- NormalizeData(T_2025_NI, verbose = F)



T_2023_SA <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/2023_SA_aging_PBMC/2023_SA_aging_Tcell.rds')

T_2023_SA$identity <- 'Naive'
T_2023_SA$source <- '2023_SA'

T_2023_SA$batch2 <- T_2023_SA$source

all(rownames(T_2023_SA) %in% features_df_V3$symbol)

exp <- T_2023_SA$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

T_2023_SA <- CreateSeuratObject(counts = exp, meta.data = T_2023_SA@meta.data)

T_2023_SA <- NormalizeData(T_2023_SA, verbose = F)



T_GSE158055 <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE158055/Tcell.rds')


all(rownames(T_GSE158055) %in% features_df_V3$symbol)

exp <- T_GSE158055$RNA@counts

if(all(rownames(exp) %in% features_df_V3$symbol)){
   rownames(exp) <- features_df_V3$ensembl[match(rownames(exp), features_df_V3$symbol)] 
}

all(rownames(exp) %in% features_df_V3$ensembl)

T_GSE158055 <- CreateSeuratObject(counts = exp, 
                                               meta.data = T_GSE158055@meta.data)

T_GSE158055 <- NormalizeData(T_GSE158055, verbose = F)



T_GSE161918 <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/GSE161918//GSE161918_Tcell.rds')


exp <- T_GSE161918$RNA@counts

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

T_GSE161918 <- CreateSeuratObject(counts = exp, 
                                               meta.data = T_GSE161918@meta.data)

T_GSE161918 <- NormalizeData(T_GSE161918, verbose = F)





obj_list <- c('T_2023_Imm','T_2022_NA', 
'T_2024_AIDA_V1', 'T_2025_NI', 'T_2023_SA',
             'T_GSE158055','T_this_study', 'T_GSE161918')

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

# saveRDS(exp_mean_list, file = 'exp_mean_list.rds')

if(exists(x = 'T_exp_mean_orig')){
    T_exp_mean_orig <- T_exp_mean_orig
    logi <- sapply(exp_mean_list, function(x){
        all(colnames(x) %in% colnames(T_exp_mean_orig))
    })
    obj_names <- names(exp_mean_list)[!logi]
    if(length(obj_names) > 0){
        for(obj_name in obj_names){
            df_iter <- data.frame(exp_mean_list[[obj_name]], check.rows = F, check.names = F)
            T_exp_mean_orig <- merge(T_exp_mean_orig, df_iter, by = 'row.names', all = T)
            rownames(T_exp_mean_orig) <- T_exp_mean_orig$Row.names
            T_exp_mean_orig$Row.names <- NULL
            T_exp_mean_orig[is.na(T_exp_mean_orig)] <- 0
        }           
    }
 
}else{
    for(i in 1:length(exp_mean_list)){
        if(i==1){
            T_exp_mean_orig <- data.frame(exp_mean_list[[i]], check.rows = F, check.names = F)
        }else{
            df_iter <- data.frame(exp_mean_list[[i]], check.rows = F, check.names = F)
            T_exp_mean_orig <- merge(T_exp_mean_orig, df_iter, by = 'row.names', all = T)
            rownames(T_exp_mean_orig) <- T_exp_mean_orig$Row.names
            T_exp_mean_orig$Row.names <- NULL
            T_exp_mean_orig[is.na(T_exp_mean_orig)] <- 0
        }
    }
}
# saveRDS(T_exp_mean_orig, file = 'T_exp_mean_orig.rds')
T_exp_mean <- t(T_exp_mean_orig)

common_meta <- c('source', 'Age', 'batch')

meta_list <- lapply(names(exp_mean_list), function(obj){
    obj <- get(obj)
    meta <- obj[[common_meta]]
})
meta_list <- do.call(rbind, meta_list)




T_meta <- meta_list
T_meta <- T_meta[!duplicated(T_meta$batch), c('source','Age','batch')]
rownames(T_meta) <- T_meta$batch

# saveRDS(T_meta, file = 'T_meta.rds')

identical(rownames(T_exp_mean),rownames(T_meta))

samples_used <- unique(meta_list$batch)
samples_tbl <- table(meta_list$batch)
samples_used <- intersect(samples_used,
                          names(samples_tbl)[samples_tbl >= 50])

names(samples_tbl)[samples_tbl < 50]

modeling_samples <- samples_used
modeling_samples <- setdiff(modeling_samples, 
                            T_meta$batch[T_meta$source=='2019_PNAs' & T_meta$Age < 100])




source_for_gene_detected <- c('2023_Imm', '2024_AIDA_V1', '2025_NI')
if(exists(x = 'lm_list')){
    lm_list <- lm_list
}else{
    lm_list <- list()
}
for(source in setdiff(source_for_gene_detected, names(lm_list))){
    logi <- T_meta$source == source & T_meta$batch %in% modeling_samples
    tmp_exp <- T_exp_mean[logi, ]
    tmp_exp <- tmp_exp[, colSums(tmp_exp) > 0]
    data_used <- cbind(T_meta[logi, ], tmp_exp)    
    lm_res <- future_lapply(colnames(tmp_exp), function(gene){
        formula <- as.formula(paste0('`', gene, '`', " ~ Age"))
        # 建立线性回归模型
        model <- lm(formula, data = data_used)
        pval <- summary(model)$coefficients[, 4][2]
        slope = coef(model)[2]
        #print(which(colnames(T_exp_pseudobulk_orig)==gene))
        df <- data.frame(pval = pval, gene = gene, slope = slope)         
    })
    lm_res <- do.call(rbind, lm_res)
    lm_res$source <- source
    lm_res$padj <- p.adjust(lm_res$pval, method = 'fdr')
    lm_list[[source]] <-   lm_res  
}



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



source_for_modeling <- c('2024_AIDA_V1',  '2025_NI')

unique(T_meta$source)

frail_cohort <- unique(T_2022_NA$batch[T_2022_NA$age_group=='frail'])
#covid19_cohort1 <- T_meta$batch[T_meta$source %in% c('2022_NC_COVID19_aging')]
centenarians_cohort <- T_meta$batch[T_meta$Age >=100]
covid19_cohort1 <- T_GSE158055$batch[!T_GSE158055$CoVID.19.severity %in% c('control')]
covid19_cohort2 <- T_GSE161918$batch[!T_GSE161918$Class %in% c('HC')]
T_2023_Imm_batch <- unique(T_2023_Imm$batch)

ind_exclude_sample <- union(c(frail_cohort,covid19_cohort1,centenarians_cohort,covid19_cohort2), 
c('S-HC013','S-HC014',T_2023_Imm_batch))

save(frail_cohort, centenarians_cohort, covid19_cohort1, covid19_cohort2,T_2023_Imm_batch, 
    ind_exclude_sample, file = 'T_aging_clock_samples_used.RData')


logi <- T_meta$source %in% source_for_modeling#
batch_for_modeling <- intersect(unique(T_meta$batch[logi]), modeling_samples)
batch_for_modeling <- setdiff(batch_for_modeling, centenarians_cohort)
# batch_for_modeling <- setdiff(batch_for_modeling, covid19_cohort2)
set.seed(1234)
train_samples <- sample(batch_for_modeling, size = 0.7 * length(batch_for_modeling))
test_samples <- setdiff(batch_for_modeling, train_samples)

logi <- T_meta$batch %in% train_samples
train_x <- log1p(T_exp_mean[logi, gene_used])
train_y <- T_meta$Age[logi]

logi <- T_meta$batch %in% test_samples
test_x <- log1p(T_exp_mean[logi, gene_used])
test_y <- T_meta$Age[logi]



ind_samples <- setdiff(modeling_samples, batch_for_modeling)

logi <- T_meta$batch %in% ind_samples
ind_x <- log1p(T_exp_mean[logi, gene_used])
ind_y <- T_meta$Age[logi]




alpha <- 0.5
seed <- 1
set.seed(seed)
cv_model <- cv.glmnet(x = train_x, y = train_y, alpha = alpha)
final_model <- glmnet(train_x, train_y, alpha = alpha, lambda = cv_model$lambda.min)
train_plot_df <- data.frame(predict(final_model, newx = train_x), train_y)
colnames(train_plot_df) <- c('pred', 'true')    
test_plot_df <- data.frame(predict(final_model, newx = test_x), test_y)
colnames(test_plot_df) <- c('pred', 'true')
#AIDA_plot_df <- data.frame(predict(final_model, newx = AIDA_x), AIDA_y)
#colnames(AIDA_plot_df) <- c('pred', 'true')  
ind_plot_df <- data.frame(predict(final_model, newx = ind_x), ind_y)
colnames(ind_plot_df) <- c('pred', 'true')  

# logi <- !rownames(ind_plot_df) %in% c(#frail_cohort, 
#                                          covid19_cohort1,centenarians_cohort,
#                                          covid19_cohort2,covid19_cohort3 ) & 
#                !rownames(ind_plot_df) %in% c('S-HC013','S-HC014')
#ind_plot_df <- ind_plot_df[logi,]
train_val <- mean(abs(train_plot_df$pred-train_plot_df$true))
test_val <- mean(abs(test_plot_df$pred-test_plot_df$true))
#AIDA_val <- mean(abs(AIDA_plot_df$pred-AIDA_plot_df$true)) 
ind_val <- mean(abs(ind_plot_df$pred-ind_plot_df$true))   

train_cor <- cor(train_plot_df$pred,train_plot_df$true)
test_cor <- cor(test_plot_df$pred,test_plot_df$true)
#AIDA_cor <- cor(AIDA_plot_df$pred,AIDA_plot_df$true)
ind_cor <- cor(ind_plot_df$pred,ind_plot_df$true)
train_plot_df$source <- T_meta$source[match(rownames(train_plot_df), T_meta$batch)]
ind_plot_df$source <- T_meta$source[match(rownames(ind_plot_df), T_meta$batch)]
df <- data.frame(train_val = train_val, test_val = test_val, 
                 #AIDA_val = AIDA_val, 
                 ind_val= ind_val, 
                 train_cor = train_cor, test_cor = test_cor, 
                 #AIDA_cor = AIDA_cor,
                 ind_cor =ind_cor,
                 seed = seed)

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

options(repr.plot.width=5, repr.plot.height=5)
ggplot(train_plot_df, 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5))+
labs(x = 'Actual age', y='Predicted age', title = 'Training cohort')
ggsave('T_aging_clock_train_scatter.svg', width = 5, height = 5)

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

options(repr.plot.width=5, repr.plot.height=5)
ggplot(healthy_ind_plot_df[ ], 
       aes(true, pred_adj))+
geom_point(mapping = aes())+
ggpubr::stat_cor(size = 13.75/.pt)+
stat_smooth(method = 'lm')+
d_theme_w(size =13.75,italic_text = '')+
theme(plot.title = element_text(hjust=0.5),
     plot.margin = margin(r = 2.4, unit = 'mm'))+
labs(x = 'Actual age', y='Predicted age', title = 'Healthy cohort')
ggsave('T_aging_clock_test_scatter.svg', width = 5, height = 5)

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

save(train_plot_df, ind_plot_df,healthy_ind_plot_df, file = 'T_aging_clock_plot_df.RData')

saveRDS(exp_mean_list, file = 'T_exp_mean_list.rds')
saveRDS(T_exp_mean_orig, file = 'T_exp_mean_orig.rds')
saveRDS(T_meta, file = 'T_meta.rds')
saveRDS(lm_list, file = 'T_lm_list.rds')

saveRDS(meta_list, file = 'T_meta_list.rds')

combined_features_df[grepl('ENSG00000113302', combined_features_df$ensembl),]




