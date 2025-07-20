library(Seurat)
library(Matrix)
library(openxlsx)
source('my_nonparametric_test.r')

PB_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/PB/combined_result/PB_step_2.rds')


SM <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/SM/humanSM.rds')
TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')

get_Age_Gene_test <- function(
    data,
    idents_var, 
    age_var,
    clusters, 
    file = 'Age_Gene_test_list.rds', 
    xlsx_file = 'Age_Gene_test_list.xlsx', 
    ...
    
){
    if(file.exists(file)){
        Age_Gene_test_list <- readRDS(file)
    }
    else{
        Age_Gene_test_list <- list()
    }
    exists_clu <- intersect(names(Age_Gene_test_list), clusters)
    clusters <- unique(setdiff(clusters, names(Age_Gene_test_list)))
    if(length(exists_clu) > 0){
        print(paste0(paste(exists_clu, collapse = ', '), ' had been calculated, ignored these clusters'))
        print(paste0('Caculated the rest clusters: ', paste(clusters, collapse = ', ')))        
    }

    new_Age_Gene_test_list <- lapply(clusters, function(cluster){
        data_cluster <- data[, idents_var %in%  cluster]
        idents_var_cluster <- idents_var[idents_var %in% cluster]
        age_var_cluster <- age_var[idents_var %in% cluster]
        result_cluster <- my_nonparametric_test(data = data_cluster, group_var = age_var_cluster,...)
        result_cluster$class <- 'other'
        result_cluster$class[result_cluster[, grep('variable', colnames(result_cluster))[1]] %in% TF] <- 
        'transcription factor'
        result_cluster$class[result_cluster[, grep('variable', colnames(result_cluster))[1]] %in% SM] <- 
        'surface marker'
        return(result_cluster)
    })
    names(new_Age_Gene_test_list) <- clusters
    Age_Gene_test_list <- c(Age_Gene_test_list, new_Age_Gene_test_list)
    
    
    saveRDS(Age_Gene_test_list, file = file)
    
    wb <- createWorkbook()
    for(i in names(Age_Gene_test_list)){
        addWorksheet(wb, sheetName = i)
        writeData(wb, sheet = i, rowNames = F, x = Age_Gene_test_list[[i]])
    }
    saveWorkbook(wb, file = xlsx_file, overwrite = TRUE) 
    
    return(Age_Gene_test_list)
}

PB_step_2$new_idents <- as.character(PB_step_2$identity)
####combined naive and RTE because of the lack of RTE in centenarian samples
PB_step_2$new_idents[grepl('CD8.*TN|CD8.*RTE',PB_step_2$new_idents)] <- 'CD8+TN'
PB_step_2$new_idents[grepl('CD4.*TN|CD4.*RTE',PB_step_2$new_idents)] <- 'CD4+TN'
PB_step_2_identity <- as.character(PB_step_2$new_idents)
PB_step_2_age <- PB_step_2$age_3

Age_Gene_test_list <- 
get_Age_Gene_test(data = expm1(PB_step_2$RNA@data), clusters = unique(PB_step_2_identity), combined_rest_group=T,
                  idents_var = PB_step_2_identity, age_var = PB_step_2_age)
