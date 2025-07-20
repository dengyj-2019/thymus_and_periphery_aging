library(Seurat)
library(Matrix)

library(openxlsx)

library(ggplot2)
library(plyr)
library(reshape2)
source('my_nonparametric_test.r')


filtered_combined_Niche_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_Niche_2.rds')

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
        result_cluster <- my_nonparametric_test(data = data_cluster, group_var = age_var_cluster, ...)
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

filtered_combined_Niche_2_identity <- as.character(Idents(filtered_combined_Niche_2))

filtered_combined_Niche_2_age <- filtered_combined_Niche_2$Age_group


Age_Gene_test_list_new <- 
get_Age_Gene_test(data = expm1(filtered_combined_Niche_2$RNA@data), clusters = unique(filtered_combined_Niche_2_identity),
                  file = 'Age_Gene_test_list_combined.rds', 
                  xlsx_file = 'Age_Gene_test_list_combined.xlsx', 
                  combined_rest_group=T,
                  idents_var = filtered_combined_Niche_2_identity, age_var = filtered_combined_Niche_2_age)
