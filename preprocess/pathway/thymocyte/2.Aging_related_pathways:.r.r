library(Seurat)
library(Matrix)
library(openxlsx)
source('my_nonparametric_test.r')

thymus_step_2 <- 
readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/combined_result/thymus_step_2.rds')



go <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/c5.all.v7.0.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(go), function(x) {go[[x]][1]})
go <- sapply(1:length(go), function(x) {go[[x]][-c(1,2)]})
names(go) <- geneset_names

reactome <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/ReactomePathways.gmt"), split = '\t')
geneset_names <- sapply(1:length(reactome), function(x) {reactome[[x]][1]})
reactome <- sapply(1:length(reactome), function(x) {reactome[[x]][-c(1,2)]})
names(reactome) <- geneset_names

kegg <- strsplit(readLines("/data1/02.private/dengyj/analysis/database/GSEA/c2.cp.kegg.v7.0.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(kegg), function(x) {kegg[[x]][1]})
kegg <- sapply(1:length(kegg), function(x) {kegg[[x]][-c(1,2)]})
names(kegg) <- geneset_names

database <- c(go, kegg, reactome)


thymus_step_2_es <- readRDS('/data1/02.private/dengyj/analysis/thymus/Enrichment/ssGSEA/thymus/thymus_step_2_es.rds')

get_Age_GS_test <- function(
    data,
    idents_var, 
    age_var,
    clusters, 
    file = 'Age_GS_test_list.rds', 
    xlsx_file = 'Age_GS_test_list.xlsx',
    database, 
    genes, 
    ...
){
    if(file.exists(file)){
        Age_GS_test_list <- readRDS(file)
    }
    else{
        Age_GS_test_list <- list()
    }
    exists_clu <- intersect(names(Age_GS_test_list), clusters)
    clusters <- unique(setdiff(clusters, names(Age_GS_test_list)))
    if(length(exists_clu) > 0){
        print(paste0(paste(exists_clu, collapse = ', '), ' had been calculated, ignored these clusters'))
        print(paste0('Caculated the rest clusters: ', paste(clusters, collapse = ', ')))        
    }

    new_Age_GS_test_list <- lapply(clusters, function(cluster){
        data_cluster <- data[, idents_var %in%  cluster]
        idents_var_cluster <- idents_var[idents_var %in% cluster]
        age_var_cluster <- age_var[idents_var %in% cluster]
        result_cluster <- my_nonparametric_test(data = data_cluster, group_var = age_var_cluster, ...)
        result_cluster$genes <- get_pathway_gene(
            pathways = result_cluster[, grep('variable$', colnames(result_cluster))[1]], 
            database = database,
            genes = genes
        )
        return(result_cluster)
    })
    names(new_Age_GS_test_list) <- clusters
    Age_GS_test_list <- c(Age_GS_test_list, new_Age_GS_test_list)
    
    
    saveRDS(Age_GS_test_list, file = file)
    
    wb <- createWorkbook()
    for(i in names(Age_GS_test_list)){
        addWorksheet(wb, sheetName = i)
        writeData(wb, sheet = i, rowNames = F, x = Age_GS_test_list[[i]])
    }
    saveWorkbook(wb, file = xlsx_file, overwrite = TRUE) 
    
    return(Age_GS_test_list)
}

thymus_step_2_identity <- as.character(thymus_step_2$identity)
###combined ETP and TP cells due to the limited numbers during aging
thymus_step_2_identity[grep('ETP|TP', thymus_step_2_identity)] <- 'ETP_TP' 
thymus_step_2_age <- thymus_step_2$age_3


Age_GS_test_list <- 
get_Age_GS_test(data = thymus_step_2_es, clusters = unique(thymus_step_2_identity), database = database, 
                genes = rownames(thymus_step_2$RNA@data)[rowSums(thymus_step_2$RNA@data) > 0],
                combined_rest_group = T,
                idents_var = thymus_step_2_identity, age_var = thymus_step_2_age)

