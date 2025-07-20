library(Seurat)
library(openxlsx)
library(future)
library(future.apply)
library(Matrix)
library(openxlsx)

#################################################################正式计算
setwd('/data1/02.private/dengyj/analysis/thymus/crosstalk/2502222/')###working dir
source('support/calculation.R')###code for calculation
load('support/Human_crosstalk_DB_20_12_20.RData')###ligand-receptor database
complex_names <- Hs_complex_info_V2$complex_name
Hs_complex_info_V2 <- lapply(complex_names, function(complex){
  component <- as.character(Hs_complex_info_V2[Hs_complex_info_V2$complex_name %in% complex,
                                         grep('component', colnames(Hs_complex_info_V2))])
  component <- component[component != 'None']
})
names(Hs_complex_info_V2) <- complex_names



filtered_combined_Niche_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_combined_Niche_2.rds')


ETP_TP <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/ETP_TP.rds')

filtered_combined_Niche_2[['RNA']] <- filtered_combined_Niche_2[['symbol']]
DefaultAssay(filtered_combined_Niche_2) <- 'RNA'
filtered_combined_Niche_2[['symbol']] <- NULL

combined_thymocytes <-merge(filtered_combined_Niche_2, ETP_TP)

combined_thymocytes$identity <- Idents(combined_thymocytes)

#saveRDS(combined_thymocytes, file = 'combined_thymocytes.rds')

crosstalk <- pre_prossessing(expm1(combined_thymocytes$RNA@data), database = Hs_DB_V2,
                             identity = combined_thymocytes$identity, condition = combined_thymocytes$Age_group, 
                             complex_info = Hs_complex_info_V2)




plan('multicore', workers = 4)

options(future.globals.maxSize = 100 * 1024^3)
crosstalk <- get_permutation2(crosstalk,  multiprocess=T,iteration = 1000,  block_size = 50, 
                              slot_used = 'condition', 
                              ligand_pct_threshold = 0.2, receptor_pct_threshold = 0.1, style = 'product')

plan('sequential')

saveRDS(crosstalk, file ='crosstalk.rds')