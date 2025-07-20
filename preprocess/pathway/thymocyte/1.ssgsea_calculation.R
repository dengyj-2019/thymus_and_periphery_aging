library(Seurat)
library(future)
library(openxlsx)
library(Matrix)
source('fast_ssgsea.R')


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

thymus_step_2 <- readRDS('/data1/02.private/dengyj/analysis/thymus/clustering_10_20/combined_result/thymus_step_2.rds')

exp <- thymus_step_2$RNA@data[rowSums(thymus_step_2$RNA@data) > 0, ]


options(future.globals.maxSize = 64 * 1024^3)
plan('multisession', workers = 32)
thymus_step_2_es <- ssgsea(exp = exp , genesets = database, 
                      min_size = 10, max_size = 500, exp_filter = F, normalization = F)
plan('sequential')

saveRDS(thymus_step_2_es, file = 'thymus_step_2_es.rds')

