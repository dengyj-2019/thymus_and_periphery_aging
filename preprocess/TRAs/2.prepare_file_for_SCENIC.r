library(Seurat)
library(SCopeLoomR)

filtered_TRA_TEC <- readRDS('/data1/02.private/dengyj/analysis/thymus/Integrated/Niche/filtered_TRA_TEC.rds')


loom <- build_loom('filtered_TRA_TEC.loom',dgem=filtered_TRA_TEC$symbol@data)
close_loom(loom)
