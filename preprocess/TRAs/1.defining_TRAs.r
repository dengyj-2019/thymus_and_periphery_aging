################prepare files for TRAs
library(openxlsx)
library(data.table)
library(future)
load('/data1/02.private/dengyj/analysis/database/specie_conversion/specie_conversion.rds')
setwd('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1')

cal_tau <- function(val){
  nor_val <- 1-val/max(val)
  tau <- sum(nor_val)/(length(nor_val)-1)
  return(tau)
}



features_df <- read.table('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1//features.tsv')[, c('V1', 'V2')]
colnames(features_df) <- c('ensembl', 'symbol')

data_consensus <- data.frame(fread('/data1/02.private/dengyj/analysis/thymus/TRAs_HPA_V21.1/consensus/rna_tissue_consensus.tsv'))

data_consensus$nTPM[data_consensus$nTPM < 1] <- 0
data_consensus$nTPM <- log10(data_consensus$nTPM+1)
tau_consensus <- pbapply::pbsapply(unique(data_consensus$Gene), function(i){
  val3 <- data_consensus[data_consensus$Gene == i, 'nTPM']
  tau3 <- cal_tau(val3)
  return(c(tau3))
})

options(repr.plot.width=6, repr.plot.height=6)
cutoff <- 0.8
plot(density(tau_consensus[!is.na(tau_consensus)]))
abline(v = cutoff, col = 'red')
sum(tau_consensus[!is.na(tau_consensus)] >= cutoff)
TRAs_consensus <- unique(features_df[features_df$ensembl %in% 
                                       names(tau_consensus[!is.na(tau_consensus)])[tau_consensus[!is.na(tau_consensus)] >= cutoff], 'symbol'])

tissue_consensus <- future.apply::future_sapply(unique(data_consensus$Gene), function(x){
  df_iter <- data_consensus[data_consensus$Gene %in% x, ]
  return(df_iter[which.max(df_iter$nTPM), 'Tissue'])
})

TRAs_consensus_df <- data.frame('ensembl' = names(tau_consensus), 
                           'symbol' =  features_df[match(names(tau_consensus), features_df$ensembl) , 'symbol'], 
                           tau = tau_consensus, tissue = tissue_consensus)





save(TRAs_consensus,file = 'TRAs_symbol.RData')
save(tau_consensus, file = 'TRAs_tau.RData')
save(tissue_consensus, file = 'TRAs_tissue.RData')
save(TRAs_consensus_df, file = 'TRAs_df.RData')



