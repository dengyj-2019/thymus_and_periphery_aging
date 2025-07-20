library(openxlsx)
library(future)
library(future.apply)
library(Matrix)
library(stringr)
library(psych)
library(AUCell)
library(pbapply)

setClass("CrossTalk",slots=list(data="dgCMatrix",nor_data = 'dgCMatrix', LR_data = 'list',#grn_data = 'dgCMatrix',
                                database="data.frame",complex_info = 'list', TF_target = "data.frame", 
                                regulon = 'list',GS_data = 'matrix',
                                identity = 'factor', condition = 'factor', 
                                p_adj_gene_idents = 'list', p_adj_gene_con = 'list',
                                p_adj_GS_idents = 'list', p_adj_GS_con = 'list',
                                inter_intensity_idents = 'list', inter_intensity_con = 'list',
                                pct_idents = 'list', pct_con = 'list', 
                                intra_intensity_idents = 'list', intra_intensity_con = 'list',
                                intensity_idents = 'list', intensity_con = 'list', 
                                activation_mat_idents = 'list', activation_mat_con = 'list', 
                                permutation_idents = 'matrix',permutation_con = 'matrix',
                                #pboot_idents = 'list', pboot_con = 'list', 
                                sig_db_idents = 'data.frame', sig_db_con = 'data.frame', 
                                condition_net = 'data.frame',
                                pnull_idents = 'list', pnull_con = 'list', 
                                L_R_matrix_normalization = 'logical', 
                                inter_p_threshold = 'numeric', intra_p_threshold = 'numeric'))

setMethod(
  f = "show",
  signature = "CrossTalk",
  definition = function(object) {
    cat("An object of class 'Crosstalk'", "\n")
    cat(nrow(object@LR_data$ligand_data), 'ligands and', nrow(object@LR_data$receptor_data),
        'receptors across', ncol(object@data), 'samples')
  }
)

pre_prossessing <- function(
  data,##########需要输入已经进行了normalization的数据
  database, 
  TF_target = NULL, 
  identity = NULL, 
  condition = NULL,
  complex_info, 
  interaction_used=NULL, #############传入自定义配受体互作
  perform_normalize = T, 
  inter_p_threshold = 0.05#, 
  #build_grn_data = F###########不推荐用TRUE，因为过于占用内存
){
  if(class(data) != 'dgCMatrix'){
    if(class(data) == 'matrix'){
      data <- Matrix::Matrix(data, sparse = T)
    } else{
      stop('Data must be a "matrix" or "dgCMatrix" object')
    }
  }
  data <- data[rowSums(data) > 0, ]###去除所有细胞都不表达的基因
  
  intercelllular_colnames <- c('interaction_name',  'ligand',  'receptor')
  intracelllular_colnames <- c('downstream_TF')
  ###检测是否database具有合适的列名
  if(any(!intercelllular_colnames %in% colnames(database))){
    stop(paste0(intercelllular_colnames[!intercelllular_colnames %in% colnames(database)], 
                ' neccessary for intercellular signaling analysis was/were not in the column names of database, please check!'))
  }
  if(any(!intracelllular_colnames %in% colnames(database))){
    warning(paste0(intracelllular_colnames[!intracelllular_colnames %in% colnames(database)], 
                   ' neccessary for intracellular signaling analysis was/were not in the column names of database, please check!'))
  }
  #####去除重复的interaction
  database <- database[!duplicated(database$interaction_name), ]
  #####这一步防止complex和基因的名字重复导致后面的计算混乱，所以会去除与基因名重复的complex
  if(any(names(complex_info) %in% rownames(data))){
    warning(paste0('Rename elements in complex_info whose names also exit in rownames of data: ', 
                   paste(names(complex_info)[names(complex_info) %in% rownames(data)], collapse = ', ')))
    names(complex_info)[names(complex_info) %in% rownames(data)] <- 
      paste0('complex_', names(complex_info)[names(complex_info) %in% rownames(data)])
  }
  #####去除重读的complex
  complex_info <- complex_info[!duplicated(names(complex_info))]
  #####要求复合物至少含两个或以上组件
  complex_info <- complex_info[sapply(complex_info, length) > 1]
  
  ####过滤database
  db_to_use <- get_db_to_use(data =data, database = database, complex_info = complex_info, interaction_used = interaction_used)
  if(nrow(db_to_use) == 0){
    stop('The database with row 0 after prossesing, please check!')
  }
  ####仅保留database中配受体的表达矩阵
  data_to_use <- get_gene_to_use(data, db_to_use, complex_info = complex_info)
  if(nrow(data_to_use) == 0 | ncol(data_to_use) == 0){
    stop('The data was empty after prossesing, please check!')
  }
  ####
  if(!is.null(TF_target)){
    ####要求TF_target具有合适的列名
    TF_target_colnames <- c('TF', 'target')
    if(any(!TF_target_colnames %in% colnames(TF_target))){
      stop(paste0(TF_target_colnames[!TF_target_colnames %in% colnames(TF_target)], 
                  'neccessary for intracellular signaling analysis was/were not in the column names of TF_target, please check!'))
    }
    ####仅保留在表达矩阵中存在的TF和target
    TF_target <- TF_target[TF_target$TF %in% rownames(data) & TF_target$target %in% rownames(data), ]
    if(nrow(TF_target) == 0){
      warning('The TF_target slot necessary for intracellular signaling analysis with row 0 after filtering!')
    }
  }else{
    TF_target <- data.frame('TF' = NA, 'target'=NA)
  }
  object <- new("CrossTalk",
                data=data_to_use,database=db_to_use, complex_info = complex_info, 
                TF_target = TF_target)
  # if(build_grn_data){##############
  #   object@grn_data <- data[rownames(data) %in% union(TF_target$TF, TF_target$target), ]
  # }
  ####检测identity slot且identity slot是必须提供的
  if(!is.null(identity)){
    if(length(identity) != ncol(data)){
      stop('The length of identity must be identical with number of cells in data!')
    }
    if(class(identity) != 'factor'){
      identity <- as.factor(identity)
    }
    identity <- droplevels(identity)
    if(any(str_detect(levels(identity), '\\||~'))){
      stop('Special strings ("~" or "|") were detected in identity, please change the these strings for further analysis!')
    }
    object@identity <- identity
  }else{
    stop('Please provide the identity slot that is necessary for downstream analysis!')
  }
  ####检测condition slot,  condition slot不是必须提供的
  if(!is.null(condition)){
    if(length(condition) != ncol(data)){
      stop('The length of condition must be identical with number of cells in data!')
    }
    if(class(condition) != 'factor'){
      condition <- as.factor(condition)
    }
    condition <- droplevels(condition)
    if(any(str_detect(levels(condition), '\\||~'))){
      stop('Special strings ("~" or "|") were detected in condition, please change the these strings for further analysis!')
    }
    object@condition <- condition
  }else{
    warning('Please provide the condition slot that is needed for downstream analysis!')
  }
  ####构建对象的normalized data，进行（min-max normalization)
  mat <- t(apply(object@data, 1, function(x) x/(max(x)-min(x))))
  mat <- mat
  object@nor_data <- Matrix::Matrix(mat, sparse = T)
  object <- build_LR_matrix(object= object, normalize = perform_normalize)
  ####设置用于下游构建regulon的表面信号互作阈值
  object@inter_p_threshold <- inter_p_threshold
  return(object)
}

get_db_to_use <- function(
  data,###行是基因，列是细胞
  database,
  complex_info,
  interaction_used = NULL
){
  if(!is.null(interaction_used)){
    db_logic <- database$interaction_name %in% interaction_used
    if(any(db_logic)){
      database <- database[db_logic, ]
    }else{
      stop('The selected interactions passed by "interaction_used" parameter does not exit in database!')
    }
  }
  keep <- 
    sapply(1:nrow(database), FUN = 
             function(x){
               if(database$ligand[x] %in% names(complex_info)){
                 component <- complex_info[[database$ligand[x]]]
                 if(all(component %in% rownames(data))){#这一步要求complex的所有component都在表达矩阵中
                   ligand_keep <- TRUE
                 } else{
                   ligand_keep <- FALSE
                 }
               } else if(database$ligand[x] %in% rownames(data)){
                 ligand_keep <- TRUE
               } else{
                 ligand_keep <- FALSE
               }
               
               if(database$receptor[x] %in% names(complex_info)){
                 component <- complex_info[[database$receptor[x]]]
                 if(all(component %in% rownames(data))){#这一步要求complex的所有component都在表达矩阵中
                   receptor_keep <- TRUE
                 } else{
                   receptor_keep <- FALSE
                 }
               } else if(database$receptor[x] %in% rownames(data)){
                 receptor_keep <- TRUE
               } else{
                 receptor_keep <- FALSE  
               }
               
               return(ligand_keep & receptor_keep)
             })
  database_to_use <- database[keep, ]
  database_to_use$ligand_class <- database_to_use$receptor_class <- 'simple'
  database_to_use$ligand_class[database_to_use$ligand %in% names(complex_info)] <- 'complex'
  database_to_use$receptor_class[database_to_use$receptor %in% names(complex_info)] <- 'complex'
  return(database_to_use)
}




get_gene_to_use <- function(
  data, 
  database, 
  complex_info
){
  gene_to_use <- 
    sapply(1:nrow(database), FUN = 
             function(x){
               ligand_gene <- 
                 if(database$ligand_class[x] %in% 'simple'){
                   database$ligand[x]
                 } else{
                   complex_info[[database$ligand[x]]]
                 }
               receptor_gene <- 
                 if(database$receptor_class[x] %in% 'simple'){
                   database$receptor[x]
                 } else{
                   complex_info[[database$receptor[x]]]
                 }
               return(c(ligand_gene, receptor_gene))
             }
    )
  data_to_use <- data[rownames(data) %in% unlist(gene_to_use),]
  return(data_to_use)
}


build_LR_matrix <- function(
  object,
  normalize = T
){
  if(normalize){
    data <- object@nor_data
    object@L_R_matrix_normalization <- TRUE
  }else{
    data <- object@data
    object@L_R_matrix_normalization <- FALSE
  }
  database <- object@database
  complex_info <- object@complex_info
  #######################################################
  ligand <- database$ligand
  names(ligand) <- database$ligand_class
  ligand <- ligand[!duplicated(ligand)]
  ligand_data <- 
    sapply(1:length(ligand), FUN = 
             function(x){
               ligand_iter <- ligand[x]
               if(names(ligand_iter) == 'simple'){
                 L_data <- data[ligand_iter, ]
               }else{
                 L_data <- 
                   get_complex_exp(
                     data = data, 
                     complex = ligand_iter, 
                     complex_info = complex_info)
               }
               return(L_data)
             })
  ligand_data <- t(ligand_data)
  rownames(ligand_data) <- unname(ligand)
  ligand_data <- Matrix::Matrix(ligand_data, sparse = T)
  #########################################################  
  receptor <- database$receptor
  names(receptor) <- database$receptor_class
  receptor <- receptor[!duplicated(receptor)]
  receptor_data <- 
    sapply(1:length(receptor), FUN = 
             function(x){
               receptor_iter <- receptor[x]
               if(names(receptor_iter) == 'simple'){
                 R_data <- data[receptor_iter, ]
               }else{
                 R_data <- 
                   get_complex_exp(
                     data = data, 
                     complex = receptor_iter, 
                     complex_info = complex_info)
               }
               return(R_data)
             })
  receptor_data <- t(receptor_data)
  rownames(receptor_data) <- unname(receptor)
  receptor_data <- Matrix::Matrix(receptor_data, sparse = T)
  object@LR_data <- 
    list(ligand_data = ligand_data, 
         receptor_data = receptor_data)
  return(object)
}

get_complex_exp <- function(
  data, 
  complex, 
  complex_info
){
  complex_component <- unique(complex_info[[complex]])
  complex_component <- intersect(complex_component, rownames(data))
  if(length(complex_component) > 0){
    data_component <- data[complex_component, ,drop = F]
    complex_val <- colMeans(data_component)
    complex_val[apply(data_component, 2, function(x){any(x == 0)})] <- 0
  }else{
    complex_val <- rep(0, times = ncol(data))
  }
  return(complex_val)
}


sig_test <- function(
  object, 
  data = NULL,
  group_var = NULL,
  features = NULL, 
  test_use = wilcox.test,
  multiprocess = T,
  p_adj_method = 'BH',
  alternative = "greater",
  slot_used = c('identity', 'condition'),
  do_log = T, 
  type = c('gene', 'gene_set')
){
  if(multiprocess){
    if(future::nbrOfWorkers() > 1){
      my_sapply <- future.apply::future_sapply
    }else{
      warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
      my_sapply <- pbapply::pbsapply
    }
  }else{
    my_sapply <- pbapply::pbsapply
  }
  type <- match.arg(type)
  slot_used  <- match.arg(slot_used)
  data <- 
    if(is.null(data)){
      if(type == 'gene'){
        object@data
      }else{
        if(slot_used == 'identity'){
          object@GS_data[intersect(rownames(object@GS_data), names(object@regulon$identity)), ]
        }else{
          object@GS_data[intersect(rownames(object@GS_data), names(object@regulon$condition)), ]
        }
      }
    } else {
      data
    }
  if(class(data) != 'matrix'){
    data <- as.matrix(data)
  }
  ####做不做无所谓，因为用的是秩和检验
  if(do_log){
    if(type == 'gene'){
      data <- log1p(data)
    }
  }
  features <- 
    if(is.null(features)){
      rownames(data)
    }else{
      intersect(features, rownames(data))
    }
  data <- data[features, ]
  if(nrow(data) == 0 | ncol(data) == 0){
    stop('The data was empty after prossesing, please check!')
  }
  group_var <- 
    if(is.null(group_var)){
      if(slot_used == 'identity'){
        object@identity
      }else if(slot_used == 'condition'){
        object@condition
      }
    }else{
      group_var
    }
  if(length(group_var) != ncol(data)){
    stop(paste('The length of group_var must be identical with number of cells in data!', 
               'You should check the identity and condition slot or the group_var and data you provided.', collapse = '\n'))
  }
  if(class(group_var) != 'factor'){
    group_var <- factor(group_var)
  }
  group_var <- droplevels(group_var)
  #p_val_mat <- matrix(nrow = length(features), ncol = nlevels(group_var))
  p_val_mat <- 
    my_sapply(1:nlevels(group_var), FUN = 
                function(i){
                  cell.use1 <- which(group_var %in% levels(group_var)[i])
                  cell.use2 <- setdiff(1:length(group_var), cell.use1)
                  data1 <- data[, cell.use1, drop = FALSE]
                  data2 <- data[, cell.use2, drop = FALSE]
                  pvalues <- unlist(
                    x = sapply(X = 1:nrow(data1),
                               FUN = function(x) {
                                 ###进行单边检验，表达显著地高于其他分群
                                 return(test_use(data1[x, ], data2[x, ], alternative = alternative)$p.value)
                               }
                    )
                  )
                  pvalues[is.na(pvalues)] <- 1####避免均为0以及样本数小于3的影响
                  return(pvalues)
                })
  
  p_adj_mat <- apply(p_val_mat, 2 , function(x) p.adjust(x, method = p_adj_method))
  colnames(p_adj_mat) <- levels(group_var)
  rownames(p_adj_mat) <- features
  
  if(slot_used == 'identity'){
    if(type == 'gene'){
      object@p_adj_gene_idents <- list(p_adj_mat)
    }else if(type =='gene_set'){
      object@p_adj_GS_idents <- list(p_adj_mat)
    }
    return(object)
  }else if(slot_used == 'condition'){
    return(p_adj_mat)
  }
}



sig_condition_test <- function(
  object,
  alternative = "two.sided",
  type = c('gene', 'gene_set'),
  slot_used = 'condition',
  ...
){
  slot_used <- match.arg(slot_used)
  type <- match.arg(type)
  identity_var <- droplevels(object@identity)
  if(type == 'gene'){
    data <- object@data
    if(length(identity_var) != ncol(data) | length(object@condition) != ncol(data)){
      stop('The length of identity and condition slot must be identical with number of cells in data, please check!')
    }
    p_adj_mat_list <- 
      lapply(1:nlevels(identity_var), FUN = 
               function(i){
                 data_iter <- data[, identity_var %in% levels(identity_var)[i]]
                 condition_var <- object@condition[identity_var %in% levels(identity_var)[i]]
                 condition_var <- droplevels(condition_var)
                 if(nlevels(condition_var) > 1){
                   p_adj_mat_iter <- 
                     sig_test(data =  data_iter, object = object, type = type, slot_used = slot_used,
                              group_var = condition_var, alternative = alternative, ...)
                   return(p_adj_mat_iter)
                 }###############如果nlevels(condition_var) <=1 函数返还NULL对象
               })
    names(p_adj_mat_list) <- levels(identity_var)
    object@p_adj_gene_con <- p_adj_mat_list
  }else if(type == 'gene_set'){
    data <- object@GS_data[intersect(rownames(object@GS_data), names(object@regulon$condition)), ]
    if(length(identity_var) != ncol(data) | length(object@condition) != ncol(data)){
      stop('The length of identity and condition slot must be identical with number of cells in data, please check!')
    }
    p_adj_mat_greater_list <- 
      lapply(1:nlevels(identity_var), FUN = 
               function(i){
                 data_iter <- data[, identity_var %in% levels(identity_var)[i]]
                 condition_var <- object@condition[identity_var %in% levels(identity_var)[i]]
                 condition_var <- droplevels(condition_var)
                 if(nlevels(condition_var) > 1){
                   p_adj_mat_iter <- 
                     sig_test(data =  data_iter, object = object, type = type, slot_used = slot_used,
                              group_var = condition_var, alternative = 'greater', ...)
                   return(p_adj_mat_iter)
                 }
               })
    p_adj_mat_less_list <- 
      lapply(1:nlevels(identity_var), FUN = 
               function(i){
                 data_iter <- data[, identity_var %in% levels(identity_var)[i]]
                 condition_var <- object@condition[identity_var %in% levels(identity_var)[i]]
                 condition_var <- droplevels(condition_var)
                 if(nlevels(condition_var) > 1){
                   p_adj_mat_iter <- 
                     sig_test(data =  data_iter, object = object, type = type, slot_used = slot_used,
                              group_var = condition_var, alternative = 'less', ...)
                   return(p_adj_mat_iter)
                 }
               })
    names(p_adj_mat_greater_list) <-  names(p_adj_mat_less_list) <- levels(identity_var)
    object@p_adj_GS_con[['greater']] <- p_adj_mat_greater_list
    object@p_adj_GS_con[['less']] <- p_adj_mat_less_list
  }
  return(object)
}


build_sig_db <- function(
  object,
  slot_used = c('identity', 'condition'),
  p_threshold = 0.05
){
  slot_used <- match.arg(slot_used)
  p_adj_slot <- if(slot_used == 'identity') {'p_adj_gene_idents'} else{'p_adj_gene_con'}
  if(length(slot(object, name = p_adj_slot)) == 0){
    stop('Please call p-value calculation via "sig_test" or "sig_condition_test" function first !')
  } 
  if(slot_used == 'identity'){
    L_p_adj_mat <- object@p_adj_gene_idents[[1]]
    object@sig_db_idents <- 
      get_sig_db(object, 
                 L_p_adj_mat = L_p_adj_mat, 
                 R_p_adj_mat = L_p_adj_mat, 
                 p_threshold = p_threshold, 
                 slot_used = slot_used)
  }else if(slot_used == 'condition'){
    idents <- names(object@p_adj_gene_con)
    idents_num <- length(idents)
    sig_db <- 
      lapply(1:idents_num, FUN = 
               function(i){
                 L_p_adj_mat = object@p_adj_gene_con[[i]]
                 if(!is.null(L_p_adj_mat)){####如果某个cluster中只存在一种处理条件时，该对象为NULL
                   sig_db_iter <- 
                     get_sig_db(object,
                                L_p_adj_mat = L_p_adj_mat, 
                                R_p_adj_mat = L_p_adj_mat, 
                                p_threshold = p_threshold,
                                slot_used = slot_used)
                   if(nrow(sig_db_iter) > 0){
                     sig_db_iter$identity <- idents[i]
                     return(sig_db_iter)
                   }
                 }
               }
      )
    sig_db <- do.call(rbind, sig_db)
    object@sig_db_con <- sig_db
    con_net <- find_condition_net(object)
    object@condition_net <- con_net
  }
  return(object)
}

get_sig_db <- function(
  object,
  L_p_adj_mat = NULL,
  R_p_adj_mat = NULL, 
  p_threshold = 0.05,
  slot_used = c('identity', 'condition')
){
  data <- object@data
  database <- object@database
  complex_info <- object@complex_info
  slot_used <- match.arg(slot_used)
  sig_pairs <- 
    sapply(1:nrow(database), FUN = 
             function(x){
               ligand_gene <- 
                 if(database$ligand_class[x] == 'simple'){
                   database$ligand[x]
                 } else{
                   complex_info[[database$ligand[x]]]
                 }
               #ligand_gene <- ligand_gene[ligand_gene != 'None']
               receptor_gene <- 
                 if(database$receptor_class[x] == 'simple'){
                   database$receptor[x]
                 } else{
                   complex_info[[database$receptor[x]]]
                 }
               #receptor_gene <- receptor_gene[receptor_gene != 'None']
               sig_ligand_group <- unique(as.character(unlist(apply(L_p_adj_mat[ligand_gene, ,drop = F],1, function(x) names(x)[x <= p_threshold]))))
               sig_receptor_group <- unique(as.character(unlist(apply(R_p_adj_mat[receptor_gene, , drop =F],1, function(x) names(x)[x <= p_threshold]))))
               if(slot_used == 'identity'){
                 if(length(sig_ligand_group) != 0 & length(sig_receptor_group) != 0){
                   sig_pairs_iter <- paste0(get_pairs(sig_ligand_group, sig_receptor_group), collapse = '; ')
                 } else{
                   sig_pairs_iter <- 'None'
                 }
               } else if(slot_used == 'condition'){
                 if(!(length(sig_ligand_group) == 0 & length(sig_receptor_group) == 0)){
                   if(length(sig_ligand_group) != 0 & length(sig_receptor_group) != 0){
                     sig_pairs_iter <- paste0(get_pairs(sig_ligand_group, sig_receptor_group), collapse = '; ')
                     #sig_pairs_iter <- 'L~R'
                   } else if(length(sig_ligand_group) == 0 & length(sig_receptor_group) != 0){
                     #sig_pairs_iter <- 'L~None'
                     sig_ligand_group <- 'None'
                     sig_pairs_iter <- paste0(get_pairs(sig_ligand_group, sig_receptor_group), collapse = '; ')
                   } else if(length(sig_ligand_group) != 0 & length(sig_receptor_group) == 0){
                     #sig_pairs_iter <- 'None~R'
                     sig_receptor_group <- 'None'
                     sig_pairs_iter <- paste0(get_pairs(sig_ligand_group, sig_receptor_group), collapse = '; ')
                   }
                 }else{
                   sig_pairs_iter <- 'None'
                 }
               }
               return(sig_pairs_iter)
             }
    )
  database$sig_pairs <- sig_pairs
  sig_database <- database[database$sig_pairs != 'None', ]
  return(sig_database)
}


get_pairs <- function(
  x,
  y,
  sep = '~'
){
  pairs <- c()
  for(i in x){
    for(j in y){
      pairs <- c(pairs, paste(i, j, sep = sep))
    }
  }
  return(unique(pairs))
}


find_condition_net <- function(####这个函数用于构建需要计算的互作网络，以减少置换的计算量
  object
){
  database <- object@sig_db_con
  identity_var <- object@identity
  idents_num <- nlevels(identity_var)
  data <- 
    lapply(1:nrow(database), FUN = 
             function(i){
               if(stringr::str_detect(database$sig_pairs[i], 'None~')){
                 mat <- 
                   matrix(NA, ncol = 7, nrow = idents_num,
                          dimnames = list(1:idents_num,
                                          #col_names = 
                                          c('interaction_name', 
                                            'ligand',
                                            'receptor',
                                            'ligand_class', 
                                            'receptor_class',
                                            'ligand_cluster', 
                                            'receptor_cluster')))
                 mat[,'interaction_name'] <- database$interaction_name[i]
                 mat[,'ligand'] <- database$ligand[i]
                 mat[,'receptor'] <- database$receptor[i]
                 mat[,'ligand_class'] <- database$ligand_class[i]
                 mat[,'receptor_class'] <- database$receptor_class[i]
                 mat[,'ligand_cluster'] <- database$identity[i]
                 mat[,'receptor_cluster'] <- levels(identity_var)
               } else if(stringr::str_detect(database$sig_pairs[i], '~None')){
                 mat <- 
                   matrix(NA, ncol = 7, nrow = idents_num,
                          dimnames = list(1:idents_num,
                                          #col_names = 
                                          c('interaction_name', 
                                            'ligand',
                                            'receptor',
                                            'ligand_class', 
                                            'receptor_class',
                                            'ligand_cluster', 
                                            'receptor_cluster')))
                 mat[,'interaction_name'] <- database$interaction_name[i]
                 mat[,'ligand'] <- database$ligand[i]
                 mat[,'receptor'] <- database$receptor[i]
                 mat[,'ligand_class'] <- database$ligand_class[i]
                 mat[,'receptor_class'] <- database$receptor_class[i]
                 mat[,'ligand_cluster'] <- levels(identity_var)
                 mat[,'receptor_cluster'] <- database$identity[i]
               } else{
                 mat <- 
                   matrix(NA, ncol = 7, nrow = (2*idents_num-1),
                          dimnames = list(1:(2*idents_num-1),
                                          #col_names = 
                                          c('interaction_name', 
                                            'ligand',
                                            'receptor',
                                            'ligand_class', 
                                            'receptor_class',
                                            'ligand_cluster', 
                                            'receptor_cluster')))
                 mat[,'interaction_name'] <- database$interaction_name[i]
                 mat[,'ligand'] <- database$ligand[i]
                 mat[,'receptor'] <- database$receptor[i]
                 mat[,'ligand_class'] <- database$ligand_class[i]
                 mat[,'receptor_class'] <- database$receptor_class[i]
                 mat[1:idents_num, 'ligand_cluster'] <- 
                   database$identity[i]
                 mat[1:idents_num, 'receptor_cluster'] <- 
                   levels(identity_var)
                 mat[(idents_num+1):(2*idents_num-1), 'ligand_cluster'] <- 
                   setdiff(levels(identity_var), database$identity[i])
                 mat[(idents_num+1):(2*idents_num-1), 'receptor_cluster'] <- 
                   database$identity[i]
               }
               return(mat)
             }
    )
  data <- do.call(rbind, data)
  data <- as.data.frame(data)
  data$pairs_pairs <- paste(data$interaction_name, data$ligand_cluster, data$receptor_cluster, sep = '~')
  data <- data[!duplicated(data$pairs_pairs), ]
  return(data)
}


get_permutation2 <- function(
  object, 
  slot_used = c('identity', 'condition'),
  style = c('product', 'distance', 'mean'), 
  normalize = T,
  percent = T, 
  ligand_pct_threshold = 0.1,
  receptor_pct_threshold = 0.1,
  seed_used = 1234,
  iteration = 1000,######迭代的次数
  multiprocess = T,####是否使用多线程来运行get_permutation2函数
  multiprocess_to_get_intensity = F,####使用多线程来运行get_intercellular_score函数，建议是不需要
  block = T,#####为防止内存耗尽，较大数据集可以使用这个参数用于对迭代进行分割
  block_size = 30#,#####将多少次迭代作为一个迭代块
  #cross_condition = F
){
  if(multiprocess){
    if(future::nbrOfWorkers() > 1){
      my_lapply <- future.apply::future_lapply
    }else{
      warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
      my_lapply <- pbapply::pblapply
    }
  }else{
    my_lapply <- pbapply::pblapply
  }
  slot_used <- match.arg(slot_used)
  if(normalize != object@L_R_matrix_normalization){
    object <- build_LR_matrix(object= object, normalize = normalize)
  }
  object <- do.call(get_real_intensity, args = list(object = object, slot_used = slot_used, style = style, 
                                                    percent = percent, 
                                                    normalize = normalize, 
                                                    ligand_pct_threshold = ligand_pct_threshold, 
                                                    receptor_pct_threshold = receptor_pct_threshold, 
                                                    #real_call = T, 
                                                    multiprocess_to_get_intensity = multiprocess_to_get_intensity#, 
                                                    #cross_condition = cross_condition
  ))
  ligand_data <- object@LR_data$ligand_data
  receptor_data <- object@LR_data$receptor_data
  
  if(slot_used == 'identity'){
    database <- 
      if(nrow(object@sig_db_idents) > 0){
        object@sig_db_idents
      }else{
        object@database
      }
    if(length(object@identity) != ncol(ligand_data) | length(object@identity) != ncol(receptor_data)){
      stop('The length of identity slot must be identical with number of cells in ligand_data or receptor_data!')
    }
    ligand_group_var <- receptor_group_var <- object@identity
    intensity_mat <- object@inter_intensity_idents
    permutation <- get_permutation_loc(object = object, iteration = iteration, seed_used = seed_used, slot_used = slot_used)
    #object@permutation_idents <- permutation
    object@pnull_idents <- lapply(intensity_mat, function(x){
      matrix(0, nrow = nrow(x), ncol = ncol(x), 
             dimnames = list(rownames(x), colnames(x)))
    }) 
    names(object@pnull_idents) <- names(intensity_mat)
    if(block){
      block_size <- block_size
    } else{
      block_size <- iteration
    }
    for(i in 1:ceiling(iteration/block_size)){
      time_start <- Sys.time()
      iter_num <- 
        if(block_size*i > iteration){
          (block_size*(i-1)+1):iteration
        } else{
          (block_size*(i-1)+1):(block_size*i)
        }
      Pboot <- 
        my_lapply(iter_num, FUN = 
                    function(x){
                      ligand_group_var_iter <- ligand_group_var[permutation[, x]]
                      receptor_group_var_iter <- receptor_group_var[permutation[, x]]
                      intensity_mat_iter <- 
                        get_intercellular_score(object,  
                                                ligand_data = ligand_data,
                                                receptor_data = receptor_data,
                                                ligand_group_var = ligand_group_var_iter,
                                                receptor_group_var = receptor_group_var_iter,
                                                database = database,
                                                percent = percent, 
                                                ligand_pct_threshold = ligand_pct_threshold,
                                                receptor_pct_threshold = receptor_pct_threshold,
                                                real_call = F, 
                                                multiprocess = multiprocess_to_get_intensity,
                                                style = style)
                      #print(paste0('Finish iteration :', x))
                      return(intensity_mat_iter)
                    })
      Pnull <- 
        lapply(Pboot, FUN = 
                 function(iter){
                   p_mat_iter <- 
                     lapply(names(iter), FUN=
                              function(x){
                                (iter[[x]] - intensity_mat[[x]] >= 0) * 1
                              }
                     )
                   names(p_mat_iter) <- names(iter)
                   #print(paste0('Finish iteration :', x))
                   return(p_mat_iter)
                 })
      Pblock <- 
        lapply(names(intensity_mat), FUN = 
                 function(db_row){
                   list_mat_iter <- 
                     lapply(Pnull, FUN = function(x) x[[db_row]])
                   return(Reduce('+', list_mat_iter))
                 })
      object@pnull_idents <- mapply('+', object@pnull_idents, Pblock, SIMPLIFY = F)
      time_end <- Sys.time()
      time_used <- as.double(difftime(time_end,time_start, units = 'secs'))
      print(paste0('Finish block : ', i, ' ; total block : ', 
                   ceiling(iteration/block_size), ' ; ', floor(time_used/60), ' mins ', 
                   round(time_used%%60,1), ' secs used.'))
    }
    object@pnull_idents <- 
      lapply(names(object@pnull_idents), FUN = function(x){
        mat <- object@pnull_idents[[x]]/iteration
        mat[object@pct_idents[[x]] == 0] <- 1
        return(mat)
        #x <- Matrix::Matrix(x, sparse = T)
      })
    names(object@pnull_idents) <- names(intensity_mat)
    return(object)
  }else if(slot_used == 'condition'){
    database <- 
      if(nrow(object@sig_db_con) > 0){
        object@sig_db_con
      }else{
        object@database
      }
    ligand_group_var <- receptor_group_var <- object@identity
    intensity_mat <- object@inter_intensity_con
    permutation <- get_permutation_loc(object = object, iteration = iteration, seed_used = seed_used, slot_used = slot_used)
    #object@permutation_con <- permutation
    object@pnull_con[['greater']] <- object@pnull_con[['less']] <- 
      lapply(intensity_mat, function(x){
        intensity_mat_iter <- lapply(x, function(mat){
          matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                 dimnames = list(rownames(mat), colnames(mat)))
        })
        names(intensity_mat_iter) <- names(x)
        return(intensity_mat_iter)
      })
    names(object@pnull_con[['greater']]) <- names(object@pnull_con[['less']]) <- names(intensity_mat)
    
    if(block){
      block_size <- block_size
    } else{
      block_size <- iteration
    }
    for(i in 1:ceiling(iteration/block_size)){
      time_start <- Sys.time()
      iter_num <- 
        if(block_size*i > iteration){
          (block_size*(i-1)+1):iteration
        } else{
          (block_size*(i-1)+1):(block_size*i)
        }
      Pboot <- 
        my_lapply(iter_num, FUN = 
                    function(x){
                      #new_ligand_group_var <- ligand_group_var[permutation[, x]]
                      #new_receptor_group_var <- receptor_group_var[permutation[, x]]
                      condition_iter <- object@condition[permutation[, x]]
                      intensity_mat_list <- lapply(levels(object@condition), function(y){                        
                        ligand_group_var_iter <-  droplevels(ligand_group_var[condition_iter %in% y])
                        receptor_group_var_iter <- droplevels(receptor_group_var[condition_iter %in% y])
                        ligand_data_iter <- ligand_data[, condition_iter %in% y, drop = F]
                        receptor_data_iter <- receptor_data[, condition_iter %in% y, drop = F]
                        result <- 
                          get_intercellular_score(object,  
                                                  ligand_data = ligand_data_iter,
                                                  receptor_data = receptor_data_iter,
                                                  ligand_group_var = ligand_group_var_iter,
                                                  receptor_group_var = receptor_group_var_iter,
                                                  database = database, 
                                                  percent = percent, 
                                                  ligand_pct_threshold = ligand_pct_threshold,
                                                  receptor_pct_threshold = receptor_pct_threshold,
                                                  style = style,
                                                  real_call = F, 
                                                  multiprocess=multiprocess_to_get_intensity)
                        return(result)
                      })
                      
                      names(intensity_mat_list) <- levels(object@condition)
                      
                      new_intensity_mat_list <- lapply(database$interaction_name, function(db){
                        result_con <- lapply(names(intensity_mat_list), function(con){
                          return(intensity_mat_list[[con]][[db]])
                        })
                        names(result_con)<- names(intensity_mat_list)
                        return(result_con)
                      })
                      names(new_intensity_mat_list) <- database$interaction_name
                      return(new_intensity_mat_list)
                    })
      
      Pnull <- 
        lapply(Pboot, FUN = 
                 function(iter){
                   greater <- 
                     lapply(names(iter), FUN=
                              function(db){
                                result <- lapply(names(iter[[db]]), function(con){
                                  return((iter[[db]][[con]] - intensity_mat[[db]][[con]] >= 0) * 1)
                                })
                                names(result) <- names(iter[[db]])
                                return(result)
                              }
                     )
                   less <- 
                     lapply(names(iter), FUN=
                              function(db){
                                result <- lapply(names(iter[[db]]), function(con){
                                  return((iter[[db]][[con]] - intensity_mat[[db]][[con]] <= 0) * 1)
                                })
                                names(result) <- names(iter[[db]])
                                return(result)
                              }
                     )
                   names(greater) <- names(less) <- names(iter)
                   return(list(greater = greater, less = less))
                 })
      
      Pblock_greater <- 
        lapply(names(intensity_mat), function(db){
          result_con <- lapply(names(intensity_mat[[db]]),function(con){
            list_mat_iter <- 
              lapply(Pnull, FUN = function(iter){
                result_iter <- iter[['greater']][[db]][[con]]
                return(result_iter)
              })
            return(Reduce('+', list_mat_iter))
          })
          names(result_con) <- names(intensity_mat[[db]])
          return(result_con)
        })
      names(Pblock_greater) <- names(intensity_mat)
      
      object@pnull_con[['greater']] <- lapply(names(Pblock_greater), function(db){
        conditions <- names(Pblock_greater[[db]])####这一步使得条件能够在下一步运算中对齐
        mapply('+', object@pnull_con[['greater']][[db]][conditions], Pblock_greater[[db]][conditions], SIMPLIFY = F)
      })
      names(object@pnull_con[['greater']]) <- names(Pblock_greater)
      
      Pblock_less <- 
        lapply(names(intensity_mat), function(db){
          result_con <- lapply(names(intensity_mat[[db]]),function(con){
            list_mat_iter <- 
              lapply(Pnull, FUN = function(iter){
                result_iter <- iter[['less']][[db]][[con]]
                return(result_iter)
              })
            return(Reduce('+', list_mat_iter))
          })
          names(result_con) <- names(intensity_mat[[db]])
          return(result_con)
        })
      names(Pblock_less) <- names(intensity_mat)
      
      object@pnull_con[['less']] <- lapply(names(Pblock_less), function(db){
        conditions <- names(Pblock_less[[db]])####这一步使得条件能够在下一步运算中对齐
        mapply('+', object@pnull_con[['less']][[db]][conditions], Pblock_less[[db]][conditions], SIMPLIFY = F)
      })
      names(object@pnull_con[['less']]) <- names(Pblock_less)
      
      time_end <- Sys.time()
      time_used <- as.double(difftime(time_end,time_start, units = 'secs'))
      print(paste0('Finish block : ', i, ' ; total block : ', 
                   ceiling(iteration/block_size), ' ; ', floor(time_used/60), ' mins ', 
                   round(time_used%%60,1), ' secs used.'))
    }
    
    pct_con_logi <- get_condition_logi(object@pct_con)
    #pct_con_logi <- object@pct_con
    
    object@pnull_con[['greater']] <- 
      lapply(names(Pblock_greater), FUN = function(db){
        mat <- lapply(names(object@pnull_con[['greater']][[db]]), function(con){
          mat_iter <- object@pnull_con[['greater']][[db]][[con]]/iteration
          mat_iter[pct_con_logi[[db]][[con]] == 0] <- 1
          return(mat_iter)
        })
        names(mat) <- names(object@pnull_con[['greater']][[db]])
        return(mat)
        #x <- Matrix::Matrix(x, sparse = T)
      })
    names(object@pnull_con[['greater']]) <- names(Pblock_greater)
    
    object@pnull_con[['less']] <- 
      lapply(names(Pblock_less), FUN = function(db){
        mat <- lapply(names(object@pnull_con[['less']][[db]]), function(con){
          mat_iter <- object@pnull_con[['less']][[db]][[con]]/iteration
          mat_iter[pct_con_logi[[db]][[con]] == 0] <- 1
          return(mat_iter)
        })
        names(mat) <- names(object@pnull_con[['less']][[db]])
        return(mat)
        #x <- Matrix::Matrix(x, sparse = T)
      })
    names(object@pnull_con[['less']]) <- names(Pblock_less)
    
    return(object)
  }
}


get_real_intensity <- function(
  object, 
  slot_used = c('identity', 'condition'),
  style = c('product', 'distance', 'mean'), 
  normalize = T,
  percent = T, 
  ligand_pct_threshold = 0.1,
  receptor_pct_threshold = 0.1,
  multiprocess_to_get_intensity = F#,
  #real_call = T
  #cross_condition = F
){
  if(normalize != object@L_R_matrix_normalization){
    object <- build_LR_matrix(object= object, normalize = normalize)
  }
  slot_used <- match.arg(slot_used)
  ligand_data <- object@LR_data$ligand_data
  receptor_data <- object@LR_data$receptor_data
  if(slot_used == 'identity'){
    database <- 
      if(nrow(object@sig_db_idents) > 0){
        object@sig_db_idents
      }else{
        object@database
      }
    if(length(object@identity) != ncol(ligand_data) | length(object@identity) != ncol(receptor_data)){
      stop('The length of identity slot must be identical with number of cells in ligand_data or receptor_data!')
    }
    ligand_group_var <- receptor_group_var <- object@identity
    intensity_mat_list <- 
      get_intercellular_score(object,  
                              ligand_data = ligand_data,
                              receptor_data = receptor_data,
                              ligand_group_var = ligand_group_var,
                              receptor_group_var = receptor_group_var,
                              database = database, 
                              percent = percent, 
                              ligand_pct_threshold = ligand_pct_threshold,
                              receptor_pct_threshold = receptor_pct_threshold,
                              style = style,
                              real_call = T, 
                              multiprocess=multiprocess_to_get_intensity)##
    object@inter_intensity_idents <- intensity_mat_list$intensity_mat
    object@pct_idents <- intensity_mat_list$pct_mat
  } else if(slot_used == 'condition'){
    database <- 
      if(nrow(object@sig_db_con) > 0){
        object@sig_db_con
      }else{
        object@database
      }
    if(length(object@identity) != ncol(ligand_data) | 
       length(object@identity) != ncol(receptor_data) |
       length(object@identity) != length(object@condition)){
      stop('The length of identity and condition slot must be identical with number of cells in ligand_data or receptor_data!')
    }
    ligand_group_var <- receptor_group_var <- object@identity
    intensity_mat_list <- lapply(levels(object@condition), function(x){
      ligand_group_var_iter <-  droplevels(ligand_group_var[object@condition %in% x])
      receptor_group_var_iter <- droplevels(receptor_group_var[object@condition %in% x])
      ligand_data_iter <- ligand_data[, object@condition %in% x, drop = F]
      receptor_data_iter <- receptor_data[, object@condition %in% x, drop = F]
      result <- 
        get_intercellular_score(object,  
                                ligand_data = ligand_data_iter,
                                receptor_data = receptor_data_iter,
                                ligand_group_var = ligand_group_var_iter,
                                receptor_group_var = receptor_group_var_iter,
                                database = database, 
                                percent = percent, 
                                ligand_pct_threshold = ligand_pct_threshold,
                                receptor_pct_threshold = receptor_pct_threshold,
                                style = style,
                                real_call = T, 
                                multiprocess=multiprocess_to_get_intensity)
      return(result)
    })
    names(intensity_mat_list) <- levels(object@condition)
    intensity_mat <- lapply(database$interaction_name, function(db){
      result_con <- lapply(names(intensity_mat_list), function(con){
        return(intensity_mat_list[[con]][['intensity_mat']][[db]])
      })
      names(result_con)<- names(intensity_mat_list)
      return(result_con)
    })
    names(intensity_mat) <- database$interaction_name
    
    pct_mat <- lapply(database$interaction_name, function(db){
      result_con <- lapply(names(intensity_mat_list), function(con){
        return(intensity_mat_list[[con]][['pct_mat']][[db]])
      })
      names(result_con)<- names(intensity_mat_list)
      return(result_con)
    })
    names(pct_mat) <- database$interaction_name
    
    object@inter_intensity_con <- intensity_mat
    object@pct_con <- pct_mat
  }
  return(object)
}



get_intercellular_score <- function(
  object,
  ligand_data,
  receptor_data,
  ligand_group_var,
  receptor_group_var,
  database,
  style = c('product', 'distance', 'mean'),
  percent = T, 
  ligand_pct_threshold = 0.1,
  receptor_pct_threshold = 0.1,
  real_call = F, 
  multiprocess = F
){
  if(multiprocess){
    if(future::nbrOfWorkers() > 1){
      my_lapply <- future.apply::future_lapply
    }else{
      warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
      my_lapply <- lapply
    }
  }else{
    my_lapply <- lapply
  } 
  
  dataL <- ligand_data[match(database$ligand, rownames(ligand_data)), ,drop =F]
  dataR <- receptor_data[match(database$receptor, rownames(receptor_data)), ,drop =F]
  
  dataLavg <- my_group_mean(dataL, ligand_group_var)
  dataRavg <- my_group_mean(dataR, receptor_group_var)
  
  dataLavg2 <- my_group_pct(dataL, ligand_group_var)
  dataRavg2 <- my_group_pct(dataR, receptor_group_var)
  
  if(percent){
    dataLavg <- dataLavg * dataLavg2
    dataRavg <- dataRavg * dataRavg2
  } 
  style <- match.arg(style)
  intensity_mat <- 
    my_lapply(1:nrow(database), FUN = 
                function(x){
                  return(
                    #Matrix::Matrix(
                    get_cross(dataLavg[x, ,drop =T], 
                              dataRavg[x, ,drop =T], 
                              style = style)#, sparse = T)
                  )####行名对应配体，列名对应受体
                })
  names(intensity_mat) <- database$interaction_name
  if(real_call){
    pct_mat <- 
      my_lapply(1:nrow(database), FUN = 
                  function(x){
                    return(
                      #Matrix::Matrix(
                      get_cross_pct(dataLavg2[x, ,drop =T], 
                                    dataRavg2[x, ,drop =T], 
                                    ligand_pct_threshold, 
                                    receptor_pct_threshold)#, sparse = T)
                    )####行名对应配体，列名对应受体
                  })
    names(pct_mat) <- database$interaction_name
    return(list(intensity_mat=intensity_mat, pct_mat = pct_mat))
  }else{
    return(intensity_mat)
  }
}


get_condition_logi <- function(
  logi_mat_list = list()#######仅包含0或1的逻辑判断矩阵， 0为FALSE，1为TRUE
){
  result <- lapply(names(logi_mat_list), function(db){
    result_con <- mat_loop(logi_mat_list[[db]])
    names(result_con) <- names(logi_mat_list[[db]])
    return(result_con)
  })
  names(result) <- names(logi_mat_list)
  return(result)
}

mat_loop <- function(
  mat_list= list()
){
  new_list <- mat_list
  for(i in 1:length(mat_list)){
    mat1 <- mat_list[[i]]
    for(j in setdiff(1:length(mat_list), i)){
      mat2 <- mat_list[[j]]
      mat1 <- mat_logi(mat1 =  mat1, mat2 = mat2)
    }
    new_list[[i]] <- mat1
  }
  return(new_list)
}

mat_logi <- function(
  mat1, 
  mat2
){
  extra_mat1_col <- setdiff(colnames(mat1), colnames(mat2))
  if(length(extra_mat1_col) > 0){
    extra_col <- matrix(0, nrow = nrow(mat2), ncol = length(extra_mat1_col), 
                        dimnames = list(rownames(mat2), extra_mat1_col))
    mat2 <- cbind(mat2, extra_col)
  }
  extra_mat1_row <- setdiff(rownames(mat1), rownames(mat2))
  if(length(extra_mat1_row) > 0){
    extra_row <- matrix(0, ncol = ncol(mat2), nrow = length(extra_mat1_row), 
                        dimnames = list(extra_mat1_row, colnames(mat2)))
    mat2 <- rbind(mat2, extra_row)
  }
  new_mat2 <- mat2[rownames(mat1), colnames(mat1), drop=F]
  mat1[new_mat2 == 1] <- 1
  return(mat1)
}




build_GRN <- function(
  object, 
  data = NULL, 
  inter_p_threshold = 0.05,####指定胞间互作p值到什么程度才会计算相应的胞内TF-regulon
  intra_p_threshold = 0.05,
  regulon_p_threshold = 0.05, ######计算TF-targets基因相关性时的p值阈值
  regulon_cor_threshold = 0.1, ######计算TF-targets基因相关性时的相关性系数阈值
  top_gene_to_keep=1 ###########保留top多少的targets
  #min_Gs_size = 5, ###########计算AUC矩阵时，regulon应该至少包含多少基因
  #only_TF_target = F,####################是否根据TF_target数据库进行过滤，用于大数据集
  #nCores = 1,  ###########计算AUC矩阵时，使用的核心数量
  #max_cells = 10000,
  #ssgsea_arg_list = list(min_size = 5), 
  #do_sig_test = T, 
  #sig_test_arg_list = list(multiprocess = FALSE), ##########传递给sig_test的参数，用于在identity之间计算显著的regulon
  #do_sig_condition_test = T,
  #sig_condition_test_arg_list = list(multiprocess = FALSE), ##########传递给sig_condition_test的参数，用于在condition之间计算显著的regulon
  #do_build_activity_matrix = T,
  #do_get_intracellular_score = T
){
  if(is.null(data)){
    stop('Please pass the expression matrix via data parameter!')
  }
  if(class(data) != 'dgCMatrix'){
    if(class(data) == 'matrix'){
      data <- Matrix::Matrix(data, sparse = T)
    } else{
      stop('Data must be a "matrix" or "dgCMatrix" object')
    }
  }
  data <- data[rowSums(data) > 0, ]
  
  if(is.null(inter_p_threshold)){
    inter_p_threshold <- object@inter_p_threshold
  }else{
    inter_p_threshold <- inter_p_threshold
    object@inter_p_threshold <- inter_p_threshold
  }
  database <- object@database
  clu_TF_regulon_idents <- NULL
  clu_TF_regulon_con <- NULL
  if(length(object@pnull_idents) > 0){
    clu_TF_idents <- lapply(names(object@pnull_idents), function(db){
      pnull_iter <- object@pnull_idents[[db]]
      sig_receptor_clu <- colnames(pnull_iter)[apply(pnull_iter, 2, function(x){any(x <= inter_p_threshold)})]
      if(length(sig_receptor_clu) > 0){
        TF_iter <- database[database$interaction_name %in% db, 'downstream_TF']
        TF_iter <- unlist(strsplit(TF_iter, split = ','))
        if(length(TF_iter) > 0){
          clu_TF_idents_iter <- get_pairs(sig_receptor_clu, TF_iter, sep = '~')
        }else{
          clu_TF_idents_iter <- 'none'
        }
      }else{
        clu_TF_idents_iter <- 'none'
      }
      return(clu_TF_idents_iter)
    })
    clu_TF_idents <- unique(unlist(clu_TF_idents))
    clu_TF_idents <- clu_TF_idents[clu_TF_idents != 'none']
    if(length(clu_TF_idents) > 0){
      clu_TF_regulon_idents <- get_regulons_by_cluster(
        object = object, 
        data = data,
        clu_TF = clu_TF_idents, 
        regulon_p_threshold=regulon_p_threshold, 
        regulon_cor_threshold=regulon_cor_threshold, 
        top_gene_to_keep = top_gene_to_keep
      )
    }
  }else{
    warning('The pnull_idents slot is empty, please run get_permutation2 function first!')
  }
  if(length(object@pnull_con) > 0){
    clu_TF_con <- lapply(names(object@pnull_con$greater), function(db){
      sig_receptor_clu <- lapply(names(object@pnull_con$greater[[db]]), function(con){
        greater_iter <- object@pnull_con$greater[[db]][[con]]
        less_iter <- object@pnull_con$less[[db]][[con]]
        result_con <- union(
          colnames(greater_iter)[apply(greater_iter, 2, function(z){any(z <= inter_p_threshold)})], 
          colnames(less_iter)[apply(less_iter, 2, function(z){any(z <= inter_p_threshold)})]
        )
        return(result_con)
      })
      sig_receptor_clu <- unique(unlist(sig_receptor_clu))
      if(length(sig_receptor_clu) > 0){
        TF_iter <- database[database$interaction_name %in% db, 'downstream_TF']
        TF_iter <- unlist(strsplit(TF_iter, split = ','))
        if(length(TF_iter) > 0){
          clu_TF_con_iter <- get_pairs(sig_receptor_clu, TF_iter, sep = '~')
        }else{
          clu_TF_con_iter <- 'none'
        }
      }else{
        clu_TF_con_iter <- 'none'
      }
      return(clu_TF_con_iter)
    })
    clu_TF_con <- unique(unlist(clu_TF_con))
    clu_TF_con <- clu_TF_con[clu_TF_con != 'none']
    if(length(clu_TF_con) > 0){
      clu_TF_regulon_con <- get_regulons_by_cluster(
        object = object, 
        data = data, 
        clu_TF = clu_TF_con, 
        regulon_p_threshold=regulon_p_threshold, 
        regulon_cor_threshold=regulon_cor_threshold, 
        top_gene_to_keep = top_gene_to_keep
      )
    }
  }else{
    warning('The pnull_con slot is empty, please run get_permutation2 function first!')
  }
  clu_TF_regulon <- c(clu_TF_regulon_idents, clu_TF_regulon_con)
  clu_TF_regulon <- clu_TF_regulon[!duplicated(names(clu_TF_regulon))]
  if(!is.null(clu_TF_regulon)){
    object@regulon$identity <- clu_TF_regulon[names(clu_TF_regulon) %in% names(clu_TF_regulon_idents)]
    object@regulon$condition <- clu_TF_regulon[names(clu_TF_regulon) %in% names(clu_TF_regulon_con)]
  }
  # object <- get_GS_matrix(object, data = data, ssgsea_arg_list = ssgsea_arg_list)
  # #, min_Gs_size =min_Gs_size, #nCores=nCores, only_TF_target = only_TF_target, )
  # if(do_sig_test){
  #   sig_test_arg_list <- sig_test_arg_list[!names(sig_test_arg_list) %in% c('object', 'type')]
  #   sig_test_arg_list <- c(sig_test_arg_list, list(object=object, type = 'gene_set'))
  #   object<- do.call(sig_test, sig_test_arg_list)
  # } 
  # if(do_sig_condition_test){
  #   sig_condition_test_arg_list <- sig_condition_test_arg_list[!names(sig_condition_test_arg_list) %in% c('object', 'type')]
  #   sig_condition_test_arg_list <- c(sig_condition_test_arg_list, list(object=object, type = 'gene_set'))
  #   object <- do.call(sig_condition_test, sig_condition_test_arg_list)
  # }
  # if(do_build_activity_matrix){
  #   object <- build_activity_matrix(object, intra_p_threshold= intra_p_threshold)
  # }
  # if(do_get_intracellular_score){
  #   object <- get_intracellular_score(object)
  # }
  return(object)
}


get_regulons_by_cluster <- function(
  object, 
  data = NULL, ##########表达矩阵
  clu_TF, 
  regulon_p_threshold = 0.05, 
  regulon_cor_threshold = 0.1, 
  top_gene_to_keep=1
){
  enriched_clu <- unique(sapply(strsplit(clu_TF, split = '~'), function(x) x[1]))
  clu_TF_regulon <- lapply(enriched_clu, function(i){
    data_iter <- data[, object@identity %in% i]
    data_iter <- data_iter[rowSums(data_iter) > 0, ]
    TF_target_iter <- object@TF_target
    TF_target_iter <- TF_target_iter[TF_target_iter$TF %in% rownames(data_iter) & 
                                       TF_target_iter$target %in% rownames(data_iter), ]
    data_iter <- data_iter[rownames(data_iter) %in% union(TF_target_iter$TF,TF_target_iter$target), ]
    data_iter <- as.matrix(data_iter)
    data_iter <- log1p(data_iter)###spearman做不做log处理都可以
    clu_TF_iter <- clu_TF[str_detect(clu_TF, paste0(i, '~'))]
    clu_TF_iter <- sapply(strsplit(clu_TF_iter, split = '~'), function(x) x[2])
    clu_TF_iter <- clu_TF_iter[clu_TF_iter %in% TF_target_iter$TF]####移除不包含在TF_target数据库的TF
    if(length(clu_TF_iter) > 0){
      clu_TF_regulon_iter <- lapply(clu_TF_iter, function(j){
        cor_iter <- 
          corr.test(t(data_iter[j,,drop=F]),t(data_iter[TF_target_iter[TF_target_iter$TF == j, 'target'], , drop = F]),
                    adjust = "none",use = "pairwise", method = "spearman")
        r_iter <- t(cor_iter$r);p_iter <- t(cor_iter$p)
        r_iter[is.na(r_iter)] <- 0;p_iter[is.na(p_iter)] <- 1
        r_iter <- r_iter[r_iter >= regulon_cor_threshold & p_iter <= regulon_p_threshold, ,drop = F]
        if(nrow(r_iter) > 0){
          r_iter <- r_iter[order(r_iter[, 1], decreasing = T), , drop = F]
          r_iter <- r_iter[1:ceiling(nrow(r_iter) * top_gene_to_keep), ,drop = F]
        }
        return(unique(rownames(r_iter)))
      })
      names(clu_TF_regulon_iter) <- paste0(i, '~', clu_TF_iter)
    }else{
      clu_TF_regulon_iter <- list()
    }
    return(clu_TF_regulon_iter)
  })
  clu_TF_regulon <- unlist(clu_TF_regulon, recursive = F)
  clu_TF_regulon <- clu_TF_regulon[!sapply(clu_TF_regulon, is.null)]
  return(clu_TF_regulon)
}

# get_GS_matrix <- function(
#   object,
#   data = NULL, 
#   min_Gs_size = 5, 
#   only_TF_target = F,####是否只保留regulon内的基因进行排序操作，针对数据集比较大的情况
#   max_cell = 10000,
#   nCores = 1
# ){
#   #####覆盖掉AUCell自己的getRanking函数
#   setMethod("getRanking",
#             signature="aucellResults",
#             definition = function(object) {
#               if("ranking" %in% assayNames(object)) {
#                 object@assays@data$ranking
#               }else{
#                 stop("This object does not contain a ranking.")
#               }
#             }
#   )
#   if(length(object@regulon) > 0){
#     regulon <- c(object@regulon$identity, object@regulon$condition)
#     regulon <- regulon[!duplicated(names(regulon))]
#     regulon_len <- sapply(regulon, length)
#     regulon <- regulon[regulon_len >= min_Gs_size]
#     if(only_TF_target){
#       data <- data[intersect(rownames(data), union(object@TF_target$TF, object@TF_target$target)), ]
#     }
#     result <- lapply(1:ceiling(ncol(data)/max_cell), function(i){
#       iter_num <- 
#         if(max_cell*i > ncol(data)){
#           (max_cell*(i-1)+1):ncol(data)
#         } else{
#           (max_cell*(i-1)+1):(max_cell*i)
#         }
#       ranking_data <- AUCell_buildRankings(data[, iter_num, drop=F], nCores=nCores, plotStats=F)
#       GS_data <- AUCell_calcAUC(regulon, ranking_data, nCores=nCores)
#       GS_data <- GS_data@assays@data$AUC
#       return(GS_data)
#     })
#     result <- do.call(cbind, result)
#     object@GS_data <- result
#   }else{
#     warning('The regulon slot is empty so there is nothing to do.')
#   }
#   ####恢复掉AUCell自己的getRanking函数
#   setMethod("getRanking",
#             signature="aucellResults",
#             definition = function(object) {
#               if("ranking" %in% assayNames(object)) {
#                 SummarizedExperiment::assays(object)[["ranking"]]
#               }else{
#                 stop("This object does not contain a ranking.")
#               }
#             }
#   )
#   return(object)
# }

get_GS_matrix <- function(
  object,
  data = NULL, 
  #min_Gs_size = 5, 
  #only_TF_target = F,####是否只保留regulon内的基因进行排序操作，针对数据集比较大的情况
  #max_cells = 10000, 
  ...
){
  #####覆盖掉AUCell自己的getRanking函数
  if(length(object@regulon) > 0){
    regulon <- c(object@regulon$identity, object@regulon$condition)
    regulon <- regulon[!duplicated(names(regulon))]
    
    #ssgsea_arg_list <- ssgsea_arg_list[!names(ssgsea_arg_list) %in% c('exp', 'genesets')]
    #ssgsea_arg_list <- c(ssgsea_arg_list, list(exp=data, genesets = regulon))
    result <- ssgsea(genesets = regulon, data = data,...)
    object@GS_data <- result
  }else{
    warning('The regulon slot is empty so there is nothing to do.')
  }
  ####恢复掉AUCell自己的getRanking函数
  return(object)
}

ssgsea <- function(
  data, 
  genesets, 
  min_size = 5,
  max_size = Inf, 
  alpha = 0.25, 
  data_filter = F, 
  sd_cutoff = 0.05,
  normalization = F, 
  multiprocess = T)
{
  if(multiprocess){
    if(future::nbrOfWorkers() > 1){
      my_lapply <- future.apply::future_lapply
    }else{
      warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
      my_lapply <- pbapply::pblapply
    }
  }else{
    my_lapply <- pbapply::pblapply
  }
  if(class(data) != 'matrix'){
    if(class(data) != 'dgCMatrix'){
      stop('The expression matrix you privided must be a matrix or dgCMatrix object')
    }
  }
  
  if(data_filter){
    val <- apply(data, 2, sd)
    data <- data[val >= sd_cutoff, ,drop = F]
  }
  
  filtered_genesets <- lapply(genesets, function(gs){
    which(rownames(data) %in% gs)
  })
  names(filtered_genesets) <- names(genesets)
  gs_size <- sapply(filtered_genesets, length)
  filtered_genesets <- filtered_genesets[gs_size >= min_size & gs_size <= max_size]
  
  es = my_lapply(1:ncol(data), function(i) {
    if(class(data) != 'matrix'){###调用Matrix包，否则在多线程时会报错
      require(Matrix)
    }
    gene_rank <- as.integer(rank(data[, i]))  ###gene_rank中，基因值越大，秩越大
    gene_rank_loc = order(gene_rank, decreasing = TRUE)###类似于GSEA的rank_list，表达越高的基因越靠前
    
    es_iter = sapply(filtered_genesets, function(geneset) {
      
      pos_loc = gene_rank_loc %in% geneset###
      pos_neg = !pos_loc
      
      rank_score_alpha  = (gene_rank[gene_rank_loc] * pos_loc) ^ alpha
      
      step_cdf_pos = cumsum(rank_score_alpha)    / sum(rank_score_alpha)
      step_cdf_neg = cumsum(pos_neg) / sum(pos_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      return(sum(step_cdf_diff))
    })
    return(unlist(es_iter))
  })
  es <- do.call(cbind, es)
  if(normalization){
    es <- es/diff(range(es))
  }
  rownames(es) <- names(filtered_genesets)
  colnames(es) <- colnames(data)
  return(es)
}



build_activity_matrix <- function(
  object,  
  intra_p_threshold = 0.05
){
  database <- object@database
  inter_p_threshold <- object@inter_p_threshold
  object@intra_p_threshold <- intra_p_threshold
  if(length(object@pnull_idents) > 0){
    if(length(object@p_adj_GS_idents) > 0){
      active_df_idents <- lapply(names(object@pnull_idents), function(db){
        pnull_iter <- object@pnull_idents[[db]]
        active_df <- matrix(0, nrow = nrow(pnull_iter), ncol = ncol(pnull_iter), 
                            dimnames = list(rownames(pnull_iter), colnames(pnull_iter)))
        sig_receptor_clu <- colnames(pnull_iter)[apply(pnull_iter, 2, function(x){any(x <= inter_p_threshold)})]
        if(length(sig_receptor_clu) > 0){
          TF_iter <- database[database$interaction_name %in% db, 'downstream_TF']
          TF_iter <- unlist(strsplit(TF_iter, split = ','))
          if(length(TF_iter) > 0){
            clu_TF_idents_iter <- get_pairs(sig_receptor_clu, TF_iter, sep = '~')
            clu_TF_idents_iter <- clu_TF_idents_iter[clu_TF_idents_iter %in% 
                                                       rownames(object@p_adj_GS_idents[[1]])]
            if(length(clu_TF_idents_iter) > 0){
              sig_intracellular_clu <- sapply(clu_TF_idents_iter, function(i){
                clu_iter <- unlist(strsplit(i, split = '~'))[1]
                mat_iter <- object@p_adj_GS_idents[[1]][i, ,drop = F]
                if(mat_iter[, clu_iter] <= intra_p_threshold){
                  return(clu_iter)
                }
              }, simplify = F)
              sig_intracellular_clu <- unlist(sig_intracellular_clu)
            } else{
              sig_intracellular_clu <- NULL
            }
          }else{
            sig_intracellular_clu <- NULL
          }
        }else{
          sig_intracellular_clu <- NULL
        }
        intracellular_logi <- rep(colnames(active_df) %in% sig_intracellular_clu, times = nrow(pnull_iter))
        intracellular_logi_mat <- matrix(intracellular_logi, nrow = nrow(pnull_iter), byrow = T, 
                                         dimnames = list(rownames(pnull_iter), colnames(pnull_iter)))
        active_df[pnull_iter <= inter_p_threshold & intracellular_logi_mat] <- 1
        return(active_df)
      })
    }else{
      warning('The p_adj_GS_idents slot is empty, only use pnull_idents to build activity matrix!')
      active_df_idents <- lapply(names(object@pnull_idents), function(db){
        pnull_iter <- object@pnull_idents[[db]]
        active_df <- matrix(0, nrow = nrow(pnull_iter), ncol = ncol(pnull_iter), 
                            dimnames = list(rownames(pnull_iter), colnames(pnull_iter)))
        active_df[pnull_iter <= inter_p_threshold] <- 1
        return(active_df)
      })
    }
    names(active_df_idents) <- names(object@pnull_idents)
    object@activation_mat_idents <- active_df_idents       
  }else{
    warning('The pnull_idents slot is empty, please run get_permutation2 function first!')
  }
  if(length(object@pnull_con) > 0){
    if(length(object@p_adj_GS_con) > 0){
      active_df_con <- lapply(names(object@pnull_con$greater), function(db){
        result_con <- lapply(names(object@pnull_con$greater[[db]]), function(con){
          pnull_iter_greater <- object@pnull_con$greater[[db]][[con]]
          pnull_iter_less <- object@pnull_con$less[[db]][[con]]
          active_df <- matrix(0, nrow = nrow(pnull_iter_greater), ncol = ncol(pnull_iter_greater), 
                              dimnames = list(rownames(pnull_iter_greater), colnames(pnull_iter_greater)))
          sig_receptor_clu_greater <- colnames(pnull_iter_greater)[apply(pnull_iter_greater, 2, function(x){any(x <= inter_p_threshold)})]
          if(length(sig_receptor_clu_greater) > 0){
            TF_iter <- database[database$interaction_name %in% db, 'downstream_TF']
            TF_iter <- unlist(strsplit(TF_iter, split = ','))
            if(length(TF_iter) > 0){
              clu_TF_con_iter <- get_pairs(sig_receptor_clu_greater, TF_iter, sep = '~')
              clu_TF_con_iter <- clu_TF_con_iter[clu_TF_con_iter %in% 
                                                   rownames(object@p_adj_GS_con$greater[[1]])]
              if(length(clu_TF_con_iter) >0){
                sig_intracellular_clu_greater <- sapply(clu_TF_con_iter, function(i){
                  clu_iter <- unlist(strsplit(i, split = '~'))[1]
                  mat_iter <- object@p_adj_GS_con$greater[[clu_iter]][i, ,drop = F]
                  if(any(mat_iter[, con, drop = F] <= intra_p_threshold)){
                    return(clu_iter)
                  }
                }, simplify = F)
                sig_intracellular_clu_greater <- unlist(sig_intracellular_clu_greater)
              }else{
                sig_intracellular_clu_greater <- NULL
              }
            }else{
              sig_intracellular_clu_greater <- NULL
            }
          }else{
            sig_intracellular_clu_greater <- NULL
          }
          intracellular_logi_greater <- rep(colnames(active_df) %in% sig_intracellular_clu_greater, times = nrow(pnull_iter_greater))
          intracellular_logi_greater <- matrix(intracellular_logi_greater, nrow = nrow(pnull_iter_greater), byrow = T, 
                                               dimnames = list(rownames(pnull_iter_greater), colnames(pnull_iter_greater)))
          
          sig_receptor_clu_less <- colnames(pnull_iter_less)[apply(pnull_iter_less, 2, function(x){any(x <= inter_p_threshold)})]
          if(length(sig_receptor_clu_less) > 0){
            TF_iter <- database[database$interaction_name %in% db, 'downstream_TF']
            TF_iter <- unlist(strsplit(TF_iter, split = ','))
            if(length(TF_iter) > 0){
              clu_TF_con_iter <- get_pairs(sig_receptor_clu_less, TF_iter, sep = '~')
              clu_TF_con_iter <- clu_TF_con_iter[clu_TF_con_iter %in% 
                                                   rownames(object@p_adj_GS_con$less[[1]])]
              if(length(clu_TF_con_iter) >0){
                sig_intracellular_clu_less <- sapply(clu_TF_con_iter, function(i){
                  clu_iter <- unlist(strsplit(i, split = '~'))[1]
                  mat_iter <- object@p_adj_GS_con$less[[clu_iter]][i, ,drop = F]
                  if(any(mat_iter[, con, drop = F] <= intra_p_threshold)){
                    return(clu_iter)
                  }
                }, simplify = F)
                sig_intracellular_clu_less <- unlist(sig_intracellular_clu_less)
              }else{
                sig_intracellular_clu_less <- NULL
              }
            }else{
              sig_intracellular_clu_less <- NULL
            }
          }else{
            sig_intracellular_clu_less <- NULL
          }
          
          intracellular_logi_less <- rep(colnames(active_df) %in% sig_intracellular_clu_less, times = nrow(pnull_iter_less))
          intracellular_logi_less <- matrix(intracellular_logi_less, nrow = nrow(pnull_iter_less), byrow = T, 
                                            dimnames = list(rownames(pnull_iter_less), colnames(pnull_iter_less)))
          
          active_df[(pnull_iter_greater <= inter_p_threshold & intracellular_logi_greater) |
                      (pnull_iter_less <= inter_p_threshold & intracellular_logi_less)] <- 1
          return(active_df)
        })
        names(result_con) <- names(object@pnull_con$greater[[db]])#names(object@pnull_con$greater)
        return(result_con)
      })
    }else{
      warning('The p_adj_GS_con slot is empty, please run sig_condition_test function first!')
      active_df_con <- lapply(names(object@pnull_con$greater), function(db){
        result_con <- lapply(names(object@pnull_con$greater[[db]]), function(con){
          pnull_iter_greater <- object@pnull_con$greater[[db]][[con]]
          pnull_iter_less <- object@pnull_con$less[[db]][[con]]
          active_df <- matrix(0, nrow = nrow(pnull_iter_greater), ncol = ncol(pnull_iter_greater), 
                              dimnames = list(rownames(pnull_iter_greater), colnames(pnull_iter_greater)))
          
          active_df[(pnull_iter_greater <= inter_p_threshold) |
                      (pnull_iter_less <= inter_p_threshold)] <- 1
          return(active_df)
        })
        names(result_con) <- names(object@pnull_con$greater[[db]])#names(object@pnull_con$greater)
        return(result_con)
      })
    }
    names(active_df_con) <- names(object@pnull_con$greater)#max_names
    object@activation_mat_con <- active_df_con
  }else{
    warning('The pnull_con slot is empty, please run get_permutation2 function first!')
  }
  return(object)
}



get_intracellular_score <- function(
  object
){
  database <- object@database
  if(length(object@activation_mat_idents) > 0){
    intra_intensity_idents <- lapply(names(object@activation_mat_idents), function(x){
      intra_active_idents_iter <- object@activation_mat_idents[[x]]
      intra_intensity_idents_iter <- matrix(0, nrow = nrow(intra_active_idents_iter), ncol = ncol(intra_active_idents_iter), 
                                            dimnames = list(rownames(intra_active_idents_iter), colnames(intra_active_idents_iter)))
      sig_receptor_clu <- colnames(intra_active_idents_iter)[apply(intra_active_idents_iter, 2, function(x){any(x == 1)})]
      if(length(sig_receptor_clu) > 0){
        TF_iter <- database[database$interaction_name %in% x, 'downstream_TF']
        TF_iter <- unlist(strsplit(TF_iter, split = ','))
        if(length(TF_iter) > 0){
          clu_TF_idents_iter <- get_pairs(sig_receptor_clu, TF_iter, sep = '~')
          clu_TF_idents_iter <- clu_TF_idents_iter[clu_TF_idents_iter %in% 
                                                     rownames(object@p_adj_GS_idents[[1]])]
          if(length(clu_TF_idents_iter) > 0){
            for(i in unique(sapply(strsplit(clu_TF_idents_iter, split = '~'), function(x){x[1]}))){
              tmp_clu_TF_idents_iter <- clu_TF_idents_iter[str_detect(clu_TF_idents_iter, i)]
              val <- mean(object@GS_data[tmp_clu_TF_idents_iter, object@identity %in% i, drop = F])
              intra_intensity_idents_iter[, i] <- val
            }
          }else{
            intra_intensity_idents_iter <- intra_intensity_idents_iter
          }
        }else{
          intra_intensity_idents_iter <- intra_intensity_idents_iter
        }
      }else{
        intra_intensity_idents_iter <- intra_intensity_idents_iter
      }
      intra_intensity_idents_iter[intra_active_idents_iter == 0] <- 0
      return(intra_intensity_idents_iter)
    })
    names(intra_intensity_idents) <- names(object@activation_mat_idents)
    object@intra_intensity_idents <- intra_intensity_idents
  }else{
    warning('The activation_mat_idents slot is empty, please run build_activity_matrix function first!')
  }
  if(length(object@activation_mat_con) > 0){
    intra_active_con <- get_condition_logi(object@activation_mat_con)
    intra_intensity_con <- lapply(names(object@activation_mat_con), function(x){
      result_y <- lapply(names(object@activation_mat_con[[x]]), function(y){
        intra_active_con_iter <- intra_active_con[[x]][[y]]#object@activation_mat_con[[x]]
        intra_intensity_con_iter <- matrix(0, nrow = nrow(intra_active_con_iter), ncol = ncol(intra_active_con_iter), 
                                           dimnames = list(rownames(intra_active_con_iter), colnames(intra_active_con_iter)))
        sig_receptor_clu <- colnames(intra_active_con_iter)[apply(intra_active_con_iter, 2, function(x){any(x == 1)})]
        if(length(sig_receptor_clu) > 0){
          TF_iter <- database[database$interaction_name %in% x, 'downstream_TF']
          TF_iter <- unlist(strsplit(TF_iter, split = ','))
          if(length(TF_iter) > 0){
            clu_TF_con_iter <- get_pairs(sig_receptor_clu, TF_iter, sep = '~')
            clu_TF_con_iter <- clu_TF_con_iter[clu_TF_con_iter %in% 
                                                 rownames(object@p_adj_GS_con$greater[[1]])]
            if(length(clu_TF_con_iter) > 0){
              for(i in unique(sapply(strsplit(clu_TF_con_iter, split = '~'), function(x){x[1]}))){
                tmp_clu_TF_con_iter <- clu_TF_con_iter[str_detect(clu_TF_con_iter, i)]
                data_iter <- object@GS_data[tmp_clu_TF_con_iter, object@identity %in% i & object@condition %in% y, drop = F]
                #condition_iter <- droplevels(object@condition[object@identity %in% i & object@condition %in% y])
                val <- mean(data_iter)
                intra_intensity_con_iter[, i] <- val
              }
            }else{
              intra_intensity_con_iter <- intra_intensity_con_iter
            }
          }else{
            intra_intensity_con_iter <- intra_intensity_con_iter
          }
        }else{
          intra_intensity_con_iter <- intra_intensity_con_iter
        }
        intra_intensity_con_iter[intra_active_con_iter == 0] <- 0
        return(intra_intensity_con_iter)
      })
      names(result_y) <- names(object@activation_mat_con[[x]])
      return(result_y)
    })
    names(intra_intensity_con) <- names(object@activation_mat_con)
    object@intra_intensity_con <- intra_intensity_con
  }else{
    warning('The activation_mat_con slot is empty, please run build_activity_matrix function first!')
  }
  return(object)
}



get_permutation_loc <- function(
  object, 
  iteration = 1000, 
  seed_used =1234, 
  slot_used = c('identity', 'condition')
){
  slot_used <- match.arg(slot_used)
  if(slot_used == 'identity'){
    set.seed(seed_used)
    permutation <- replicate(iteration, sample.int(length(object@identity), size = length(object@identity)))
  }else{
    set.seed(seed_used)
    permutation <- replicate(iteration, get_permutation_loc_con(object@identity))
    permutation[get_orig_loc(object@identity), ] <- permutation
  }
  return(permutation)
}

get_permutation_loc_con <- function(
  group_var
){
  result <- 
    unlist(
      sapply(levels(group_var), FUN = 
               function(x){
                 loc <- which(group_var == x)
                 loc_permutation <- 
                   sample(loc, length(loc))
                 #names(loc_permutation) <- x
                 return(loc_permutation)
               }, simplify = F), use.names = F
    )
  return(result)
}

get_orig_loc <- function(
  group_var
){
  result <- 
    sapply(levels(group_var), FUN = 
             function(x){
               loc <- which(group_var == x)
               return(loc)
             }, simplify = F)
  return(unlist(result,use.names=T))
}



get_cross <- function(
  a,###row
  b, ###columns
  style = 'distance'
){
  #mat <- matrix(0, nrow = length(a), ncol = length(b))
  my_fun <- switch(style,
                   'product' = function(i,j){
                     return(i*j)
                   },
                   'distance'=function(i,j){
                     return(sqrt(i^2 + j^2))
                   },
                   'mean' = function(i,j){
                     return((i+j)/2)
                   }
  )
  mat <- 
    sapply(1:length(b), FUN = 
             function(x) {
               if(b[x] == 0){
                 val <- rep(0, length(a))
               } else {
                 val <- my_fun(b[x], a)
               }
               val[a==0] <- 0
               return(val)
             }
    )
  rownames(mat) <- names(a)
  colnames(mat) <- names(b)
  return(mat)
}


get_cross_pct <- function(
  a,###row
  b, ###columns
  a_threshold = 0.1, 
  b_threshold = 0.1
){
  mat <- 
    sapply(1:length(b), FUN = 
             function(x) {
               logi_1 <- a >= a_threshold
               logi_2 <- b[x] >= b_threshold
               return(as.numeric(logi_1 & logi_2))
             }
    )
  rownames(mat) <- names(a)
  colnames(mat) <- names(b)
  return(mat)
}


triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}


my_group_mean <- function(
  data, 
  group_var
){
  result <- 
    sapply(levels(group_var), FUN = 
             function(x){
               data_iter <- data[, group_var %in% x, drop = F]
               val <- rowMeans(data_iter)
             })
  return(result)  
}


my_group_pct <- function(
  data, 
  group_var
){
  data <- (data > 0)*1
  result <- 
    sapply(levels(group_var), FUN = 
             function(x){
               data_iter <- data[, group_var %in% x, drop = F]
               val <- rowSums(data_iter)/ncol(data_iter)
             })
  return(result)  
}


# ssgesa <- function(
#   exp, 
#   genesets, 
#   min_size = 5, 
#   max_size = Inf, 
#   alpha = 0.25, 
#   exp_filter = F, 
#   sd_cutoff = 0.05,
#   normalization = F, 
#   multiprocess = T)
# {
#   if(multiprocess){
#     if(future::nbrOfWorkers() > 1){
#       my_lapply <- future.apply::future_lapply
#     }else{
#       warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
#       my_lapply <- pbapply::pblapply
#     }
#   }else{
#     my_lapply <- pbapply::pblapply
#   }
#   if(class(exp) != 'matrix'){
#     if(class(exp) != 'dgCMatrix'){
#       stop('The expression matrix you privided must be a matrix or dgCMatrix object')
#     }
#   }
#   
#   if(exp_filter){
#     val <- apply(exp, 2, sd)
#     exp <- exp[val >= sd_cutoff, ,drop = F]
#   }
#   
#   filtered_genesets <- lapply(genesets, function(gs){
#     which(rownames(exp) %in% gs)
#   })
#   names(filtered_genesets) <- names(genesets)
#   gs_size <- sapply(filtered_genesets, length)
#   filtered_genesets <- filtered_genesets[gs_size >= min_size & gs_size <= max_size]
#   
#   es = my_lapply(1:ncol(exp), function(i) {
#     gene_rank <- as.integer(rank(exp[, i]))  ###gene_rank中，基因值越大，秩越大
#     gene_rank_loc = order(gene_rank, decreasing = TRUE)###类似于GSEA的rank_list，表达越高的基因越靠前
#     
#     es_iter = sapply(filtered_genesets, function(geneset) {
#       
#       pos_loc = gene_rank_loc %in% geneset###
#       pos_neg = !pos_loc
#       
#       rank_score_alpha  = (gene_rank[gene_rank_loc] * pos_loc) ^ alpha
#       
#       step_cdf_pos = cumsum(rank_score_alpha)    / sum(rank_score_alpha)
#       step_cdf_neg = cumsum(pos_neg) / sum(pos_neg)
#       
#       step_cdf_diff = step_cdf_pos - step_cdf_neg
#       
#       return(sum(step_cdf_diff))
#     })
#     return(unlist(es_iter))
#   })
#   es <- do.call(cbind, es)
#   if(normalization){
#     es <- es/diff(range(es))
#   }
#   rownames(es) <- names(filtered_genesets)
#   colnames(es) <- colnames(exp)
#   return(es)
# }


get_sig_LR_result <-function(
  object,
  signal=c('intercellular', 'both'),
  slot_used = c('identity', 'condition'),
  inter_p_threshold = 0.05,
  return_type = c('data.frame', 'list')
){
  signal <- match.arg(signal)
  slot_used <- match.arg(slot_used)
  return_type <- match.arg(return_type)
  sig_LR_result <- NULL
  if(signal == 'intercellular'){
    if(slot_used == 'identity'){
      if(length(object@pnull_idents) >0){
        sig_LR_result <- lapply(names(object@pnull_idents), function(db){
          pnull_mat <- object@pnull_idents[[db]]
          cellnames_mat <- outer(row.names(pnull_mat), colnames(pnull_mat), function(X, Y) paste0(X, "~", Y))
          result <- cellnames_mat[pnull_mat < inter_p_threshold]
          if(length(result) > 0){
            return(result)
          }
        })
        names(sig_LR_result) <- names(object@pnull_idents)
      }else{
        warning('The pnull_idents slot is empty, please run get_permutation2 function first!')
      }
    }else if(slot_used == 'condition'){
      if(length(object@pnull_con) >0){
        sig_LR_result <- lapply(names(object@pnull_con$greater), function(db){
          pnull_mat_list <- c(object@pnull_con$greater[[db]], object@pnull_con$less[[db]])
          cellnames <- lapply(pnull_mat_list, function(x){
            cellnames_mat <- outer(row.names(x), colnames(x), function(X, Y) paste0(X, "~", Y))
            result <-  cellnames_mat[x < inter_p_threshold]
            if(length(result) > 0){
              return(result)
            }
          })
          cellnames <- unique(unlist(cellnames))
          if(length(cellnames) > 0){
            return(cellnames)
          }
        })
        names(sig_LR_result) <- names(object@pnull_con$greater)
      }else{
        warning('The pnull_con slot is empty, please run get_permutation2 function first!')
      }
    }
  }else if(signal == 'both'){
    if(slot_used == 'identity'){
      if(length(object@activation_mat_idents) >0){
        sig_LR_result <- lapply(names(object@activation_mat_idents), function(db){
          activation_mat <- object@activation_mat_idents[[db]]
          cellnames_mat <- outer(row.names(activation_mat), colnames(activation_mat), function(X, Y) paste0(X, "~", Y))
          result <- cellnames_mat[activation_mat == 1]
          if(length(result) > 0){
            return(result)
          }
        })
        names(sig_LR_result) <- names(object@activation_mat_idents)
      }else{
        warning('The activation_mat_idents slot is empty, please run build_activity_matrix function first!')
      }
    }else if(slot_used == 'condition'){
      if(length(object@activation_mat_con) >0){
        sig_LR_result <- lapply(names(object@activation_mat_con), function(db){
          activation_mat_list <- object@activation_mat_con[[db]]
          cellnames <- lapply(activation_mat_list, function(x){
            cellnames_mat <- outer(row.names(x), colnames(x), function(X, Y) paste0(X, "~", Y))
            result <-  cellnames_mat[x ==1]
            if(length(result) > 0){
              return(result)
            }
          })
          cellnames <- unique(unlist(cellnames))
          if(length(cellnames) > 0){
            return(cellnames)
          }
        })
        names(sig_LR_result) <- names(object@activation_mat_con)
      }else{
        warning('The activation_mat_con slot is empty, please run build_activity_matrix function first!')
      }
    }
  }
  sig_LR_result <- sig_LR_result[sapply(sig_LR_result, length) > 0]
  if(return_type == 'list'){
    sig_LR_result <- sig_LR_result
  }else if(return_type == 'data.frame'){
    sig_LR_result <- lapply(names(sig_LR_result), function(db){
      if(length(sig_LR_result[[db]]) > 0){
        split_str <- strsplit(sig_LR_result[[db]], split = '~')
        df <- data.frame(ligand_cluster = sapply(split_str, function(x) x[1]),
                         receptor_cluster = sapply(split_str, function(x) x[2]),
                         ligand = unlist(strsplit(db, split = '~'))[1],
                         receptor = unlist(strsplit(db, split = '~'))[2],
                         interaction_name = db
        )
        return(df)
      }
    })
    sig_LR_result <- do.call(rbind, sig_LR_result)
    sig_LR_result$clu_pairs <- paste0(sig_LR_result$ligand_cluster, '~', sig_LR_result$receptor_cluster)
  }
  return(sig_LR_result)
}


get_sig_clu_result <- function(
  object,
  signal=c('intercellular', 'both'),
  slot_used = c('identity', 'condition'),
  inter_p_threshold = 0.05
){
  sig_LR_result <- get_sig_LR_result(
    object=object, 
    slot_used = slot_used, 
    signal = signal, 
    return_type = 'data.frame',
    inter_p_threshold = inter_p_threshold
  )
  all_sig_clu <- unique(sig_LR_result$clu_pairs)
  sig_clu_result <- lapply(all_sig_clu, function(x){
    return(unique(sig_LR_result[sig_LR_result$clu_pairs %in% x, 'interaction_name']))
  })
  names(sig_clu_result) <- all_sig_clu
  return(sig_clu_result)
}

get_plot_df <- function(
  object, 
  slot_used = c('identity', 'condition'), 
  data_to_plot = c('pnull', 'intercellular', 'intracellular', 'both','activation'), 
  pnull_con_type = c('greater', 'less'),
  #pnull_do_log = T, 
  clusters_pairs =NULL, 
  LR_pairs = NULL, 
  do_log = T,
  pnull_min_cutoff = 0.001,
  do_scale = T,
  value_max_cutoff = 1,###输入想要cutoff的分位数
  value_min_cutoff = 0
){
  if(is.null(clusters_pairs) | is.null(LR_pairs)){
    stop('You should provide the cluster_pairs and ligand_receptor_pairs to be plot!')
  }
  slot_used <- match.arg(slot_used)
  data_to_plot <- match.arg(data_to_plot)
  
  if(slot_used == 'identity'){
    if(data_to_plot == 'pnull'){
      if(length(object@pnull_idents) > 0){
        data <- object@pnull_idents
      }else{
        stop('The pnull_idents slot is empty, please run get_permutation2 function first!')
      }
    }else if(data_to_plot == 'intercellular'){
      if(length(object@inter_intensity_idents) > 0){
        data <- object@inter_intensity_idents
      }else{
        stop('The inter_intensity_idents slot is empty, please run get_permutation2 function first!')
      }
    }else if(data_to_plot == 'intracellular'){
      if(length(object@intra_intensity_idents) > 0){
        data <- object@intra_intensity_idents
      }else{
        stop('The intra_intensity_idents slot is empty, please run get_intracellular_score function first!')
      }
    }else if(data_to_plot == 'activation'){
      if(length(object@activation_mat_idents) > 0){
        data <- object@activation_mat_idents
      }else{
        stop('The activation_mat_idents slot is empty, please run build_activity_matrix function first!')
      }
    }else if(data_to_plot == 'both'){
      if(length(object@inter_intensity_idents) == 0){
        stop('The inter_intensity_idents slot is empty, please run get_permutation2 function first!')
      }
      if(length(object@intra_intensity_idents) == 0){
        stop('The intra_intensity_idents slot is empty, please run get_intracellular_score function first!')
      }
      data <- lapply(names(object@inter_intensity_idents), function(x){
        return(object@inter_intensity_idents[[x]] * object@intra_intensity_idents[[x]])
      })
      names(data) <- names(object@inter_intensity_idents)
    }
    data <- data[intersect(LR_pairs, names(data))]
    plot_df <- lapply(clusters_pairs, function(x){
      splited_clusters_pairs_iter <- unlist(strsplit(x, split = '~'))
      ligand_cluster <- splited_clusters_pairs_iter[1]
      receptor_cluster <- splited_clusters_pairs_iter[2]
      plot_df_iter <- lapply(names(data), function(y){
        df <- data.frame(val = data[[y]][ligand_cluster, receptor_cluster])
        colnames(df) <- 'value'
        df$interaction_name <- y
        df$ligand_cluster <- ligand_cluster
        df$receptor_cluster <- receptor_cluster
        df$clusters_pairs <- x
        return(df)
      })
      plot_df_iter <- do.call(rbind, plot_df_iter)
    })
    plot_df <- do.call(rbind, plot_df)
  }else if(slot_used == 'condition'){
    
    pnull_con_type <- match.arg(pnull_con_type)
    
    if(data_to_plot == 'pnull'){
      if(length(object@pnull_con) >= 0){
        data <- object@pnull_con[[pnull_con_type]]
      }else{
        stop('The pnull_con slot is empty, please run get_permutation2 function first!')
      }
    }else if(data_to_plot == 'intercellular'){
      if(length(object@inter_intensity_con) > 0){
        data <- object@inter_intensity_con
      }else{
        stop('The inter_intensity_con slot is empty, please run get_permutation2 function first!')
      }
    }else if(data_to_plot == 'intracellular'){
      if(length(object@intra_intensity_con) > 0){
        data <- object@intra_intensity_con
      }else{
        stop('The intra_intensity_con slot is empty, please run get_intracellular_score function first!')
      }
    }else if(data_to_plot == 'activation'){
      if(length(object@activation_mat_con) > 0){
        data <- object@activation_mat_con
      }else{
        stop('The activation_mat_con slot is empty, please run build_activity_matrix function first!')
      }
    }else if(data_to_plot == 'both'){
      if(length(object@inter_intensity_con) == 0){
        stop('The inter_intensity_con slot is empty, please run get_permutation2 function first!')
      }
      if(length(object@intra_intensity_con) == 0){
        stop('The intra_intensity_con slot is empty, please run get_intracellular_score function first!')
      }
      data <- lapply(names(object@inter_intensity_con), function(x){
        result <- lapply(names(object@inter_intensity_con[[x]]), function(y){
          return(object@inter_intensity_con[[x]][[y]] * object@intra_intensity_con[[x]][[y]])
        })
        names(result) <- names(object@inter_intensity_con[[x]])
        return(result)
      })
      names(data) <- names(object@inter_intensity_con)
    }
    data <- data[intersect(LR_pairs, names(data))]
    plot_df <- lapply(clusters_pairs, function(x){
      splited_clusters_pairs_iter <- unlist(strsplit(x, split = '~'))
      ligand_cluster <- splited_clusters_pairs_iter[1]
      receptor_cluster <- splited_clusters_pairs_iter[2]
      plot_df_iter <- lapply(names(data), function(y){
        plot_df_iter_con <- lapply(names(data[[y]]), function(z){
          df <- data.frame(val = data[[y]][[z]][ligand_cluster, receptor_cluster])
          colnames(df) <- 'value'
          df$interaction_name <- y
          df$ligand_cluster <- ligand_cluster
          df$receptor_cluster <- receptor_cluster
          df$clusters_pairs <- x
          df$condition = z
          return(df)
        })
        plot_df_iter_con <- do.call(rbind, plot_df_iter_con)
      })
      plot_df_iter <- do.call(rbind, plot_df_iter)
    })
    plot_df <- do.call(rbind, plot_df)
  }
  if(data_to_plot == 'pnull'){
    if(do_log){
      plot_df$value[plot_df$value < pnull_min_cutoff] <- pnull_min_cutoff
      plot_df$value <- -log10(plot_df$value)
    }
  }
  
  plot_df <- lapply(unique(plot_df$interaction_name), function(x){
    df_iter <- plot_df[plot_df$interaction_name == x, ,drop =F]
    value <- df_iter$value
    if(!data_to_plot %in% c('pnull', 'activation') & do_scale){
      value <- (value -mean(value))/sd(value)
    }
    if(!data_to_plot %in% c('pnull', 'activation')){
      tmp_value_max_cutoff <- quantile(value, value_max_cutoff)
      tmp_value_min_cutoff <- quantile(value, value_min_cutoff)
      
      value[value > tmp_value_max_cutoff] <- tmp_value_max_cutoff
      value[value < tmp_value_min_cutoff] <- tmp_value_min_cutoff
    }
    
    df_iter$value <- value
    return(df_iter)
  })
  plot_df <- do.call(rbind,plot_df)
  
  value_name <- data_to_plot
  if(value_name == 'pnull' & slot_used == 'condition'){
    value_name <- paste0(pnull_con_type, '_', value_name)
  }
  colnames(plot_df)[colnames(plot_df) %in% 'value'] <- value_name
  #plot_df$value_type <- data_to_plot
  return(plot_df)
}

get_crosstalk_result <- function(
  object,
  p_threshold = 0.05,
  multiprocess = F,
  slot_used = c('identity', 'condition')  
){
  slot_used <- match.arg(slot_used)  
  if(multiprocess & future::nbrOfWorkers() > 1){
    my_lapply <- future.apply::future_lapply
  } else{
    my_lapply <- pbapply::pblapply
  }  
  if((slot_used == 'identity') & (length(object@pnull_idents) > 0)){
    identity <- object@identity
    database <- object@database 
    result_idents <- 
      my_lapply(levels(identity), 
                function(idents){
                  result_idents_iter <-
                    lapply(1:length(object@pnull_idents), 
                           function(x){
                             L_R_pairs <- unlist(strsplit(names(object@pnull_idents)[x], split = '~'), 
                                                 use.names= F)
                             pnull_df <- object@pnull_idents[[x]]
                             intensity_df <- object@inter_intensity_idents[[x]]
                             intensity_df_2 <- intensity_df/(max(intensity_df)-min(intensity_df))
                             p_ligand <- pnull_df[, idents]
                             p_receptor <- pnull_df[idents, ]
                             sig_ligand_clu <- names(p_ligand[p_ligand <= p_threshold])
                             sig_receptor_clu <- names(p_receptor[p_receptor <= p_threshold])
                             pnull_val <- 
                               c(pnull_df[sig_ligand_clu, idents],
                                 pnull_df[idents, sig_receptor_clu])
                             intensity_val <- 
                               c(intensity_df[sig_ligand_clu, idents],
                                 intensity_df[idents, sig_receptor_clu])
                             intensity_val_2 <- 
                               c(intensity_df_2[sig_ligand_clu, idents],
                                 intensity_df_2[idents, sig_receptor_clu])                                 
                             ligand_clu <- 
                               c(sig_ligand_clu,
                                 rep(idents,length(sig_receptor_clu)))
                             receptor_clu <- 
                               c(rep(idents,length(sig_ligand_clu)),
                                 sig_receptor_clu)
                             Interaction_name <- names(object@pnull_idents)[x]  
                             Ligand_class <- 
                               unique(database[database$interaction_name %in%
                                                 Interaction_name, 'ligand_class'])
                             Receptor_class <- 
                               unique(database[database$interaction_name %in%
                                                 Interaction_name, 'receptor_class'])                               
                             if(length(pnull_val) > 0){
                               data <- data.frame(
                                 'Interaction_name' = names(object@pnull_idents)[x],
                                 'Ligand' = L_R_pairs[1],
                                 'Receptor' = L_R_pairs[2],
                                 'Ligand_class' = Ligand_class,
                                 'Receptor_class' = Receptor_class,  
                                 'Ligand_cluster' = ligand_clu,
                                 'Receptor_cluster' = receptor_clu,
                                 'P_val' = pnull_val,
                                 'Intensity' = intensity_val,
                                 'Intensity_nor' = intensity_val_2,
                                 check.names = F  
                               )
                             } 
                           })
                  result_idents_iter <- do.call(rbind, result_idents_iter)
                  if(class(result_idents_iter)=='data.frame' && nrow(result_idents_iter) > 0){
                    result_idents_iter <- dplyr::distinct(result_idents_iter)
                  }
                  return(result_idents_iter)
                }) 
    names(result_idents) <- levels(identity)
    return(result_idents)
  }else if((slot_used == 'condition') & (length(object@pnull_con) > 0)){
    identity <- object@identity
    condition <- object@condition
    group_var <- as.factor(paste(object@identity, object@condition, sep = '|'))
    database <- object@database 
    col_names <-  colnames(object@inter_intensity_con[[1]])
    row_names <-  rownames(object@inter_intensity_con[[1]])
    result_con <- 
      my_lapply(levels(identity), 
                function(idents){
                  result_con_iter <-
                    lapply(1:length(object@pnull_con), 
                           function(x){
                             data_iter <- NULL
                             sig_ligand_clu <- NULL
                             sig_receptor_clu <- NULL
                             idents_col_mode <- get_pairs(idents, levels(condition), sep = '|')
                             idents_col_mode <- idents_col_mode[idents_col_mode %in% col_names]
                             idents_row_mode <- get_pairs(idents, levels(condition), sep = '|')
                             idents_row_mode <- idents_row_mode[idents_row_mode %in% row_names]
                             L_R_pairs <- unlist(strsplit(names(object@pnull_con)[x], split = '~'), 
                                                 use.names= F)
                             pnull_df <- object@pnull_con[[x]]
                             intensity_df <- object@inter_intensity_con[[x]]
                             intensity_df_2 <- intensity_df/(max(intensity_df)-min(intensity_df))
                             p_ligand <- pnull_df[, idents_col_mode, drop = F]
                             p_receptor <- pnull_df[idents_col_mode, , drop =F]
                             sig_ligand_var <- 
                               rownames(p_ligand[
                                 apply(p_ligand, 1, 
                                       function(x){any(x <= p_threshold)})
                                 , ])                             
                             sig_receptor_var <- 
                               colnames(p_receptor[
                                 ,apply(p_receptor, 2, 
                                        function(x){any(x <= p_threshold)})
                               ])
                             if(length(sig_ligand_var) > 0){
                               sig_ligand_clu <- 
                                 sapply(strsplit(sig_ligand_var, split = '\\|'), 
                                        function(x) x[x %in% identity])  
                               ligand_mat <- 
                                 sapply(sig_ligand_clu, FUN = 
                                          function(i){
                                            idents_con <- get_pairs(i,levels(condition), sep = '|')
                                            idents_con <- idents_con[idents_con %in% row_names]
                                            val_iter <- 
                                              sapply(levels(condition), FUN=
                                                       function(x){
                                                         a <- idents_con[grep(x, idents_con)]
                                                         b <- idents_col_mode[grep(x, idents_col_mode)]     
                                                         if(length(a) !=0 & length(b)!= 0){
                                                           val <- pnull_df[a,b]
                                                           val2 <- intensity_df[a,b]
                                                           val3 <- intensity_df_2[a,b]
                                                         } else{
                                                           val <- NA
                                                           val2 <- NA
                                                           val3 <- NA
                                                         }
                                                         return(c(val, val2,val3))
                                                       }
                                              )
                                            val_iter_loc <- 1:length(val_iter)
                                            val_iter <- c(val_iter[val_iter_loc %% 3 == 1], 
                                                          val_iter[val_iter_loc %% 3 == 2],
                                                          val_iter[val_iter_loc %% 3 == 0])
                                            return(val_iter)
                                          }
                                 )
                               ligand_mat <- t(ligand_mat)
                               rownames(ligand_mat) <- sig_ligand_clu
                               colnames(ligand_mat) <- 
                                 c(get_pairs('P_val:', levels(condition), sep = ''),
                                   get_pairs('Intensity:', levels(condition), sep = ''), 
                                   get_pairs('Intensity_nor:', levels(condition), sep = ''))
                               ligand_mat <- data.frame(ligand_mat, check.names = F)
                               data_iter <- ligand_mat
                               ligand_clu <- c(sig_ligand_clu)
                               receptor_clu <- 
                                 c(rep(idents,length(sig_ligand_clu)))
                             }
                             
                             if(length(sig_receptor_var) > 0){
                               sig_receptor_clu <- 
                                 sapply(strsplit(sig_receptor_var, split = '\\|'), 
                                        function(x) x[x %in% identity])  
                               receptor_mat <- 
                                 sapply(sig_receptor_clu, FUN = 
                                          function(i){
                                            idents_con <- get_pairs(i,levels(condition), sep = '|')
                                            idents_con <- idents_con[idents_con %in% col_names]
                                            val_iter <- 
                                              sapply(levels(condition), FUN=
                                                       function(x){
                                                         a <- idents_con[grep(x, idents_con)]
                                                         b <- idents_row_mode[grep(x, idents_row_mode)]     
                                                         if(length(a) !=0 & length(b)!= 0){
                                                           val <- pnull_df[b,a]
                                                           val2 <- intensity_df[b,a]
                                                           val3 <- intensity_df_2[b,a]
                                                         } else{
                                                           val <- NA
                                                           val2 <- NA
                                                           val3 <- NA
                                                         }
                                                         return(c(val, val2, val3))
                                                       }
                                              )
                                            val_iter_loc <- 1:length(val_iter)
                                            val_iter <- c(val_iter[val_iter_loc %% 3 == 1], 
                                                          val_iter[val_iter_loc %% 3 == 2],
                                                          val_iter[val_iter_loc %% 3 == 0])
                                            return(val_iter)
                                          }
                                 )
                               receptor_mat <- t(receptor_mat)
                               rownames(receptor_mat) <- sig_receptor_clu
                               colnames(receptor_mat) <- 
                                 c(get_pairs('P_val:', levels(condition), sep = ''),
                                   get_pairs('Intensity:', levels(condition), sep = ''),
                                   get_pairs('Intensity_nor:', levels(condition), sep = ''))
                               receptor_mat <- data.frame(receptor_mat, check.names = F)
                               if(!is.null(sig_ligand_clu)){
                                 data_iter <- rbind(data_iter, receptor_mat)
                                 ligand_clu <- 
                                   c(ligand_clu,rep(idents,length(sig_receptor_clu)))
                                 receptor_clu <- 
                                   c(receptor_clu,sig_receptor_clu)
                               }else{
                                 data_iter <- receptor_mat
                                 ligand_clu <- 
                                   c(rep(idents,length(sig_receptor_clu)))
                                 receptor_clu <- sig_receptor_clu
                               }
                             }  
                             
                             Interaction_name <- names(object@pnull_con)[x]  
                             Ligand_class <- 
                               unique(database[database$interaction_name %in%
                                                 Interaction_name, 'ligand_class'])
                             Receptor_class <- 
                               unique(database[database$interaction_name %in%
                                                 Interaction_name, 'receptor_class'])                               
                             if(!is.null(data_iter)){
                               data <- cbind(
                                 'Interaction_name' = names(object@pnull_con)[x],
                                 'Ligand' = L_R_pairs[1],
                                 'Receptor' = L_R_pairs[2],
                                 'Ligand_class' = Ligand_class,
                                 'Receptor_class' = Receptor_class,  
                                 'Ligand_cluster' = ligand_clu,
                                 'Receptor_cluster' = receptor_clu,
                                 data_iter
                               )
                             } 
                           })
                  result_con_iter <- do.call(rbind, result_con_iter)
                  if(class(result_con_iter)=='data.frame' && nrow(result_con_iter) > 0){
                    result_con_iter <- dplyr::distinct(result_con_iter)
                  }
                  return(result_con_iter)
                }) 
    names(result_con) <- levels(identity)
    return(result_con)
  }
}


get_DEGs <- function(
  object,
  matrix,
  group_var = NULL,
  ident.1 = NULL,
  cluster = '',
  pseudocount.use = 0
){
  #matrix <- as.matrix(matrix)
  cells.1 <- colnames(matrix)[group_var == ident.1]
  cells.2 <- setdiff(colnames(matrix), cells.1)
  features <- rownames(matrix)
  sum1 <- length(cells.1)
  sum2 <- length(cells.2)
  pct.1 <- round(rowSums(matrix[features, cells.1, drop = F] > 0) /sum1, digits = 3)
  pct.2 <- round(rowSums(matrix[features, cells.2, drop = F] > 0) /sum2, digits = 3)
  mean.fxn <- function(x) {return(log2(mean(x) + pseudocount.use))}
  data.1 <- apply(matrix[features, cells.1, drop = F], 1, FUN = mean.fxn)
  data.2 <- apply(matrix[features, cells.2, drop = F], 1, FUN = mean.fxn)
  data.3 <- 2^(data.1) - pseudocount.use
  data.4 <- 2^(data.2) - pseudocount.use
  avg.logFC <- data.1 - data.2
  p_val_adj <- if(cluster == ''){
    object@p_adj_gene_idents[[1]][features, ident.1]
  }else{
    object@p_adj_gene_con[[cluster]][features, ident.1]
  }
  df <- data.frame(
    avg_logFC = avg.logFC, 
    mean_umi_1 = data.3,
    mean_umi_2 = data.4,
    pct.1=pct.1,
    pct.2=pct.2,
    p_val_adj = p_val_adj,
    row.names = features
  )
  colnames(df) <- if(cluster == ''){
    paste0(ident.1, '_', colnames(df))
  }else{
    paste0(cluster, '_', ident.1, '_', colnames(df))
  }
  return(df)
}


get_DEGs_result <- function(
  object, 
  pseudocount.use = 0
){
  if(length(object@pnull_idents) > 0){
    matrix <- as.matrix(object@data)
    idents_exp <- 
      pbapply::pblapply(levels(object@identity), FUN = 
                          function(i){
                            result <- get_DEGs(
                              object = object, 
                              matrix = matrix, 
                              group_var = object@identity,
                              ident.1 = i, 
                              pseudocount.use = pseudocount.use)
                          })
    idents_exp <- do.call(cbind, idents_exp)
  }
  if(length(object@pnull_con) > 0){
    con_exp <- 
      pbapply::pblapply(levels(object@identity), FUN = 
                          function(i){
                            matrix <- as.matrix(object@data[, object@identity == i, drop = F])
                            group_var <- droplevels(object@condition[object@identity == i, drop = F])
                            if(nlevels(group_var) > 1){
                              result <- 
                                lapply(levels(object@condition), FUN = 
                                         function(j){
                                           result_iter <- get_DEGs(
                                             object = object, 
                                             matrix = matrix,
                                             group_var = object@condition[object@identity == i],
                                             cluster = i,
                                             ident.1 = j,
                                             pseudocount.use = pseudocount.use
                                           )
                                         })
                              result <- do.call(cbind, result)
                              return(result)
                            }
                          }
      )
    logi <- !sapply(con_exp, is.null)
    con_exp <- do.call(cbind, con_exp[logi])
  }
  return(list(idents_exp = idents_exp, con_exp = con_exp))
}
