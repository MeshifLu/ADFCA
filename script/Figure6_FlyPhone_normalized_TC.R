set.seed(123)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(future.apply)
  library(Seurat)
  library(RColorBrewer)
  library(reshape2)
  library(network)
  library(igraph)
})

# cellN_n700 -> filter cell types
cellN700_Df <- read.csv('/data/tcl/Project/ADFCA/Analysis/Figure/Figure3/cellN/celltype_cellN_n700.csv')
cellN700_v <- cellN700_Df$annotation

# ### Ligand-receptor information
# LR_pairs <- read.csv(file = "/data/tcl/github/FlyPhoneDB_mod/annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt", sep = "\t")

### Ligand-receptor information with manual annotations
LR_pairs <- read.csv(file = "/data/tcl/github/FlyPhoneDB_mod/annotation/Ligand_receptor_pair_TCL_20231212.txt", sep = "\t")


####################################
# Input and cluster means
####################################

### For local test

# setwd('/data/tcl/github/FlyPhoneDB_mod')
plan(multiprocess, workers = 30)

# temp
# genotype='control_inTau'
# genotype='hTau'
# genotype='control_inAB42'

genotype_V = c('control_inTau', 'hTau', 'control_inAB42', 'AB42')
# # genotype_V = c('control_inTau', 'hTau')
# genotype_V = c('hTau')

for (genotype in genotype_V) {
  print(genotype)
}

for (genotype in genotype_V) {
  start_time <- Sys.time()
  
  ### create folders
  output_dir <- paste0('~/Project/ADFCA/Analysis/Figure/Figure6/', genotype)
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  if (!dir.exists(paste0(output_dir, "/heatmap"))) {dir.create(paste0(output_dir, "/heatmap"), recursive = TRUE)}
  if (!dir.exists(paste0(output_dir, "/dotplot"))) {dir.create(paste0(output_dir, "/dotplot"), recursive = TRUE)}
  if (!dir.exists(paste0(output_dir, "/circleplot"))) {dir.create(paste0(output_dir, "/circleplot"), recursive = TRUE)}
  
  ### load serurat files
  seuratObj <- readRDS(paste0("~/Project/ADFCA/Analysis/Figure/Figure6/seuratObj_", genotype, ".rds"))
  seuratObj <- subset(seuratObj, subset = adfca_annotation %in% cellN700_v) # only use cell types with more than 700 nuclei for the analysis
  seuratObj
  
  ### load metadata
  metadata <- seuratObj@meta.data
  Idents(seuratObj) <- "adfca_annotation"
  
  ### load expression data
  exprMat <- GetAssayData(seuratObj, assay = "RNA", slot = "data")
  exprMat <- as.matrix(exprMat)
  ### modify metadata and save it to cellInfo
  cellInfo <- metadata
  cellInfo$celltype <- as.character(cellInfo$adfca_annotation) # celltype annotations are stored
  ### subset cells # not necessary?
  exprMat <- exprMat[ , row.names(cellInfo)]
  
  
  ### start real score
  gene_list <- unique(c(LR_pairs$Gene_secreted, LR_pairs$Gene_receptor))
  common_genes <- intersect(gene_list, row.names(exprMat))
  
  LR_pairs <- subset(LR_pairs, Gene_secreted %in% common_genes & Gene_receptor %in% common_genes)
  
  exprMat <- as.matrix(exprMat) # cause a memory problem => subset essential data
  exprMat <- t(exprMat)
  df_Ligand <- exprMat[ , unique(LR_pairs$Gene_secreted)] # expression of ligand genes
  df_Receptor <- exprMat[ , unique(LR_pairs$Gene_receptor)] # expression of receptor genes
  
  celltype_df_Ligand <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Ligand)
  celltype_df_Receptor <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Receptor)
  # celltype_df_Ligand[1:3, 1:3]
  # celltype_df_Receptor[1:3, 1:3]
  
  # average Ligand counts by each celltype
  df_group_by_celltype_Ligand <- celltype_df_Ligand %>%
    group_by(celltype) %>%
    summarise_all(mean) %>%
    as.data.frame()
  # df_group_by_celltype_Ligand[1:3, 1:3]
  row.names(df_group_by_celltype_Ligand) <- df_group_by_celltype_Ligand$celltype
  df_group_by_celltype_Ligand$celltype <- NULL
  df_group_by_celltype_Ligand <- t(df_group_by_celltype_Ligand)
  # df_group_by_celltype_Ligand[1:3, 1:3]
  # str(df_group_by_celltype_Ligand)
  # write.csv(df_group_by_celltype_Ligand,
  #           file = "../Data/2021-02-15_df_group_by_celltype_Ligand.csv")
  
  # average Receptor counts by each celltype
  df_group_by_celltype_Receptor <- celltype_df_Receptor %>%
    group_by(celltype) %>%
    summarise_all(mean) %>%
    as.data.frame()
  # df_group_by_celltype_Receptor[1:3, 1:3]
  row.names(df_group_by_celltype_Receptor) <- df_group_by_celltype_Receptor$celltype
  df_group_by_celltype_Receptor$celltype <- NULL
  df_group_by_celltype_Receptor <- t(df_group_by_celltype_Receptor)
  # df_group_by_celltype_Receptor[1:3, 1:3]
  # str(df_group_by_celltype_Receptor)
  # write.csv(df_group_by_celltype_Receptor,
  #           file = "../Data/2021-02-15_df_group_by_celltype_Receptor.csv")
  
  ####################################
  # Interaction score
  ####################################
  
  ligand_avg <- df_group_by_celltype_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
  receptor_avg <- df_group_by_celltype_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()
  # write.csv(ligand_avg, file = "../Data/2021-02-15_ligand_avg.csv")
  # write.csv(receptor_avg, file = "../Data/2021-02-15_receptor_avg.csv")
  
  # x <- sort(unique(cellInfo$celltype)) %>% as.data.frame()
  
  interaction_list <- list()
  LR_pairs_one <- LR_pairs # combine
  
  for (i in sort(unique(cellInfo$celltype)) ) {
    # for (i in c("aEC1", "aEC2")) {
    # print("")
    # print(i)
    # print("")
    LR_pairs_combine <- LR_pairs # combine
    for (j in sort(unique(cellInfo$celltype)) ) {
      # for (j in c("aEC3", "aEC4") ) {
      print(paste0(i, ">", j, ', Sys.time: ', Sys.time()))
      LR_pairs_tmp <- LR_pairs
      LR_pairs_tmp[[paste0(i, ">", j, "_score")]] <- log1p(ligand_avg[[i]]) * log1p(receptor_avg[[j]])
      
      # permutation -------------------------------------------------------------
      # start permutatioin
      permutation_times <- 1000
      y <- future_lapply(1:permutation_times, function(ii) {
        # for (ii in 1:permutation_times) {
        #   LR_pairs_tmp[[paste0("permute", ii)]] <- local({
        # print(i)
        # sample Ligand
        cellInfo_sample_Ligand <- cellInfo
        # str(cellInfo_sample_Ligand)
        cellInfo_sample_Ligand$celltype <- sample(cellInfo_sample_Ligand$celltype)
        
        celltype_df_sample_Ligand <- cbind(cellInfo_sample_Ligand[ , c("celltype"), drop = FALSE], df_Ligand)
        
        df_group_by_celltype_sample_Ligand <- celltype_df_sample_Ligand %>%
          group_by(celltype) %>%
          summarise_all(mean) %>%
          as.data.frame()
        # str(df_group_by_celltype_sample_Ligand)
        row.names(df_group_by_celltype_sample_Ligand) <- df_group_by_celltype_sample_Ligand$celltype
        df_group_by_celltype_sample_Ligand$celltype <- NULL
        df_group_by_celltype_sample_Ligand <- t(df_group_by_celltype_sample_Ligand)
        # str(df_group_by_celltype_sample_Ligand)
        
        # sample Receptor
        cellInfo_sample_Receptor <- cellInfo
        cellInfo_sample_Receptor$celltype <- sample(cellInfo_sample_Receptor$celltype)
        
        celltype_df_sample_Receptor <- cbind(cellInfo_sample_Receptor[ , c("celltype"), drop = FALSE], df_Receptor)
        
        df_group_by_celltype_sample_Receptor <- celltype_df_sample_Receptor %>%
          group_by(celltype) %>%
          summarise_all(mean) %>%
          as.data.frame()
        # str(df_group_by_celltype_sample_Receptor)
        row.names(df_group_by_celltype_sample_Receptor) <- df_group_by_celltype_sample_Receptor$celltype
        df_group_by_celltype_sample_Receptor$celltype <- NULL
        df_group_by_celltype_sample_Receptor <- t(df_group_by_celltype_sample_Receptor)
        # df_group_by_celltype_sample_Receptor[1:2, 1:2]
        # str(df_group_by_celltype_sample_Receptor)
        
        ####################################
        # Interaction score
        ####################################
        ligand_avg_tmp <- df_group_by_celltype_sample_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
        receptor_avg_tmp <- df_group_by_celltype_sample_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()
        # LR_pairs_tmp <- LR_pairs
        # score <- ligand_avg[[paste0("Ligand_cluster", i)]] * receptor_avg[[paste0("Receptor_cluster", j)]]
        # colnames(LR_pairs_tmp)[ncol(LR_pairs_tmp)] <- paste0("permute", i)
        tmp <- log1p(ligand_avg_tmp[[i]]) * log1p(receptor_avg_tmp[[j]])
        tmp
      }, future.seed = TRUE)
      # }
      
      df <- data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE))
      df <- t(df)
      LR_pairs_tmp <- cbind(LR_pairs_tmp, df)
      # dim(LR_pairs_tmp)
      # dim(df)
      
      LR_pairs_tmp$result <- rowSums(sapply(LR_pairs_tmp[, 13:ncol(LR_pairs_tmp)], function(x) x > LR_pairs_tmp[[paste0(i, ">", j, "_score")]]))
      # head(LR_pairs_tmp[ , c("PM1_nonhemo_interaction_score", "result")])
      LR_pairs_tmp[[paste0(i, ">", j, "_pvalues")]] <- LR_pairs_tmp$result / permutation_times
      # write.csv(LR_pairs_tmp, file = "LR_pairs_tmp.csv")
      
      LR_pairs_tmp <- LR_pairs_tmp[ , c(1:12, ncol(LR_pairs_tmp))]
      # LR_pairs_tmp[LR_pairs_tmp$PM1_nonhemo_interaction_score]
      LR_pairs_tmp[LR_pairs_tmp[[paste0(i, ">", j, "_score")]] == 0, paste0(i, ">", j, "_pvalues")] <- 1
      
      LR_pairs_combine <- cbind(LR_pairs_combine,
                                LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
      )
      LR_pairs_one <- cbind(LR_pairs_one,
                            LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
      )
    }
    
    interaction_list[[i]] <- LR_pairs_combine
    
  }
  
  write.csv(LR_pairs_one, file = paste0(output_dir, "/", genotype, "_interaction_list_TC.csv"))
  
  end_time <- Sys.time()
  total_time <- as.character(end_time - start_time)
  
  time_Df <- data.frame(timeGroup = c('startTime', 'endTime', 'totalTime'), 
                        time = c(as.character(start_time), as.character(end_time), total_time ))
  write.csv(time_Df, file = paste0(output_dir, "/", genotype, "_runningTime_TC.csv"), row.names = FALSE)
  
}

