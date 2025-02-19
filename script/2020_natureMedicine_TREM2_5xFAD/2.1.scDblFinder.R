# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder", force = TRUE, update = FALSE)

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(scater)

  library(DropletUtils)
  library(scDblFinder)
  library(BiocParallel)
  library(tidyseurat)
})

theme_set(theme_cowplot())
projectFolder='~/Project/ADFCA/'
setwd(projectFolder)
getwd()

############# Input files
projectName='ADFCA'
projectAnalysisFoler='Analysis'
dir.create(projectAnalysisFoler, showWarnings = FALSE)

subprojectFolder = paste0(projectFolder, 'Analysis/Figure/Figure7/2020_natureMedicine_TREM2_5xFAD')
subprojectFolder

rawDataFolder = '/data/tcl/dataset/mouseAtlas/2020_natureMedicine_TREM2_5xFAD/dataset'




# crSummaryFoler=file.path(projectFolder, 'cellrangerSummary')
# crMtxFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allMatrix'))
# crMtxRawFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allMatrixRaw'))
# crH5Folder=file.path(projectFolder, paste0('cellrangerSummary', '/allH5'))
# crH5RawFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allH5raw'))

scdblfinderFoler= paste(subprojectFolder, '2.scDblFinder', sep = '/')
dir.create(scdblfinderFoler, showWarnings = FALSE)
metadataFolder <- paste0(scdblfinderFoler, '/metadata')
dir.create(metadataFolder, showWarnings = FALSE)
picFolderPath <- file.path(scdblfinderFoler, 'pic')
dir.create(picFolderPath, showWarnings = FALSE)

mtxFile_L = list.files(rawDataFolder, pattern = '*_matrix.mtx.gz')
mtxFile_L

# mtxFile="GSM4160643_WT_Cor_matrix.mtx.gz"
# fileName=strsplit(mtxFile, '_matrix.mtx.gz')[[1]]
# print(fileName)
# 
# obs=paste0(rawDataFolder, '/', fileName, '_barcodes.tsv.gz')
# varDf=paste0(rawDataFolder, '/', fileName, '_genes.tsv.gz')
# 
# expression_matrix <- ReadMtx(mtx = paste0(rawDataFolder, '/', mtxFile), cells = obs, features = varDf)
# seurat_object <- CreateSeuratObject(counts = expression_matrix)

############## 
set.seed(123)

for (mtxFile in mtxFile_L) {
  print(mtxFile)
  fileName=strsplit(mtxFile, '_matrix.mtx.gz')[[1]]
  print(fileName)

  print('Input data into Seurat')
  obs=paste0(rawDataFolder, '/', fileName, '_barcodes.tsv.gz')
  varDf=paste0(rawDataFolder, '/', fileName, '_genes.tsv.gz')
  
  expression_matrix <- ReadMtx(mtx = paste0(rawDataFolder, '/', mtxFile), cells = obs, features = varDf)
  srat <- CreateSeuratObject(counts = expression_matrix)
  
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = '^mt-')
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat <- ScaleData(srat, vars.to.regress = c('nCount_RNA'))


  print('Pre-decontX: RunUMAP')
  srat    <- RunPCA(srat, npcs = 50, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:50, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  srat    <- RunUMAP(srat, dims = 1:50, verbose = F)


  print('Calculate doublet rate')
  # Calculate doublet rate
  srat
  doubletRate = 0.01*((ncol(srat))*0.0008 + 0.0527 )
  doubletRate
  doubletCellN=as.integer( ncol(srat) * doubletRate )
  print(paste0('Doublet rates: ', doubletRate))

  print('Transfer seurat data to sce')
  # Transfer seurat data to sce
  sce <- scDblFinder(GetAssayData(srat, slot="counts"), clusters=Idents(srat), dbr=doubletRate, BPPARAM=MulticoreParam(6))
  # port the resulting scores back to the Seurat object:
  srat$scDblFinder.score <- sce$scDblFinder.score
  srat$scDblFinder.class <- sce$scDblFinder.class
  table(sce$scDblFinder.class)


  p1 = FeaturePlot(object = srat, features = c('scDblFinder.score'))
  p2 = UMAPPlot(srat, reduction = 'umap', group.by = 'scDblFinder.class' )
  pdf(file = file.path(picFolderPath,
                       paste0('UMAP_scDblFinder_', fileName, '.pdf')),
      width=10, height=5)
  print(p1 | p2)
  dev.off()

  metadataOutPath=paste0(metadataFolder, '/', fileName, '.csv')
  write.csv(srat@meta.data, file = metadataOutPath, )

}







# head(srat@meta.data)
# 
# srat@assays$RNA[1:10, 1:50]
# srat@assays$RNA@counts[1:10, 1:50]
# 
# Idents(srat) <- "scDblFinder.class"
# DimPlot(srat, reduction = "umap")
# DimPlot(data.filt, group.by = "scDblFinder.score")
# UMAPPlot(srat, group.by = 'scDblFinder.score')
# DimPlot(srat, reduction = "umap", split.by = 'scDblFinder.class')
# FeaturePlot(pbmc, features = c('scDblFinder.score', ))
# 
# FeaturePlot(object = srat, features = c('scDblFinder.score'))
# UMAPPlot(srat, reduction = 'umap', group.by = 'scDblFinder.class' )
# 
# sce <- scDblFinder(sce, samples="sample_id", BPPARAM=MulticoreParam(3))
# 
# 
# 
# 
# 
# sort(srat$scDblFinder.score, decreasing = TRUE)
# 
# 
# metadata(sce)
# colData(sce)

