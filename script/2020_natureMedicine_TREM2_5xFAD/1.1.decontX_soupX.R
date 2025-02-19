suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(scater)
  library(celda)
  library(DropletUtils)
  # library(future)
  
  library(SoupX)
})

# plan("multiprocess", workers = 4)
# plan()

theme_set(theme_cowplot())
projectFolder='~/Project/ADFCA/'
setwd(projectFolder)
getwd()


############# Input files
projectName='ADFCA'
projectAnalysisFoler='Analysis'
dir.create(projectAnalysisFoler, showWarnings = FALSE)

nm2020_Folder = paste0(projectFolder, 'Analysis/Figure/Figure4/2020_natureMedicine_TREM2_5xFAD')
nm2020_Folder

rawDataFolder = '/data/tcl/dataset/mouseAtlas/2020_natureMedicine_TREM2_5xFAD/dataset'


crSummaryFoler=file.path(projectFolder, 'cellrangerSummary')
crMtxFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allMatrix'))
crMtxRawFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allMatrixRaw'))
crH5Folder=file.path(projectFolder, paste0('cellrangerSummary', '/allH5'))
crH5RawFolder=file.path(projectFolder, paste0('cellrangerSummary', '/allH5raw'))

dataFolder_L = dir(crMtxFolder, include.dirs = FALSE)
dataFolder_L


############# Output files
decontxFolderPath=paste0(projectAnalysisFoler, '/1.decontx')
dir.create(decontxFolderPath, showWarnings = FALSE)
decontX_outs=paste0(decontxFolderPath,'/decontX_outs')
dir.create(decontX_outs, showWarnings = FALSE)
decontX_outs

soupxFolderPath=paste0(projectAnalysisFoler, '/1.soupx')
dir.create(soupxFolderPath, showWarnings = FALSE)
soupX_outFolder=file.path(soupxFolderPath, 'soupX_outs')
dir.create(soupX_outFolder, showWarnings = FALSE)
soupX_outFolder


############# Main pipeline

### 1. Use SoupX for removing ambient RNAs
# dataFolder=dataFolder_L[1]
for (dataFolder in dataFolder_L) {
  print(dataFolder)
  
  dataFolderPath=file.path(decontxFolderPath, dataFolder)
  dir.create(dataFolderPath, showWarnings = FALSE)
  picFolderPath=paste0(dataFolderPath, '/pic')
  dir.create(picFolderPath, showWarnings = FALSE) # Folder for the output pic
  
  soupDataFolderPath=file.path(soupxFolderPath, dataFolder)
  dir.create(soupDataFolderPath, showWarnings = FALSE)
  soupPicFolderPath=paste0(soupDataFolderPath, '/pic')
  dir.create(soupPicFolderPath, showWarnings = FALSE) # Folder for the output pic
  
  
  print(c('decontxFolderPath: ', decontxFolderPath))
  print(c('dataFolderPath: ', dataFolderPath))
  print(c('picFolderPath: ', picFolderPath))
  
  print(c('soupxFolderPath: ', soupxFolderPath))
  print(c('soupDataFolderPath: ', soupDataFolderPath))
  print(c('soupPicFolderPath: ', soupPicFolderPath))
  
  ##################################
  ### Input data into Seurat, without removing ambient RNAs
  print('Input data into Seurat')
  filt.matrix <- Read10X_h5(
    paste0(crH5Folder, '/', dataFolder, ".h5"),
    use.names = T)
  raw.matrix  <- Read10X_h5(
    paste0(crH5RawFolder, '/', dataFolder, ".h5"),
    use.names = T)
  # str(raw.matrix)
  # str(filt.matrix)

  srat  <- CreateSeuratObject(counts = filt.matrix)
  # Idents(srat) <- srat@meta.data$orig.ident

  ### Basic features
  print('Pre-decontX: basic cellQC')
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = '^mt:')

  pdf(file = file.path(picFolderPath,
                       paste0('nFeature_nCount_ptMt_', dataFolder, '.pdf')),
      width=10, height=5)
  print(VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()

  file.copy(from=file.path(picFolderPath, paste0('nFeature_nCount_ptMt_', dataFolder, '.pdf')),
            to=file.path(soupPicFolderPath, paste0('nFeature_nCount_ptMt_', dataFolder, '.pdf')),
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)

  ### Normalize data using sctransform
  srat    <- SCTransform(srat, verbose = F)

  # ### Normalize data using logNormalize
  # srat <- NormalizeData(srat)
  # srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  # srat <- ScaleData(srat, vars.to.regress = c('nCount_RNA'))

  
  # Idents(srat) <- srat@meta.data$seurat_clusters
  print('Pre-decontX: RunUMAP')
  srat    <- RunPCA(srat, npcs = 50, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
  srat    <- RunTSNE(srat, dims = 1:50, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:50, verbose = F)
  srat    <- FindClusters(srat, verbose = T)

  pdf(file = file.path(picFolderPath,
                       paste0('UMAP_tSNE_', dataFolder, '.pdf')),
      width=15, height=6)
  print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsne", label = TRUE))
  dev.off()

  file.copy(from=file.path(picFolderPath, paste0('UMAP_tSNE_', dataFolder, '.pdf')),
            to=file.path(soupPicFolderPath,  paste0('UMAP_tSNE_', dataFolder, '.pdf')),
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)


  ############################
  ### DecontX analysis

  print('DecontX analysis')
  sce <- read10xCounts(paste(crMtxFolder, dataFolder, sep = '/'))
  raw.sce <- read10xCounts(paste(crMtxRawFolder,  dataFolder, sep = '/'))

  sce.decontX <- decontX(x = sce, background = raw.sce)


  # Cluster labels on UMAP
  umap <- reducedDim(sce.decontX, "decontX_UMAP")

  pdf(file = file.path(picFolderPath,
                       paste0('decontX_contamination_', dataFolder, '.pdf')),
      width=15, height=6)
  print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) +
          plotDecontXContamination(sce.decontX))
  dev.off()


  decontX_outEach=paste0(decontX_outs, '/', dataFolder)
  # decontX_outEach
  if (dir.exists(decontX_outEach)) {
    unlink( decontX_outEach, recursive = TRUE ) #remove folder
  }
  DropletUtils:::write10xCounts(decontX_outEach, round(decontXcounts(sce.decontX)),
                                gene.symbol = rowData(sce)$Symbol,
                                barcodes = colData(sce)$Barcode)
  colData(sce)

  ### Post-decontX
  print('Post-decontX: cellQC')
  sratDecontx  <-
    Read10X(decontX_outEach) %>%
    CreateSeuratObject(project = strsplit(dataFolder, '_')[[1]][2], min.cells = 3, min.features = 200)

  sratDecontx[["percent.mt"]] <- Seurat::PercentageFeatureSet(sratDecontx, pattern = '^mt:')
  
  # sratDecontx <- NormalizeData(sratDecontx)
  # sratDecontx <- FindVariableFeatures(sratDecontx, selection.method = "vst", nfeatures = 2000)
  # sratDecontx <- ScaleData(sratDecontx, vars.to.regress = c('nCount_RNA'))
  
  ### sctransform
  sratDecontx    <- SCTransform(sratDecontx, verbose = F)
  # ### Normalize data using logNormalize
  # sratDecontx <- NormalizeData(sratDecontx)
  # sratDecontx <- FindVariableFeatures(sratDecontx, selection.method = "vst", nfeatures = 2000)
  # sratDecontx <- ScaleData(sratDecontx, vars.to.regress = c('nCount_RNA'))

  print('Post-decontX: RunUMAP')
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)

  pdf(file = file.path(picFolderPath,
                       paste0('UMAP_tSNE_', dataFolder, '_Decontx.pdf')),
      width=15, height=6)
  print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsne", label = TRUE) )
  dev.off()


  #########################################
  ### soupX

  # 1. Run suopX
  sc = SoupChannel(raw.matrix, filt.matrix)
  sc

  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- setDR(sc, umap)

  # head(meta)
  # rownames(srat@assays$RNA)


  pdf(file = file.path(soupPicFolderPath,
                       paste0(dataFolder, '_rho_Soupx.pdf')),
      width=8, height=6)
  print(sc  <- autoEstCont(sc, forceAccept = TRUE))
  dev.off()


  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20)
  adj.matrix  <- adjustCounts(sc, roundToInt = T)

  soupX_outData=file.path(soupX_outFolder, dataFolder)
  if (dir.exists(soupX_outData)) {
    unlink( soupX_outData, recursive = TRUE ) #remove folder
  }
  # adj.matrix
  # sc$toc
  DropletUtils:::write10xCounts( soupX_outData,
                                 adj.matrix)

  # 2. Plot soupX-filtered umap
  print('Post-SoupX: cellQC')
  sratSoupx  <-
    Read10X(soupX_outData) %>%
    CreateSeuratObject(project = dataFolder, min.cells = 3, min.features = 200)

  sratSoupx[["percent.mt"]] <- Seurat::PercentageFeatureSet(sratSoupx, pattern = '^mt:')
  print(VlnPlot(sratSoupx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)) #
  print('Post-SoupX: RunUMAP')

  ### sctransform
  sratSoupx    <- SCTransform(sratSoupx, verbose = F)
  # ### Normalize data using logNormalize
  # sratSoupx <- NormalizeData(sratSoupx)
  # sratSoupx <- FindVariableFeatures(sratSoupx, selection.method = "vst", nfeatures = 2000)
  # sratSoupx <- ScaleData(sratSoupx, vars.to.regress = c('nCount_RNA'))
  
  # # Idents(sratSoupx) <- sratSoupx@meta.data$seurat_clusters
  # print('Pre-decontX: RunUMAP')
  # sratSoupx    <- RunPCA(sratSoupx, npcs = 50, verbose = F)
  # sratSoupx    <- RunUMAP(sratSoupx, dims = 1:50, verbose = F)
  # sratSoupx    <- RunTSNE(sratSoupx, dims = 1:50, verbose = F)
  # sratSoupx    <- FindNeighbors(sratSoupx, dims = 1:50, verbose = F)
  # sratSoupx    <- FindClusters(sratSoupx, verbose = T)

  print('Post-soupX: RunUMAP')
  sratSoupx    <- RunPCA(sratSoupx, verbose = F)
  sratSoupx    <- RunUMAP(sratSoupx, dims = 1:50, verbose = F)
  sratSoupx    <- RunTSNE(sratSoupx, dims = 1:50, verbose = F)
  sratSoupx    <- FindNeighbors(sratSoupx, dims = 1:50, verbose = F)
  sratSoupx    <- FindClusters(sratSoupx, verbose = T)
  

  pdf(file = file.path(soupPicFolderPath,
                       paste0('UMAP_tSNE_', dataFolder, '_Soupx.pdf')),
      width=15, height=6)
  print( DimPlot(sratSoupx, reduction = "umap", label = TRUE) + DimPlot(sratSoupx, reduction = "tsne", label = TRUE) )
  dev.off()


  # sc = SoupChannel(raw.matrix, filt.matrix, calcSoupProfile = FALSE)
  # sc = estimateSoup(sc)
  #
  # sc = autoEstCont(sc)
  #
  #
  # sc = SoupChannel(raw.matrix, filt.matrix, calcSoupProfile = FALSE)
  # sc$



  rdataFile=file.path(decontxFolderPath,
                      paste('1', dataFolder, 'Decontx_Soupx.RData', sep = '.'))
  save.image(file = rdataFile )



}








