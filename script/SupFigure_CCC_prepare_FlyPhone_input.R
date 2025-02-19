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
  library(dior)
})

setwd('~/Project/ADFCA/')

### Prepare folder

figureCCCFolder='Analysis/Figure/FigureCCC'
dir.create(figureCCCFolder, showWarnings = FALSE)


### Load h5 file, (h5 file was generated using scDIOR)
seuratObj <- dior::read_h5("adataProcess/v0.6/adfca_headbody_v0.6_woHarmony_raw.h5")
seuratObj

metadata <- seuratObj@meta.data
seurat_concat <- CreateSeuratObject(counts = seuratObj@assays$umi_counts, project = "ADFCA")
seurat_concat@meta.data <- metadata
seurat_concat@assays$RNA@data <- seuratObj@assays$RNA@data

# seurat_concat@assays$RNA@counts[1:10,1:10]
# seurat_concat@assays$RNA@data[1:10,1:10]
# seurat_concat

# control
seurat_gt <- subset(x = seurat_concat, subset = genotype  == 'control')
seurat_gt
saveRDS(seurat_gt, file = paste0(figureCCCFolder, '/seuratObj_', 'control', '.rds'))
seurat_gt <- subset(x = seurat_concat, subset = genotype  == 'control' & age %in% c('10', '20'))
seurat_gt
saveRDS(seurat_gt, file = paste0(figureCCCFolder, '/seuratObj_', 'control_inAB42', '.rds'))
seurat_gt <- subset(x = seurat_concat, subset = genotype  == 'control' & age %in% c('20', '30'))
seurat_gt
saveRDS(seurat_gt, file = paste0(figureCCCFolder, '/seuratObj_', 'control_inTau', '.rds'))

# AB42
seurat_gt <- subset(x = seurat_concat, subset = genotype  == 'AB42')
saveRDS(seurat_gt, file = paste0(figureCCCFolder, '/seuratObj_', 'AB42', '.rds'))
# hTau
seurat_gt <- subset(x = seurat_concat, subset = genotype  == 'hTau')
saveRDS(seurat_gt, file = paste0(figureCCCFolder, '/seuratObj_', 'hTau', '.rds'))

