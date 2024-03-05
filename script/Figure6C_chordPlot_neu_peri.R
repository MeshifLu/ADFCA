# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder", force = TRUE, update = FALSE)

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(circlize)
})

theme_set(theme_cowplot())
projectFolder='/data/tcl/Project/ADFCA'
setwd(projectFolder)
getwd()

############# Input files
projectName='ADFCA'
projectAnalysisFoler='Analysis'
dir.create(projectAnalysisFoler, showWarnings = FALSE)

figure6Folder <- paste(projectFolder, projectAnalysisFoler, 'Figure', 'Figure6', sep = '/')
figure6Folder

############# 

set.seed(999)

############# load matrix

neuToPer_Df = as.matrix(read.csv(paste0(figure6Folder, '/neuToPer_Df.csv'), row.names = 1))
neuToPer_Df_T = as.matrix(read.csv(paste0(figure6Folder, '/neuToPer_Df_T.csv'), row.names = 1))

#############  neuToPer_Df

#create a chord diagram but without labeling 
circos.par(start.degree = 270 ) # clock.wise = FALSE # track.margin=c(0,0)
# chordDiagram(as.matrix( neuToPer_Df), annotationTrack = "grid", preAllocateTracks = 1, big.gap = 30)
chordDiagram( neuToPer_Df, annotationTrackHeight = c(0.03, 0.01), annotationTrack = "grid", preAllocateTracks = 1, big.gap = 20)
# chordDiagram( neuToPer_Df, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))

#add the labels and axis
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  # print(xlim, ylim)
  # print(sector.name)
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 3, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.75) # adj => how close 
  
  #print axis
  circos.axis(h = "top", labels.cex = 0.35, major.tick.percentage = 0.5, major.at = seq(0,500,100), labels.facing = 'clockwise',
              sector.index = sector.name, track.index = 2)
}, bg.border = NA)


#saving the plot (high definition)
# dev.copy(png, paste0(figure6Folder, '/pic/neuToPer_Df.png'), width=16, height=16, units="in", res=600)
dev.copy(pdf, paste0(figure6Folder, '/pic/neuToPer_Df.pdf'), width=16, height=16, ) #units="in", res=500
dev.off()
circos.clear()


