suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(circlize)
  library(ggrepel)
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


############# Plot 


# degAllE05Gene_Df = read.csv('Analysis/Figure/Figure3/degSummary/degLong_Df_genotype_compareToControl_cellN700.csv')
# degAllE05Gene_Df_tau = read.csv('Analysis/Figure/Figure3/degSummary/degAllE05Gene_Df_hTau.csv')

degN_neuron_Df = read.csv('Analysis/Figure/Figure6/degN_neuron_Df.csv')
# degN_neuron_Df_tau = 

degN_neuron_Df_tau = degN_neuron_Df %>% dplyr::filter(hTau >0 )
degN_neuron_Df_tau$gene = factor(degN_neuron_Df_tau$gene, levels = degN_neuron_Df_tau$gene)

pdf(file = file.path('Analysis/Figure/Figure6/pic/',
                     paste0('lineplot_deFequency_interval_hTau.pdf')),
    width=4.5, height=3)
p <- degN_neuron_Df_tau %>%
  ggplot(aes(x=gene, y=hTau )) + 
  geom_point()+
  geom_line()+
  xlab('Genes differentially expressed in hTau strain')+
  ylab('Number of cell types') +
  theme(axis.text.x=element_blank()) 
# p + geom_text(data = subset(degN_neuron_Df_tau, hTau>19), aes(gene,hTau,label=gene))

p+  geom_text_repel(data = subset(degN_neuron_Df_tau, hTau>19), aes(gene,hTau,label=gene), 
                    min.segment.length = 0, seed = 42, box.padding = 0.5) +
  geom_point(color = ifelse(degN_neuron_Df_tau$hTau > 19, "red", "black"))
dev.off()
