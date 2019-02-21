#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pShaEP = args[1]

#pShaEP = "/home/aborrel/imitanib/results/MDanalysis/CHEMBL3617738_2hyy_MD/RMSDs/ligand/ligShaEP"
#namex = strsplit(pRMSD, "/")[[1]][7]

dShaEP = read.table(pShaEP, sep = "\t", header = TRUE)
rownames(dShaEP) = dShaEP[,1]



a = ggplot(dShaEP, aes(x = Time))+
  geom_line(aes(y = ESPscore, color = "ESP score"), size = 1.5)+
  geom_line(aes(y = Shape, color = "Shape score"), size = 1.5) +
  scale_fill_manual("Protein")+
  labs(y = "Scores scaled") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 15))
 
print(a)
ggsave(paste(pShaEP, ".png", sep=""), width = 12,height = 12, dpi = 300)




