#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pRMSD = args[1]

#pRMSD = "/home/aborrel/imitanib/results/MDanalysis/CHEMBL3617738_2hyy_MD/RMSDs/protein/protRMSD"

#namex = strsplit(pRMSD, "/")[[1]][7]

dRMSD = read.table(pRMSD, sep = "\t", header = TRUE)
rownames(dRMSD) = dRMSD[,1]



a = ggplot(dRMSD, aes(x = Time))+
  geom_line(aes(y = RMSDall, color = "RMSD All atoms"), size = 1)+
  geom_line(aes(y = RMSDC, color = "RMSD Calpha"), size = 1) +
  scale_fill_manual("Protein")+
  labs(y = "RMSD (Ang.)") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 15))
 
print(a)
ggsave(paste(pRMSD, "_RMSD.png", sep=""), width = 12,height = 12, dpi = 300)


b = ggplot(dRMSD, aes(x = Time))+
  geom_line(aes(y = Dmax, color = "Dmax"), size = 1.5) +
  scale_fill_manual("Protein")+
  labs(y = "Distance (Ang.)") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 15))

print(b)
ggsave(paste(pRMSD, "_Dmax.png", sep=""), width = 12,height = 12, dpi = 300)


