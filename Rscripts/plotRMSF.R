#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pRMSF = args[1]

#pRMSF = "/home/aborrel/imitanib/results/MDanalysis/CHEMBL3617738_2hyy_MD/RMSDs/ligand/ligRMSF"
#pRMSF = "/home/aborrel/imitanib/results/MDanalysis/CHEMBL3617738_2hyy_MD/RMSDs/residues/resRMSD"


dRMSD = read.table(pRMSF, sep = "\t", header = TRUE)
rownames(dRMSD) = dRMSD[,1]



a = ggplot(dRMSD, aes(x = NameRes))+
  geom_line(aes(y = all, color = "RMSF All atoms"), size = 1)+
  geom_line(aes(y = Ca, color = "RMSF Calpha"), size = 1) +
  scale_fill_manual("Protein")+
  labs(y = "RMSF (Ang.)") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 15))

print(a)
ggsave(paste(pRMSF, "_RMSF.png", sep=""), width = 30,height = 12, dpi = 300)


b = ggplot(dRMSD, aes(x = NameRes))+
  geom_line(aes(y = Dmax, color = "Dmax")) +
  scale_fill_manual("Protein")+
  labs(y = "Distance (Ang.)") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 15))

print(b)
ggsave(paste(pRMSF, "_Dmax.png", sep=""), width = 30,height = 12, dpi = 300)


