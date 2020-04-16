#!/usr/bin/env Rscript
library(ggplot2)
library(grid)
library(gridExtra)

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_prot = args[1]
p_lig = args[2]
p_res = args[3]
chemid = args[4]
pr_out = args[5]


#p_prot = "./../../results/MD-out/CHEMBL1173498_2hyy_MD/RMSDs/protein/protRMSD"
#p_lig = "./../../results/MD-out/CHEMBL1173498_2hyy_MD/RMSDs/ligand/ligRMSD"
#p_res = "./../../results/MD-out/CHEMBL1173498_2hyy_MD/RMSDs/residues/resRMSD_BS"
#pr_out = "./../../results/MD-qualityCheck/IC50/"  
#chemid = "CHEMBL1173498"


# plot RMSD prot
dRMSDprot = read.csv(p_prot, sep = "\t")
a = ggplot(dRMSDprot, aes(x = Time))+
  geom_line(aes(y = RMSDall, color = "All atoms"), size = 0.3)+
  geom_line(aes(y = RMSDC, color = "C-alpha"), size = 0.3) +
  scale_fill_manual("Protein")+
  labs(y = "RMSD (Ang.)", x = "Time (ns)") +
  ggtitle("RMSD protein") +
  theme(axis.text.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.x =  element_text(size = 10, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 10))



# for RMSD lig
dRMSDlig = read.csv(p_lig, sep = "\t")
b = ggplot(dRMSDlig, aes(x = Time))+
  geom_line(aes(y = RMSD), size = 0.3, colour="#FF9999")+
  scale_fill_manual("Protein")+
  labs(y = "RMSD (Ang.)", x = "Time (ns)") +
  ggtitle("RMSD ligand") +
  theme(axis.text.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.x =  element_text(size = 10, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 10))




# RMSF for res
dRMSF = read.csv(p_res, sep = "\t")
rownames(dRMSF) = dRMSF[,1]

lBS = as.vector(as.double(rownames(dRMSF)[which(dRMSF$BS == 1)]))

fres = ggplot(dRMSF, aes(x = NameRes))+
  geom_line(aes(y = all, color = "All atoms"), size = 0.3)+
  geom_line(aes(y = Ca, color = "C-alpha"), size = 0.3) +
  scale_fill_manual("Protein")+
  labs(y = "RMSF (Ang.)", x="Residues") +
  ggtitle("RMSF residues") +
  geom_vline(xintercept = lBS, color=t_col("black", 80), size=1, linetype="solid")+
  theme(axis.text.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.y = element_text(size = 10, hjust = 0.5, vjust =0.1), 
        axis.title.x =  element_text(size = 10, hjust = 0.5, vjust =0.1), 
        legend.title = element_text(size = 0), legend.text = element_text (size= 10))



g = grid.arrange(
  name = chemid,
  a,b,fres,
  layout_matrix = rbind(c(1, 2),
                        c(3, 3)),
  top = textGrob(chemid,gp=gpar(fontsize=15,font=3)))

ggsave(paste(pr_out,chemid,".png", sep = ""), plot=g, device = "png", dpi = 300, units = "cm", width = 21, height = 15)


#ggsave(paste(pr_out,chemid,".pdf", sep = ""), plot=g, device = "pdf", dpi = 300, units = "cm", width = 21, height = 15)


