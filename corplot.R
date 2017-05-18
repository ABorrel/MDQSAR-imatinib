#!/usr/bin/env Rscript

library(ggplot2)

plotCor = function(din, vblack){
  
  cor.val = cor(din[,1], din[,2])
  p.val = cor.test (din[,1], din[,2])$p.value
  #print (cor.test (d_temp[,2], d_temp[,3]))
  #print (p.val)
  
  
  # coloration point considering the affinity - colounm 3
  #din = din[order(din[,3], decreasing = TRUE),]
  #colpointramp <- colorRampPalette(c("lightgreen", "red"))
  #colpoint = colpointramp(dim(din)[1])


  pIC50 = din$Aff
  qplot(Dock_score,emodel,data=din) +
  geom_point(size=2, stroke = 0.5)+
  geom_point(aes(colour = pIC50))+
  scale_colour_gradient(low="red", high="green")
  ggsave(paste(p_filin, "_ggplot.png",sep=""),width = 8,height = 7)

  legend_title="-log(Ki, Kd, IC50)"

  ggplot(din, aes(Dock_score,emodel, label = rownames(din))) +
  geom_text(size = 2, aes(colour = pIC50)) +
    #scale_fill_manual(legend_title)+
    scale_colour_gradient(low="red", high="green")
    #theme(legend.title = "-log(Ki, Kd, IC50)")
  ggsave(paste(p_filin, "text_ggplot.png",sep=""),width = 8,height = 7)


  ASSAY_CHEMBLID = din$ASSAY_CHEMBLID
  ggplot(din, aes(Dock_score,emodel)) +
  geom_point(size = 1, aes( colour = ASSAY_CHEMBLID))+
  theme(legend.position="none")
  #scale_colour_gradient(low="red", high="green")
  ggsave(paste(p_filin, "BA_ggplot.png",sep=""),width = 8,height = 7)
  
  
  ggplot(din, aes(Dock_score,emodel)) +
    geom_text(size = 0.70, aes( colour = ASSAY_CHEMBLID, label = ASSAY_CHEMBLID))+
    theme(legend.position="none")
  #scale_colour_gradient(low="red", high="green")
  ggsave(paste(p_filin, "BAtxt_ggplot.png",sep=""),width = 8,height = 7)
  

  # png
  #png(paste(p_filin, ".png", sep = ""), width = 1200, height = 1200, res = 150)
  #layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), c(4,1), c(4,0.1), TRUE)
  #plot(din[,1], din[,2], main = paste("cor =", round(cor.val,3), "p-val =",round(p.val,6) , "dim =", dim(din)[1], sep = " "), xlab = colnames (din)[1], ylab = colnames (din)[2], pch = 19, col = colpoint)
  #abline(lm (din[,2]~ din[,1]), col = "red")
  
  #iblack = which(rownames(din) == vblack)
  #points(din[iblack,1], din[iblack,2])#, vblack)
  #text(din[iblack,1], din[iblack,2]-10, vblack)
  
  #legend_image <- as.raster(matrix(colpointramp(dim(din)[1]), ncol=1))
  #par(mar = c(0,0,1,0))
  #plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '-log10 Affinity')
  #text(x=1.5, y = seq(0,1,l=5), labels = seq(round(min(din[,3]),1),round(max(din[,3])),l=5))
  #rasterImage(legend_image, -0.5, -0.5, 1,1)
  
  #dev.off ()


  # with text
  #png(paste(p_filin, "text.png", sep = ""), width = 2200, height = 2200, res = 350)
  #layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), c(4,1), c(4,0.1), TRUE)
  #plot(din[,1], din[,2], main = paste("cor =", round(cor.val,3), "p-val =",round(p.val,6) , "dim =", dim(din)[1], sep = " "), xlab = colnames (din)[1], ylab = colnames (din)[2], pch = 19, col = colpoint, type = "n")

  #text(din[,1], din[,2], rownames(din), cex = 0.3, col = colpoint)

  #abline(lm (din[,2]~ din[,1]), col = "red")

  #iblack = which(rownames(din) == vblack)
  #points(din[iblack,1], din[iblack,2])#, vblack)
  #text(din[iblack,1], din[iblack,2], vblack)

  #legend_image <- as.raster(matrix(colpointramp(dim(din)[1]), ncol=1))
  #par(mar = c(0,0,1,0))
  #plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '-log10 Affinity')
  #text(x=1.5, y = seq(0,1,l=5), labels = seq(round(min(din[,3]),1),round(max(din[,3])),l=5))
  #rasterImage(legend_image, -0.5, -0.5, 1,1)

  #dev.off ()




}


############
##  MAIN  ##
############

args <- commandArgs(TRUE)
p_filin = args[1]
pchembl = args[2]

#p_filin = "/home/aborrel/imitanib/results/analysis/dockingXPAll/ScoreVSAff.txt"
#pchembl = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_allAff.txt"

din = read.csv (p_filin, header = TRUE, sep = "\t")
rownames(din) = din[,1]
din = din[,-1]

# extract bioassay
dchembl = read.csv(pchembl, header = TRUE, sep = "\t")
rownames(dchembl) = dchembl$CMPD_CHEMBLID

ASSAY_CHEMBLID = dchembl[rownames(din), c("ASSAY_CHEMBLID")]

din = cbind(din, ASSAY_CHEMBLID)
# glevec
plotCor(din, c("CHEMBL941"))



#biofilm_3IX4_XP=read.csv(p_filin,header=T,stringsAsFactors = F, sep = "\t")
#pIC50= biofilm_3IX4_XP$Aff
#qplot(Dock_score,emodel,data=biofilm_3IX4_XP) +
#  geom_point(colour="black",size=4)+
#  geom_point(aes(colour = pIC50))+
#  scale_colour_gradient(low="red", high="green")
#ggsave(paste("file_of_your_image.png",sep=""),width = 8,height = 7)

