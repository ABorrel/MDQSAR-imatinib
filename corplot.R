#!/usr/bin/env Rscript

library(ggplot2)

plotCor = function(din, vblack){
  
  cor.val = cor(din[,1], din[,2])
  p.val = cor.test (din[,1], din[,2])$p.value
  #print (cor.test (d_temp[,2], d_temp[,3]))
  #print (p.val)
  
  
  # coloration point considering the affinity - colounm 3
  din = din[order(din[,3], decreasing = TRUE),]
  colpointramp <- colorRampPalette(c("lightgreen", "red"))
  colpoint = colpointramp(dim(din)[1])
  
  # png
  png(paste(p_filin, ".png", sep = ""), width = 1200, height = 1200, res = 150)
  layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), c(4,1), c(4,0.1), TRUE)
  plot(din[,1], din[,2], main = paste("cor =", round(cor.val,3), "p-val =",round(p.val,6) , "dim =", dim(din)[1], sep = " "), xlab = colnames (din)[1], ylab = colnames (din)[2], pch = 19, col = colpoint)
  abline(lm (din[,2]~ din[,1]), col = "red")
  
  iblack = which(rownames(din) == vblack)
  points(din[iblack,1], din[iblack,2])#, vblack)
  text(din[iblack,1], din[iblack,2], vblack)
  
  legend_image <- as.raster(matrix(colpointramp(dim(din)[1]), ncol=1))
  par(mar = c(0,0,1,0))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log10 Affinity')
  text(x=1.5, y = seq(0,1,l=5), labels = seq(round(min(din[,3]),1),round(max(din[,3])),l=5))
  rasterImage(legend_image, -0.5, -0.5, 1,1)
  
  dev.off ()
  
}


############
##  MAIN  ##
############

args <- commandArgs(TRUE)
p_filin = args[1]

#p_filin = "~/imitanib/results/analysis/docking/ScoreVSAff.txt"

din = read.csv (p_filin, header = TRUE, sep = "\t")
rownames(din) = din[,1]
din = din[,-1]

# glevec
plotCor(din, c("CHEMBL941"))



#biofilm_3IX4_XP=read.csv(p_filin,header=T,stringsAsFactors = F, sep = "\t")
#pIC50= biofilm_3IX4_XP$Aff
#qplot(Dock_score,emodel,data=biofilm_3IX4_XP) +
#  geom_point(colour="black",size=4)+
#  geom_point(aes(colour = pIC50))+
#  scale_colour_gradient(low="red", high="green")
#ggsave(paste("file_of_your_image.png",sep=""),width = 8,height = 7)

