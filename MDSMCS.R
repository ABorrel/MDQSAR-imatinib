#!/usr/bin/env Rscript


MDSbyAff = function(matin, aff, pout){
  
  matin = as.dist(matin)
  aff = as.matrix(aff)
  fit <- cmdscale(matin, eig=TRUE, k=2)
  
  #print (fit)
  aff = aff[order(aff, decreasing = TRUE),]
  dataplot = fit$points
  dataplot = dataplot[names(aff),]
  
  lcolaff <- colorRampPalette(c("red", "lightgreen"))
  lcolaff = lcolaff(dim(fit$points)[1])
  
  png(paste(pout, "MDSplot.png"), 800, 800)
  plot (dataplot[,1], dataplot[,2], xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, pch = 19, col = lcolaff)
  dev.off()
  
  png(paste(pout, "MDSplotText.png"), 800, 800)
  plot (dataplot[,1], dataplot[,2], xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n")
  text (dataplot[,1], dataplot[,2], labels = rownames(dataplot),  cex = 1.8, font = 12, col = lcolaff)
  dev.off()
  
}



###########
#  MAIN   #
###########

args = commandArgs(TRUE)
pmatrix = args[1]
paffinity = args[2]

pmatrix="/home/aborrel/imitanib/results/analysis/MCS/tanimoto"
paffinity = "/home/aborrel/imitanib/results/analysis/MCS/aff"

dmatrix = read.table(pmatrix, sep="\t", header=TRUE)
daff = read.table(paffinity, sep="\t", header=TRUE, stringsAsFactors = F)

MDSbyAff(dmatrix, daff, pmatrix)