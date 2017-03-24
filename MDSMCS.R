#!/usr/bin/env Rscript


MDSbyAff = function(matin, aff){
  
  matin = as.dist(matin)
  fit <- cmdscale(matin, eig=TRUE, k=2)
  
  lcolaff <- colorRampPalette(c("white", "red"))
  lcolaff = lcolaff[dim(daffinity)[1]]
  
  plot (fit$points[,1], fit$points[,2], xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, pch = 19, xlim = c(-1.3, 1.3), ylim = c(-1, 1))
  #text (fit$points[,1], fit$points[,2]+0.05, labels = num_desc,  cex = 1.8, col = color_desc, font = 12)
  #text (fit$points[descriptor_selected,1], fit$points[descriptor_selected,2]+0.02, labels = "*",  cex = 4, col = color_desc)
  #text (fit$points[,1], fit$points[,2], labels = name_descriptor,  cex = 2.6, col = color_desc)
  




}



###########
#  MAIN   #
###########

args = commandArgs(TRUE)
pmatrix = args[1]
paffinity = args[2]

pmatrix="/home/aborrel/imitanib/results/analysis/MCS/_tanimoto"
paffinity = "/home/aborrel/imitanib/results/analysis/MCS/_aff"

dmatrix = read.table(pmatrix, sep="\t", header=TRUE)
print(dmatrix)

daff = read.table(paffinity, sep="\t")
print(daff)

MDSbyAff(dmatrix, daff)