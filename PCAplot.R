#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
library(fastICA)

generatePCAcoords = function(din){
  
  dinScale = scale(din)
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp

  return(list(data_plot, var_cap[1]))
  
  
}


generateIPAcoords = function(din, path_result){
  
  a = fastICA(din, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)
  print (a)
  
  gifGeneration(paste(path_result, "ICA3D", sep = ""), a$S)
}


PCAcombined3plans = function(d1D, d2D, d3D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord2D = generatePCAcoords(d2D)
  lcoord3D = generatePCAcoords(d3D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord2D[[1]][,1])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  print(paste(lcoord1D[[2]], lcoord2D[[2]], lcoord3D[[2]], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D", sep = ""), coordSpace)
  
  
}


PCAcombined2plans = function(d1D, d2D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord3D = generatePCAcoords(d2D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord1D[[1]][,2])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  #print(paste(lcoord1D[[2]], lcoord2D[[2]], lcoord3D[[2]], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D2plan", sep = ""), coordSpace)
  
  
}



PCA3D = function(din, path_result){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  #col.desc = colorDesc(colnames(din))
  
  col.desc = "black"
  
  gifGeneration(paste(path_result, "PCA3D", sep = ""), data_plot)
  
}


PCAplot = function (din, path_result){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  #col.desc = colorDesc(colnames(din))
  
  col.desc = "black"
  
  print(col.desc)
  
  png (paste (path_result, "_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  colorplot = rep("black", dim(data_plot)[1])
  names(colorplot) = rownames(data_plot)
  colorplot["CHEMBL941"] = "red"
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2, col = colorplot)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  colpoint <- colorRampPalette(c("red", "lightgreen"))
  
  png (paste (path_result, "_color.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = colpoint(dim(data_plot)[1]), xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (path_result, "_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
  
  svg (file = paste (path_result, "_descriptor.svg", sep = ""), 25, 25)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 3 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 3.5)
  dev.off()
  
}


