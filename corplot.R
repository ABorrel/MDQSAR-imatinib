#!/usr/bin/env Rscript



plotCor = function(din){
  
  cor.val = cor(din[,1], din[,2])
  p.val = cor.test (din[,1], din[,2])$p.value
  #print (cor.test (d_temp[,2], d_temp[,3]))
  #print (p.val)
  
  # png
  png (paste(p_filin, ".png", sep = ""), width = 1200, height = 1200, res = 150)
  plot (din[,1], din[,2], main = paste("cor =", cor.val, "p-val =",p.val , "dim =", dim(din)[1], sep = " "), xlab = colnames (din)[1], ylab = colnames (din)[2], pch = 19)
  abline (lm (din[,2]~ din[,1]), col = "red")
  dev.off ()
  
}


############
##  MAIN  ##
############

args <- commandArgs(TRUE)
p_filin = args[1]

p_filin = "/home/aborrel/imitanib/results/analysis/docking/ScoreVSAff.txt"

din = read.csv (p_filin, header = TRUE, sep = "\t")
rownames(din) = din[,1]
din = din[,-1]

plotCor(din)
