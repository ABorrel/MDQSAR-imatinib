#!/usr/bin/env Rscript
library(factoextra)
library(ggplot2)
require(cluster)

optimalCluters = function (din, prout, met){
  print (din)
  
  p = fviz_nbclust(din, hcut, method = met, k.max = 50)
  print(p)
}




