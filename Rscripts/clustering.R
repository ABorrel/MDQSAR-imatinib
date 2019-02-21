#!/usr/bin/env Rscript
library(factoextra)
library(ggplot2)
require(cluster)
library(RootsExtremaInflections)
library(mclust)
library(fossil)

optimalCluters = function (din, prout, metcluster, metOptNB, metagregation){
  
  
  if (metcluster == "hclust"){
    p = fviz_nbclust(din, hcut, hcut_metho = metagregation, method = metOptNB, k.max = 50)
    ggsave(paste(prout, metcluster, "_" , metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)
  }else if(metcluster == "kmeans"){
    p = fviz_nbclust(din, kmeans, method = metOptNB, k.max = 50)
    ggsave(paste(prout, metcluster, "_" , metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)    
  }
  

  if(metOptNB == "wss"){
    dcluster = as.matrix(p$data)
    d = inflexi(as.double(dcluster[,1]),as.double(dcluster[,2]),1,length(dcluster[,1]),3,3,plots=FALSE)
    nboptimal = d$finfl[1]
    print(nboptimal)
    
  }else if (metOptNB ==   "silhouette"){
    nboptimal = which(p$data[,2] == max(p$data[,2]))
  }else if (metOptNB == "gap_stat"){
    dcluster = as.matrix(p$data)
    #distorigin = abs(scale(as.double(dcluster[,5]), 0)-scale(-1*as.double(dcluster[,6]), 0))
    d = inflexi(as.double(dcluster[,5]),-1*as.double(dcluster[,6]),1,length(dcluster[,1]),3,3,plots=FALSE)
    nboptimal = d$finfl[1]
  }
  
  
  if (metcluster == "hclust"){
    outclust = hcut(din, k = nboptimal, hc_method = metagregation)
  }else if(metcluster == "kmeans"){
    outclust = hkmeans(din, nboptimal)
  }
  
  # PCA with clusters
  fviz_cluster(outclust, labelsize = 5)
  ggsave(paste(prout, "PCA_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 12, width = 12)
  
  
  # dendogram fviz
  if(metcluster == "hclust"){
    fviz_dend(outclust, show_labels = FALSE, type = "circular")
    ggsave(paste(prout, "dendov1_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 12, width = 13)
  }
  
  # dendogram old
  dcluster = cbind(rownames(din),outclust$cluster)
  colnames(dcluster) = c("names", "cluster")
  dcluster = as.data.frame(dcluster)
  
  # save cluster
  write.csv(dcluster, paste(prout, "Table_", metcluster, "_",  metagregation, "_", metOptNB, ".csv", sep = ""), row.names = FALSE)

  d <- dist(din, method = "euc")
  tupgma2 <- upgma(d, method = metagregation)
  tupgma2 = groupOTU(tupgma2, nboptimal)
  t4 <- ggtree(tupgma2, layout="circular", size=1, aes(color=cluster))
  t4 <- t4 %<+% dcluster + geom_text(aes(color=cluster, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #scale_color_continuous(low='red', high='lightgreen') +
    #scale_color_manual(values=c("grey","red")) +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  print(t4)
  ggsave(paste(prout, "dendo_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)
  
}



# clustering purity 
optimunNbClusterClassPurity = function(din, daff, aggMeth, distMeth){
  
  nmax = dim(din)[1] -10
  i = 1
  
  #i = 12
  #nmax = 12
  
  dout = data.frame()
  
  while(i <= nmax){
    outclust = hcut(din, k = i, hc_method = aggMeth, hc_metric = distMeth)
    dclust = cbind(outclust$cluster, daff[names(outclust$cluster), "classAff"])
    rownames(dclust) = names(outclust$cluster)
    
    ARI = adjustedRandIndex(dclust[,2], dclust[,1])
    lpurity = computePurityCluster(dclust)
    RI = lpurity[2]
    purity = lpurity[1]
    dout[i,1] = i
    dout[i,2] = ARI
    dout[i,3] = RI
    dout[i,4] = purity
    
    #print ("=====")
    #print(RI)
    #print(purity)
    #print(ARI)
    
    #print (lpurity)
    i = i + 1
  }
  colnames(dout) = c("NB", "ARI", "RI", "purity")
  return(dout)
}



computePurityCluster = function(dclusters){
  
  nmax = max(dclusters[,1])
  lpurity = c()
  lRI = c()
  
  for(i in seq(1,nmax)){
    dclust = dclusters[which(dclusters[,1] == i),]
    if (!is.null(dim(dclust))){
      lpurity = append(lpurity, computePu(dclust))
      lRI = append(lRI, computeRI(dclust))
    }
  }
  avPur = mean(lpurity)
  avRI = mean(lRI)
  return(c(avPur, avRI))
}



computePu = function(dcluster){
  
  dcluster = na.omit(dcluster)
  tconfus = table(dcluster[,2])
  
  pur = max(tconfus)/sum(tconfus)
  return(pur)
  
}


computeRI = function(dcluster){
  
  dcluster = na.omit(dcluster)
  tconfus = table(dcluster[,2])
  irep = (which(t(tconfus) == max(tconfus))[1])-1
  lreal = rep(1, dim(dcluster)[1])
  
  criteria = perftable(dcluster[,2], lreal)
  RI = (criteria[1] + criteria[2])/(criteria[1] + criteria[2] + criteria[3] + criteria[4])
  
  return(RI)
}




perftable = function (list_predict, list_real, verbose = 0){
  
  nb_value = length (list_real)
  i = 1
  tp = 0
  fp = 0
  tn = 0
  fn = 0
  while(i <= nb_value ){
    if (list_predict[i]==1){
      if (list_predict[i] == list_real[i]){
        tp = tp + 1
      }else {
        fp = fp + 1
      }
    }else{
      if (list_predict[i] == list_real[i]){
        tn = tn + 1
      }else {
        fn = fn + 1
      }
    }
    i = i + 1
  }
  if(verbose == 1){
    print (paste ("TP : ", tp, sep = ""))
    print (paste ("TN : ", tn, sep = ""))
    print (paste ("FP : ", fp, sep = ""))
    print (paste ("FN : ", fn, sep = ""))
  }
  tableval = c(tp,tn,fp,fn)
  return (tableval)
}
