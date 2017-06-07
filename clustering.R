#!/usr/bin/env Rscript
library(factoextra)
library(ggplot2)
require(cluster)
library(RootsExtremaInflections)


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




#tree <- groupOTU(tupgma, cls)
#t4 <- ggtree(tree, layout="circular", size=0.5, branch.length="none",aes(color=group)) +
#  geom_text(aes(color=group, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=3)+
#  scale_color_manual(values=c("grey",color18)) + #theme(legend.position="none")+
#  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
#  geom_treescale(x = 30, y = 30, width = NULL, offset = NULL,
#                 color = "white", linesize = 1E-100, fontsize = 1E-100)
#print(t4)