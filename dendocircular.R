#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)


dendogramCircleAffinity = function (ddesc, daff, pfilout){
  
    
  daff = as.data.frame(daff)
  matTrans1 <- scale(ddesc)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes())
  t1 <- t1 %<+% daff + geom_tippoint(aes(color=Aff), alpha=0.65, size=2)+
    scale_color_continuous(low='red', high='lightgreen')+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  dtypeplot =  as.data.frame(daff[,3])
  rownames(dtypeplot) = rownames(daff)
  
  daffplot = cbind(daff[,2], daff[,2])
  rownames(daffplot) = rownames(daff)
  daffplot = as.data.frame(daffplot[,-1])
  
  print(dtypeplot)
  t2 = gheatmap(t1, dtypeplot, font.size =0, offset = 9, width = 0.05, colnames_offset_x = 2, colnames_offset_y = 0)
  
  #t3 = gheatmap(t2, daffplot, font.size = 2, offset = 2, width = 0.05, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")

  open_tree(t2, 15) %>% rotate_tree(15)

    ggsave(paste(pfilout, "typeAff.png", sep = ""), dpi=300, height = 11, width = 11)
  
  
    
    t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes())
    t1 <- t1 %<+% daff + geom_tippoint(aes(color=Aff), alpha=0.65, size=2)+
      scale_color_continuous(low='red', high='lightgreen')+
      theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
      geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                     color = "white", linesize = 1E-100, fontsize = 1E-100)
    
    dtypeplot =  as.data.frame(daff[,3])
    rownames(dtypeplot) = rownames(daff)
    
    daffplot = as.data.frame(cbind(daff[,2], rep(0, length(daff[,2]))))
    daffplot = as.data.frame(cbind(daffplot, daff[,2]))
    rownames(daffplot) = rownames(daff)
    
    print(dtypeplot)
    t3 = gheatmap(t1, daffplot, font.size = 0, offset = 2, width = 0.15, colnames_offset_x = 2, colnames_offset_y = 0, low = "red", high = "lightgreen")
    
    open_tree(t3, 15) %>% rotate_tree(15)
    
    ggsave(paste(pfilout, "Aff.png", sep = ""), dpi=300, height = 11, width = 11)
    
    
    
  
  #pfilout = paste(prout, "dendo_cluster.png", sep = "")
  
  #t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes(color=cluster))
  #t1 <- t1 %<+% dcluster +
  #  geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
  #theme(legend.position="right")+
  #  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  #  geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
  #                 color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  #daff = daff[,-1]
  #t2 = gheatmap(t1, daff, font.size =0 , offset = 3, width = 0.2, colnames_offset_x = 2, colnames_offset_y = -0.5, low = "red", high = "lightgreen") +
  #scale_color_continuous(low='red', high='lightgreen') +
  #  theme_tree()
  #print (t2)
  #open_tree(t2, 15) %>% rotate_tree(15)
  #ggsave(pfilout, dpi=300, height = 11, width = 11)
  
  
  
  
  
  
  
}



dendogramCircle = function(ddes, daff, pfilout){
  
  print (daff)
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=PCHEMBL_VALUE, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=2) +
    geom_tippoint(aes(color=PCHEMBL_VALUE), alpha=0.75, size=1)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
                 color = "white", linesize = 1E-100, fontsize = 1E-100)
  print(t4)
  ggsave(paste(pfilout, "aff.png", sep = ""), dpi=300, height = 8, width = 9)
  
  
  #t4 <- ggtree(tupgma2, layout="circular", size=1)
  #t4 <- t4 %<+% daff + geom_text(aes(color=ASSAY_CHEMBLID, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=2) +
  #  geom_tippoint(aes(color=ASSAY_CHEMBLID), alpha=0.75, size=1)+
    #scale_color_continuous(low='red', high='lightgreen') +
    #theme(legend.position="right")+
  #  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  #  geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
  #                 color = "white", linesize = 1E-100, fontsize = 1E-100)
  #print(t4)
  #ggsave(paste(pfilout, "BA.png", sep = ""), dpi=300, height = 8, width = 9)
  
  
}

#metabol=read.csv("your_file_with_fingerprint.csv",header=T,stringsAsFactors = F)


#dockingscore=read.csv("your_file_with_dockingscore.txt",header = T,stringsAsFactors = F)
#rownames(dockingscore)=dockingscore$Metabolite
#rownames(metabol)=metabol$Molecule.name
#dockingscore=dockingscore[metabol$Molecule.name,]



  
#dockingscoredata <- cbind(as.data.frame(row.names(metabol)),as.data.frame(dockingscore))




