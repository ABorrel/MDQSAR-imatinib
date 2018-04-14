#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)



dendogramCircleAff = function(ddes, daff, pfilout){
  
  # add column with name
  daff = cbind(rownames(daff), daff)
  colnames(daff)[2] = "pAff"
  print(dim(ddes))
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=pAff, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=2)+
    geom_tippoint(aes(color=pAff), alpha=0.75, size=1)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) 
    geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  print(t4)
  ggsave(paste(pfilout, ".png", sep = ""), dpi=300, height = 8, width = 9)

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




