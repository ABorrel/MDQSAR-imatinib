#!/usr/bin/env Rscript

library(ggplot2)

plotCorScoreDock = function(din, vblack, prout){
  
  cor.val = cor(din[,1], din[,2])
  p.val = cor.test (din[,1], din[,2])$p.value
  #print (cor.test (d_temp[,2], d_temp[,3]))
  #print (p.val)
  
  print(vblack)
  
  ltypeAff = unique(din$typeAff)
  
  minaff = min(din$Aff)
  maxaff = max(din$Aff)
  
  for(typeAff in ltypeAff){
    dtemp = din[din$typeAff == typeAff,]
    
    qplot(Dock_score,emodel,data=dtemp) +
      geom_point(size=4)+
      geom_point(size=3, aes(colour = Aff))+
      xlim(-20,0)+
      ylim(-180,0)+
      xlab("Docking Score")+
      ylab("Emodel Score")+
      scale_colour_gradient(paste("p", typeAff, sep = ""), low="red", high="green")
    ggsave(paste(prout, typeAff, "_Aff.png", sep=""), dpi=300, width = 8,height = 7)
    
    dtemp[which(dtemp$Aff == max(dtemp$Aff)),] = maxaff
    dtemp[which(dtemp$Aff == min(dtemp$Aff)),] = minaff
    
    qplot(Dock_score,emodel,data=dtemp) +
      geom_point(size=4)+
      geom_point(size=3, aes(colour = Aff))+
      xlim(-20,0)+
      ylim(-180,0)+
      xlab("Docking Score")+
      ylab("Emodel Score")+
      theme(legend.position = "none")+
      theme(text = element_text(size=19))+
      scale_colour_gradient(paste("p", typeAff, sep = ""), low="red", high="green")
    ggsave(paste(prout, typeAff, "_scale.png", sep=""), dpi=300, width = 8,height = 7)
    
  }
  
  
  qplot(Dock_score,emodel,data=din) +
    geom_point(size=4)+
    geom_point(size=3, aes(colour = Aff))+
    xlim(-20,0)+
    ylim(-180,0)+
    xlab("Docking Score")+
    ylab("Emodel Score")+
    theme(text = element_text(size=19))+
    scale_colour_gradient("pAffinity", low="red", high="green")
  ggsave(paste(prout, "all_scale.png", sep=""), dpi=300, width = 8,height = 7)
  
  
  # for extra interpretation 
  
  # text with legend  
  ggplot(din, aes(Dock_score,emodel, label = rownames(din))) +
    geom_text(size = 2, aes(colour = Aff)) +
    scale_colour_gradient(low="red", high="green")
  ggsave(paste(prout, "AllAssay_ggplot.png",sep=""),width = 8,height = 7)
  
  ggplot(din, aes(Dock_score,emodel, label = din$Aff)) +
    geom_text(size = 2, aes(colour = Aff)) +
    scale_colour_gradient(low="red", high="green")
  ggsave(paste(prout, "AllIC50_ggplot.png",sep=""),width = 8,height = 7)
  
  ASSAY_CHEMBLID = din$ASSAY_CHEMBLID
  ggplot(din, aes(Dock_score,emodel)) +
    geom_point(size = 2, aes( colour = ASSAY_CHEMBLID))+
    theme(legend.position = "none")
  ggsave(paste(prout, "BA_ggplot.png",sep=""),width = 8,height = 7)
  
  ggplot(din, aes(Dock_score,emodel, label = ASSAY_CHEMBLID)) +
    geom_text(size = 2, aes( colour = ASSAY_CHEMBLID))+
    theme(legend.position = "none")
  ggsave(paste(prout, "BAtext_ggplot.png",sep=""),width = 8,height = 7)
}




############
##  MAIN  ##
############

args <- commandArgs(TRUE)
p_filin = args[1]
pchembl = args[2]
prout = args[3]

#p_filin = "./../../research/imatinib/ScoreVSAff.txt"
#pchembl = "./../../research/imatinib/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"
#prout = "./../../research/imatinib/"

din = read.csv (p_filin, header = TRUE, sep = "\t")
rownames(din) = din[,1]
din = din[,-1]

# extract bioassay
dchembl = read.csv(pchembl, header = TRUE, sep = "\t")
rownames(dchembl) = dchembl$CMPD_CHEMBLID

ASSAY_CHEMBLID = dchembl[rownames(din), c("ASSAY_CHEMBLID")]

din = cbind(din, ASSAY_CHEMBLID)
# glevec
plotCorScoreDock(din, c("CHEMBL941"), prout)



#biofilm_3IX4_XP=read.csv(p_filin,header=T,stringsAsFactors = F, sep = "\t")
#pIC50= biofilm_3IX4_XP$Aff
#qplot(Dock_score,emodel,data=biofilm_3IX4_XP) +
#  geom_point(colour="black",size=4)+
#  geom_point(aes(colour = pIC50))+
#  scale_colour_gradient(low="red", high="green")
#ggsave(paste("file_of_your_image.png",sep=""),width = 8,height = 7)

