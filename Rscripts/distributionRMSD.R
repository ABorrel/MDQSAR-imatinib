#!/usr/bin/env Rscript
library(ggplot2)



multipleHist = function(d, prout){

  
  ggplot(d, aes(x=RMSDlig, fill = TypeAff)) + 
    geom_bar(stat = "bin")+
    #geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
    #               binwidth=.5,
    #               colour="black", fill="white") +
    labs(x = "RMSD ligand (Ang)", y = "Frequencies") + 
    #xlim (c(4, 11))+
    #  ylim(c(0, 0.65))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    #geom_density(alpha=.2, fill="#FF6666")
  ggsave(paste(prout, "RMSDlig_all.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ggplot(d, aes(x=RMSDca, fill = TypeAff)) + 
    geom_bar(stat = "bin")+
    #geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
    #               binwidth=.5,
    #               colour="black", fill="white") +
    labs(x = "RMSD-Ca (Ang)", y = "Frequencies") + 
    #xlim (c(4, 11))+
    #  ylim(c(0, 0.65))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    #geom_density(alpha=.2, fill="#FF6666")
    ggsave(paste(prout, "RMSDca_all.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ggplot(d, aes(x=RMSDall, fill = TypeAff)) + 
    geom_bar(stat = "bin")+
    #geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
    #               binwidth=.5,
    #               colour="black", fill="white") +
    labs(x = "RMSDall (Ang)", y = "Frequencies") + 
    #xlim (c(4, 11))+
    #  ylim(c(0, 0.65))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    #geom_density(alpha=.2, fill="#FF6666")
    ggsave(paste(prout, "RMSDall_all.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  
  ltype = unique(d[, "TypeAff"])
  
  for(type in ltype){
    print (type)
    dtype = d[which(d[,"TypeAff"]== type),]
    print (dtype)
    if(type == "IC50"){
      cold = "red"
    }else if(type == "Ki"){
      cold = "blue"
    }else{
      cold = "green"
    }
    xexp = "RMSD (Ang)"
    
    
    
    ggplot(dtype, aes(x=RMSDlig)) + 
        geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                       binwidth=.5,
                       colour="black", fill="white") +
      labs(x = xexp, y = "Frequencies") + 
      #xlim (c(1, 3))+
      #  ylim(c(0, 0.65))+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
        geom_density(alpha=.2, fill=cold)
    ggsave(paste(prout, "RMSDlig_", type, ".png", sep = ""), dpi = 300, width = 8, height = 7)
   
    
    
    ggplot(dtype, aes(x=RMSDca)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.1,
                     colour="black", fill="white") +
      labs(x = xexp, y = "Frequencies") + 
      #xlim (c(1, 3))+
      #  ylim(c(0, 0.65))+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill=cold)
    ggsave(paste(prout, "RMSDca_", type, ".png", sep = ""), dpi = 300, width = 8, height = 7)
    
    
    
    ggplot(dtype, aes(x=RMSDall)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.1,
                     colour="black", fill="white") +
      labs(x = xexp, y = "Frequencies") + 
      #xlim (c(1, 3))+
      #  ylim(c(0, 0.65))+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill=cold)
    ggsave(paste(prout, "RMSDall_", type, ".png", sep = ""), dpi = 300, width = 8, height = 7)
    
     
  }
  
  
}


##########
#  MAIN  #
##########

args <- commandArgs(TRUE)
paff = args[1]
prout = args[2]

#paff = "/home/aborrel/imitanib/results/CHEMBL/AffAllcurated"

daff = read.table(paff, sep = "\t", header = TRUE)
# remove kd
daff = daff[-which(daff$TypeAff == "Kd"),]

multipleHist(daff, prout)


