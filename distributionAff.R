library(ggplot2)



multipleHist = function(d, prout){

  
  ggplot(d, aes(x=Aff, fill = Type)) + 
    geom_bar(stat = "bin")+
    #geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
    #               binwidth=.5,
    #               colour="black", fill="white") +
    labs(x = "p(Affinity)", y = "Frequencies") + 
    xlim (c(4, 11))+
    #  ylim(c(0, 0.65))+
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    #geom_density(alpha=.2, fill="#FF6666")
  ggsave(paste(prout, "_all.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ltype = unique(d[, "Type"])
  
  for(type in ltype){
    print (type)
    dtype = d[which(d[,"Type"]== type),]
    print (dtype)
    if(type == "IC50"){
      cold = "red"
    }else if(type == "Ki"){
      cold = "blue"
    }else{
      cold = "green"
    }
    xexp = paste("p(", type, ")", sep = "")
    ggplot(dtype, aes(x=Aff)) + 
        geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                       binwidth=.5,
                       colour="black", fill="white") +
      labs(x = xexp, y = "Frequencies") + 
      xlim (c(4, 11))+
      #  ylim(c(0, 0.65))+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
        geom_density(alpha=.2, fill=cold)
    ggsave(paste(prout, "_", type, ".png", sep = ""), dpi = 300, width = 8, height = 7)
    
  }
  
  
}


##########
#  MAIN  #
##########

args <- commandArgs(TRUE)
paff = args[1]


paff = "/home/aborrel/imitanib/results/CHEMBL/AffAllcurated"

daff = read.table(paff, sep = "\t", header = TRUE)
multipleHist(daff, paff)


