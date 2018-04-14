


library(ggplot2)
source("performance.R")



plotRFperf = function(dpred, prout){
  
  valr2 = calR2(dpred[,2], dpred[,1])
  colnames(dpred) = c("Yreal", "Ypredict")
  
  p = ggplot(dpred, aes(Yreal, Ypredict))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=5, y=11.5, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    xlim (c(4, 12)) +
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  ggsave(paste(prout, "PerfRFregpoint_CV10.png", sep = ""), width = 6,height = 6, dpi = 300)
  
  
  #dpred = cbind(dtest[,"Aff"], vpredtest)
  #Vcluster = dcluster
  #d = cbind(dpred, Vcluster)
  #colnames(dpred) = c("NAME", "Yreal", "Ypredict", "cluster")
  #dpred = as.data.frame(dpred)
  #print(dpred)
  
  p = ggplot(dpred, aes(Yreal, Ypredict, label=rownames(dpred)))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=5, y=11.5, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    geom_text(size = 2.6, aes(label= rownames(dpred)), parse = TRUE, color="black", nudge_y = 0.06) + 
    labs(x = "pAff", y = "Predict pAff") + 
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    xlim (c(4, 12)) +
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  print(p)
  ggsave(paste(prout, "PerfRFregname_CV10.png", sep = ""), width = 8,height = 8, dpi = 300)
}

  
# performance

ppred = "./../results/analysis/QSARs/Lig-Lig2D_Ki/perfRFRegCV_10.txt"
prout = "./../results/analysis/QSARs/Lig-Lig2D_Ki/"


# test - plot for publication
dpred = read.table(ppred, header = TRUE)
dpred = as.data.frame(dpred)


plotRFperf(dpred, prout)
