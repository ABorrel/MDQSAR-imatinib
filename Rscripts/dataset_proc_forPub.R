#!/usr/bin/env Rscript
# script develop to built supporting information and add dataset split in the formated table 
library(ggplot2)
library(plyr)


barplotTrainTest = function(din, type_aff, pfilout){
  
  mu <- ddply(din, "set", summarise, grp.mean=mean(PCHEMBL_VALUE))
  
  p = ggplot(data = din, aes(x=PCHEMBL_VALUE, fill = set))+
    geom_histogram(binwidth=0.5, position="dodge")+#, width = 0.9, position=position_dodge())+
    #geom_text(aes(label=count), vjust=1.6, color="white", size=8)+
    labs(x=type_aff, y="Count")+
    theme(legend.position = "none", axis.title.y = element_text(size=20),axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))+
    guides(fill=guide_legend(title="Set"))+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=set),
               linetype="solid")
  
  ggsave(paste(pfilout, ".png", sep = ""),width = 7, height = 6, dpi = 300)
  
  pval = t.test(din$PCHEMBL_VALUE[which(din$set == "Testing")], din$PCHEMBL_VALUE[which(din$set == "Training")])  
  print(pval)

}


##########
#  MAIN  #
##########

args <- commandArgs(TRUE)
pdataset = args[1]
ptrainIC50 = args[2]
ptestIC50 = args[3]
ptrainKi = args[4]
ptestKi = args[5]



pdataset = "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/CHEMBL_dataset/tab_filtered_IC50-Ki_manualedit.csv"
ptrainIC50 = "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/QSAR_split/trainSet_IC50.csv"
ptestIC50 = "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/QSAR_split/testSet_IC50.csv"
ptrainKi = "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/QSAR_split/trainSet_Ki.csv"
ptestKi = "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/QSAR_split/testSet_Ki.csv"



dchem = read.csv(pdataset, sep = "\t")
rownames(dchem) = dchem[,1]

dtrainIC50 = read.csv(ptrainIC50, sep = ",")
rownames(dtrainIC50) = dtrainIC50[,1]

dtestIC50 = read.csv(ptestIC50, sep = ",")
rownames(dtestIC50) = dtestIC50[,1]

dtrainKi = read.csv(ptrainKi, sep = ",")
rownames(dtrainKi) = dtrainKi[,1]

dtestKi = read.csv(ptestKi, sep = ",")
rownames(dtestKi) = dtestKi[,1]

# add in the dataset table the set
set = rep("excluded", dim(dchem)[1])
names(set) = rownames(dchem)
set[intersect(rownames(dtrainIC50), names(set))] = "Training"
set[intersect(rownames(dtrainKi), names(set))] = "Training"

set[intersect(rownames(dtestIC50), names(set))] = "Testing"
set[intersect(rownames(dtestKi), names(set))] = "Testing"

dout = cbind(dchem, set)
dout_clean = dout[-which(dout$set == "excluded"),]

# draw barplot affinity for the training and the test set
# ki
dki = dout_clean[which(dout_clean$STANDARD_TYPE == "Ki"),]
barplotTrainTest(dki, "pKi", paste(dirname(pdataset), "/barplot_ki", sep = ""))

# IC50
dIC50 = dout_clean[which(dout_clean$STANDARD_TYPE == "IC50"),]
barplotTrainTest(dIC50, "pIC50", paste(dirname(pdataset), "/barplot_IC50", sep = ""))

# all
barplotTrainTest(dout_clean, "pall", paste(dirname(pdataset), "/barplot_all", sep = ""))


write.csv(dout, "C:/Users/Aborrel/research/NCSU/imatinib-MD/results/CHEMBL_dataset/ChEMBL_dataset.csv")






