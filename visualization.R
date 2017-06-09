#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
source("PCAplot.R")
source("dendocircular.R")
source("clustering.R")


################
#     MAIN     #
################

#args <- commandArgs(TRUE)
#pdesc = args[1]
#pdata = args[2] #to take affinity
#prout = args[3]
#valcor = as.double(args[4])
#plotPCA = args[5]
#corMatrix = args[6]
#histplot = args[7]
#circularDendo = args[8]

pdesc = "/home/aborrel/imitanib/results/analysis/desc/tableDesc.csv"
pdata = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"
prout = "/home/aborrel/imitanib/results/analysis/desc/"
plotPCA = 0
corMatrix = 0
histplot = 0
circularDendo = 0
optimal_clustering = 1
valcor = 0.90

# Opening
ddata = read.csv(pdata, sep = "\t", header = TRUE)
#print(dim(ddata))
daffinity = ddata[,c("CMPD_CHEMBLID", "PCHEMBL_VALUE", "ASSAY_CHEMBLID")]
#print (daffinity)
#lcolaff <- colorRampPalette(c("white", "red"))
#lcolaff = lcolaff[dim(daffinity)[1]]
orderaff = order(daffinity[,2],decreasing=F)
#print (orderaff)

daffinity = as.data.frame(daffinity[orderaff,])
rownames(daffinity) = daffinity[,1]
#daffinity = as.data.frame(daffinity[,-1])

# desc
dglobal = openData(pdesc, valcor, prout, c(1,2))
#print (dim(dglobal[[1]]))

dglobal = dglobal[[1]]
rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]
# order with affinity
dglobal = dglobal[orderaff,]
#print(rownames(dglobal))
#print(colnames(dglobal))

#d1D_data = delnohomogeniousdistribution(d1D_data, 75)

if (corMatrix == 1){
  cardMatrixCor(cor(dglobal), paste(prout, "matrixCor_", valcor, sep = ""), 6)
}

if(plotPCA == 1){
  PCAplot(dglobal, paste(prout, "global_", valcor, sep = ""))
}

if (histplot == 1){
  histDataOne(data1 = dglobal, paste(prout, "histDesc_", valcor, ".pdf", sep = ""))
}

if (optimal_clustering ==1 ){
  
  lmetclustering = c("hclust", "kmeans")
  lmetagregation = c("ward.D2", "ward.D", "complete", "single", "average")
  lmetoptimal = c("silhouette", "wss", "gap_stat")
  
  for(metclustering in lmetclustering){
    if (metclustering == "kmeans"){
      lmetagregation = c("ward.D2")
    }else{
      lmetagregation = c("ward.D2", "ward.D", "complete", "single", "average")
    }
    for (metagregation in lmetagregation){
      for(metoptimal in lmetoptimal){
        optimalCluters(dglobal, prout, metclustering, metoptimal, metagregation)
      }
    }
  }
}

if (circularDendo == 1){
  dendogramCircle(dglobal, daffinity, paste(prout, "dendo_", valcor, sep = ""))
}