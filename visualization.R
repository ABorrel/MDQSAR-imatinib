#!/usr/bin/env Rscript
source("cardMatrix.R")
source("PCAplot.R")
source("dendocircular.R")
source("clustering.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity
prout = args[3]
valcor = as.double(args[4])
maxQuantile = as.double(args[5])
plotPCA = args[6]
corMatrix = args[7]
histplot = args[8]
circularDendo = args[9]
optimal_clustering = args[10]

pdesc = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/1D2D.csv"
pdata = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
prout = "/home/borrela2/imatinib/results/analysis/DescLig2D3D/dendoAff_1D2D/"
valcor = 0.9
maxQuantile = 95
plotPCA = 0
corMatrix = 0
histplot = 0
circularDendo = 1
optimal_clustering = 0


#pdesc = "/home/aborrel/imitanib/results/analysis/desc/tableDescall.desc"
#pdata = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt"
#prout = "/home/aborrel/imitanib/results/analysis/desc/"
#plotPCA = 1
#corMatrix = 1
#histplot = 1
#circularDendo = 1
#optimal_clustering = 1
#valcor = 0.90

# Opening
daff = read.csv(pdata, sep = "\t", header = TRUE)
rownames(daff) = daff[,1]

orderaff = rownames(daff)[order(daff[,2],decreasing=F)]

daff = as.data.frame(daff[orderaff,])

# desc
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]
print("==Matrix after selection==")
rownames(dglobal) = dglobal[,1]
dglobal = as.data.frame(dglobal)

# order with affinity
dglobal = dglobal[orderaff,]
dglobal = na.omit(dglobal) # remove empty line
dglobal = dglobal[,-c(1,2)]
dglobal =  delnohomogeniousdistribution(dglobal, maxQuantile)
print("====last filtering====")
print(dim(dglobal))
rname = rownames(dglobal)
dglobal <- mapply(dglobal, FUN=as.numeric)
rownames(dglobal) = rname



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
      lmetagregation = c("ward.D2", "complete", "single", "average")
    }
    for (metagregation in lmetagregation){
      for(metoptimal in lmetoptimal){
        optimalCluters(dglobal, prout, metclustering, metoptimal, metagregation)
      }
    }
  }
}

if (circularDendo == 1){
  dendogramCircleAffinity(dglobal, daff, paste(prout, "dendo_", valcor, "-", maxQuantile, sep = ""))
}
