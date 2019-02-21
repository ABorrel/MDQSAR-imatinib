#!/usr/bin/env Rscript
source("dendocircular.R")
source("clustering.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
paff = args[2]
typeaff = args[3]
cutoff = as.double(args[4]) #to take affinity
aggMeth = args[5]
distMeth = args[6]
clustMeth = args[7]
prout = args[8]


#pdesc = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/descGlobal_clean.csv"
#paff = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/top50/Aff50.csv"
#cutoff = 7.5
#aggMeth = "ward.D2"
#distMeth = "euclidean"
#clustMeth = "hclust"
#prout = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/top50/"
#typeaff = "All"

#pdesc = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/descGlobal_clean.csv"
#paff = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
#typeaff = "Ki"
#prout = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/"
  
#pdesc = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/descGlobal_clean.csv"
#paff = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
#cutoff = 7.5
#aggMeth = "ward.D2"
#distMeth = "euclidean"
#clustMeth = "hclust"
#prout = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_Ki/"
#ypeaff = "Ki"  


# open file #
#############
daff = read.csv(paff, sep = "\t", header = TRUE)
rownames(daff) = daff[,1]
daff = as.data.frame(daff)

if(typeaff != "All"){
  daff = daff[which(daff$Type == typeaff),]
}

# define active or inactive
classAff = rep("Active",dim(daff)[1]) # inverse for dendogram
classAff[which(daff$Aff>cutoff)] = "Inactive"
daff = cbind(daff, classAff)


# desc #
########
ddesc = read.csv(pdesc, sep = ",")
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]
ddesc = ddesc[intersect(rownames(daff), rownames(ddesc)),]


# clustering #
##############

daff = as.data.frame(daff)
matTrans1 <- scale(ddesc)
d <- dist(matTrans1, method = "euc")

#cutoff distance
cE = cutoff

d = as.matrix(d)

i = 1
imax = dim(d)[1]
cliff = NULL
while (i < imax){
  j = i + 1
  while(j <=imax){
    if(d[i,j] <= cE){
      diffAff = abs(daff[rownames(d)[i],2] - daff[rownames(d)[j],2])
      if(diffAff >= 1.5){
        cliff = rbind(cliff, c(rownames(d)[i], rownames(d)[j], d[i,j], diffAff))
      }
    }
    j = j + 1
  }
  i = i + 1
}

print (cliff)

colnames(cliff) = c("ID1", "ID2", "Dist", "diffAff")
write.csv(cliff, paste(prout, "activityCliff.csv", sep = ""))

