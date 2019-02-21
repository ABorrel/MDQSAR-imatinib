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
  
#pdesc = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_All/descGlobal_clean.csv"
#paff = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
#cutoff = 7.5
#aggMeth = "ward.D2"
#distMeth = "euclidean"
#clustMeth = "hclust"
#prout = "/home/borrela2/imatinib/results/analysis/ClusteringDescType/Lig2D_All/"
  


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

# dendo 
dendogramCircleClassAff(ddesc, daff, prout)
dendogramCircleAffinity(ddesc, daff, prout)
dendogramCircle(ddesc, daff, paste(prout, "dendo.png", sep = ""))

# study purity of number of cluster
dopt = optimunNbClusterClassPurity(ddesc, daff, aggMeth, distMeth)
doptimal = dopt[which(dopt[,2] == max(dopt[,2])), ]


write.csv(dopt, file = paste(prout, "purityClust.csv"))
write.csv(doptimal, file = paste(prout, "optimalClust.csv"))



