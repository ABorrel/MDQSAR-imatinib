#!/usr/bin/env Rscript
source ("tool.R")
source("dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity or class
prout = args[3]
valcor = args[4]
maxquantile = as.double(args[5])
proptraintest = args[6]
logaff = args[7]

#pdesc = "/home/aborrel/imitanib/results/analysis/QSARs/Lig/descGlobal"
#pdata = "/home/aborrel/imitanib/results/CHEMBL/AffAllcurated"
#prout = "/home/aborrel/imitanib/results/analysis/QSARs/Lig/"
#valcor = 0.80
#logaff = 0
#maxquantile = 80
#proptraintest = 0.15



##############################
# Process descriptors matrix #
##############################

dglobal = openData(pdesc, valcor, prout)
dglobal = dglobal[[1]]

print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

##########
# filter #
##########

dglobal = delnohomogeniousdistribution(dglobal, maxquantile)
print(paste("Data after filtering: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = "\t", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]

# transform #
if(logaff == 1){
  daffinity = -log10(daffinity)
}

# merge with data descriptors and remove data remove from the manual curation
lID = intersect(rownames(daffinity), rownames(dglobal))
dglobal = dglobal[lID,]
daffinity = daffinity[lID,]

##################
# divide dataset #
##################
ltraintest = samplingDataFraction(dglobal, proptraintest)
dtrain = ltraintest[[1]]
dtest = ltraintest[[2]]


#################
# Add affinity  #
#################

# training set
Aff = daffinity[rownames(dtrain),1]
dtrainglobal = cbind(dtrain, Aff)
write.csv(dtrainglobal, paste(prout, "trainSet.csv", sep = ""))

# test set
Aff = daffinity[rownames(dtest),1]
dtestglobal = cbind(dtest, Aff)
write.csv(dtestglobal, paste(prout, "testSet.csv", sep = ""))
  
lcontrol = list(dtrainglobal, dtestglobal)
controlDatasets(lcontrol, paste(prout, "qualitySplit", sep = ""))







