#!/usr/bin/env Rscript
source ("tool.R")
source("dataManager.R")
source("dendocircular.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity or class
prout = args[3]
valcor = args[4]
maxquantile = as.double(args[5])
logaff = as.integer(args[6])

#pdesc = "/home/borrela2/imatinib/results/analysis/QSARs/Lig/descGlobal"
#pdata = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
#prout = "/home/borrela2/imatinib/results/analysis/QSARs/Lig/VisuData/"
#valcor = 0.80
#logaff = 0
#maxquantile = 80



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



##############################
# drawn circular dendrogram  #
##############################

dendogramCircleAff(dglobal, daffinity, paste(prout, "tree"))




