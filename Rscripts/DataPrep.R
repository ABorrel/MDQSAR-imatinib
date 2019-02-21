#!/usr/bin/env Rscript
source ("tool.R")
source("dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
valcor = args[2]
maxquantile = as.double(args[3])
prout = args[4]

#pdesc = "/home/borrela2/imatinib/results/analysis/QSARs/Lig-FPI-BS/descGlobal"
#valcor = 0.80
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

write.csv(dglobal, paste(pdesc, "_clean.csv", sep = ""))






