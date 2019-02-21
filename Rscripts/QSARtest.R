#!/usr/bin/env Rscript
source("tool.R")
source("MachinLearning.R")
source("performance.R")
source("dataManager.R")

library(chemmodlab)
library(rpart)
################
#     MAIN     #
################


args <- commandArgs(TRUE)
ptrain = args[1]
ptest = args[2]
pcluster = args[3]
prout = args[4]
nbCV = as.integer(args[5])



# to test
ptrain = "/home/aborrel/imitanib/results/analysis/QSARs/Lig/trainSet.csv"
ptest = "/home/aborrel/imitanib/results/analysis/QSARs/Lig/testSet.csv"
pcluster = "0"
prout = "/home/aborrel/imitanib/results/analysis/QSARs/Lig/"

# cross validation 10
nbCV = 10


# model regression #
####################
modelPCRreg = 0 
modelPLSreg = 0
modelSVMreg = 1
modelRFreg = 0
modelCartreg = 0
modelNNreg = 0
modelDLreg = 0
chemmodlabreg = 0



#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Train csv: ", ptrain, sep = ""))
print (paste("Test csv: ", ptest, sep = ""))
print (paste("Folder out: ", prout, sep = ""))
print(paste("Nb of CV: ", nbCV, sep = ""))
print("")

print("=====Machine learning=====")
print("---Regression model----")
print(paste("PCR: ", modelPCRreg, sep = ""))
print(paste("PLS: ", modelPLSreg, sep = ""))
print(paste("SVM: ", modelSVMreg, sep = ""))
print(paste("CART: ", modelCartreg, sep = ""))
print(paste("RF: ", modelRFreg, sep = ""))
print(paste("NN: ", modelNNreg, sep = ""))
print(paste("DP: ", modelDLreg, sep = ""))
print(paste("Chemmodlab: ", chemmodlabreg, sep = ""))
print("")



##############################
# Process descriptors matrix #
##############################

# training set
dtrain = read.csv(ptrain, header = TRUE)
rownames(dtrain) = dtrain[,1]
dtrain = dtrain[,-1]

# test set
dtest = read.csv(ptest, header = TRUE)
rownames(dtest) = dtest[,1]
dtest = dtest[,-1]

# cluster
if (pcluster != "0"){
  dcluster = read.csv(pcluster, header = TRUE)
  namescpd = dcluster[,1]
  dcluster = dcluster[,-1]
  names(dcluster) = namescpd
}else{
  namescpd = cbind(rownames(dtrain), rownames(dtest))
  dcluster = rep(1, length(namescpd))
  names(dcluster) = namescpd
}


print("==== Dataset ====")
print(paste("Data train: dim = ", dim(dtrain)[1], dim(dtrain)[2], sep = " "))
print(paste("Data test: dim = ", dim(dtest)[1], dim(dtest)[2], sep = " "))
print("")

# sampling data for CV #
########################
lgroupCV = samplingDataNgroup(dtrain, nbCV)
controlDatasets(lgroupCV, paste(prout, "ChecksamplingCV", nbCV, sep = ""))


##### REGRESSION MODELS #########
#################################
print("**************************")
print("*****  REGRESSSION   *****")
print("**************************")

### PCR ####
############

if (modelPCRreg == 1){
  nbCp = PCRgridCV(lgroupCV, prout)
  PCRCV(lgroupCV, nbCp, dcluster, prout)
  PCRTrainTest(dtrain, dtest, dcluster, nbCp)
}

### PLS  ####
#############

if (modelPLSreg == 1){
  nbcp = PLSCV(lgroupCV, dcluster, prout)
  PLSTrainTest(dtrain, dtest, dcluster, nbcp)
}

### SVM ###
###########

if(modelSVMreg == 1){
  vgamma = 2^(-1:1)
  vcost = 2^(2:8)
  #SVMRegCV(lgroupCV, vgamma, vcost, dcluster, prout)
  SVMRegTrainTest(dtrain, dtest, vgamma, vcost, dcluster, prout)
}

######
# RF #
######

if (modelRFreg == 1){
  vntree = c(10,50,100,200,500, 1000)
  vmtry = c(1,2,3,4,5,10,15,20, 25, 30)
  
  #RFregCV(lgroupCV, 50, 5, dcluster, prout)# for test
  parameters = RFGridRegCV(vntree, vmtry, lgroupCV,  prout)
  #RFregCV(lgroupCV, parameters[[1]], parameters[[2]], dcluster, prout)
  RFreg(dtrain, dtest, parameters[[1]], parameters[[2]], dcluster, prout)
}



############
#   CART   #
############

if(modelCartreg == 1){
  CARTRegCV(lgroupCV, dcluster, prout)
  CARTreg(dtrain, dtest, dcluster, prout)
}


#############
# CHEMMOLAB #
#############

if(chemmodlabreg == 1){
  dchem = cbind(dtrain[,dim(dtrain)[2]],dtrain[,-dim(dtrain)[2]] )
  colnames(dchem)[1] = "Aff"
  
  pdf(paste(prout, "Reg_chemmolab.pdf", sep = ""))
  fit = ModelTrain(dchem, ids = FALSE)
  CombineSplits(fit, metric = "R2")
  CombineSplits(fit, metric = "rho")
  dev.off()
}



####################
#  NEURAL NETWORK  #
####################


if(modelNNreg ==1){
  vdeacay = c(0.001, 0.01, 0.1, 1, 10, 100)
  vsize = seq(1,15)
  NNRegCV(lgroupCV, dcluster, prout)
  NNReg(dtrain, dtest, dcluster, vdeacay, vsize, prout)
}




####################
#  DEEP LEARNING   #
####################


if(modelDLreg ==1){
  #DLRegCV(lgroupCV, dcluster, prout)
  DLReg(dtrain, dtest, dcluster, prout)
}








