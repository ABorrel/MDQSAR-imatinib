#!/usr/bin/env Rscript

source("performance.R")
library("pls")
source("tool.R")
source("dataManager.R")
library (randomForest)
library (MASS)
library(rpart)
library(rpart.plot)
library(e1071)
library(ggplot2)
#library(neuralnet)
library(nnet)
library(clusterGeneration)
library(stringr)
library(reshape2)

##################
#    PCR MODELS  #
##################

# PCR - cross validation for grid
PCRgridCV = function(lfolds, prout){
  
  # performance in CV and number of cmp
  # return best number of cmp
  print(paste("== PCR in CV with ", length(lfolds), " Automatic optimization by folds ===", sep = ""))

  
  maxCp = dim(lfolds[[1]])[2] - 1
  
  vcpRMSEP = NULL
  vcpR2 = NULL
  vcpcor = NULL
  for (cp in seq(1,maxCp)){
    i = 1
    imax = length(lfolds)
    vpred = NULL
    vreal = NULL
    while(i <= imax){
      dtrain = NULL
      dtest = NULL
      for (j in seq(1:imax)){
        if (j == i){
          dtest = lfolds[[j]]
        }else{
          dtrain = rbind(dtrain, lfolds[[j]])
        }
      }
      
      modelpcr = pcr(Aff~., data=dtrain, ncomp = cp)
      predpcr = predict(modelpcr, ncomp = cp, dtest, type = "response")
      
      vpred = append(vpred, predpcr)
      vreal = append(vreal, dtest[,"Aff"])
      i = i + 1
    }
    
    #print(vpred)
    #print(vreal)
    
    corpred = cor(vreal, vpred)
    rmsepcp = vrmsep(vreal, vpred)
    R2cp = calR2(vreal, vpred)
    
    vcpcor = append(vcpcor, corpred)
    vcpRMSEP = append(vcpRMSEP,rmsepcp)
    vcpR2 = append(vcpR2, R2cp)
  }
  
  pdf(paste(prout ,"gridPCRCV", length(lfolds), ".pdf", sep = ""))
  plot(seq(1,length(vcpR2)), vcpR2, type = "l", cex = 3, main = paste("R2 by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "R2")
  plot(seq(1,length(vcpRMSEP)), vcpRMSEP, type = "l", cex = 3, main = paste("RMSEP by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "RMSEP")
  plot(seq(1,length(vcpcor)), vcpcor, type = "l", cex = 3, main = paste("Cor val by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "Cor")
  dev.off()
  
  # optimal number of component
  outcp = cbind(seq(1, length(vcpcor)), vcpcor)
  ## optimization using RMSEQ
  #dout = NULL
  #for (i in seq(1,length(vcpcor))){
  #  dout = append(dout, eucdist(0, outcp[i,1], 0, outcp[i,2]))
  #}
  #outcp = cbind(outcp, dout)
  nbCPoptimun = outcp[which(outcp[,2] == max(outcp[,2])),1]
  
  
  print(paste("Optimal component: ", nbCPoptimun, sep = ""))
  print("Perfomances Grid in CV")
  print(paste("R2=", vcpR2[nbCPoptimun], sep = ""))
  print(paste("Cor=", vcpcor[nbCPoptimun], sep = ""))
  print(paste("RMSEP=", vcpRMSEP[nbCPoptimun], sep = ""))
  print("")
  print("")
  return(nbCPoptimun)
}

# PCR - CV
PCRCV = function(lfolds, nbcomp, dcluster, prout){
  
  # performance in CV and number of cmp
  # return best number of cmp
  print(paste("==== PCR in CV with ", length(lfolds), "folds Nb compound:", nbcomp, "====",sep = ""))
  
  imax = length(lfolds)
  vpred = NULL
  vreal = NULL
  i = 1
  while(i <= imax){
    dtrain = NULL
    dtest = NULL
    for (j in seq(1:imax)){
      if (j == i){
        dtest = lfolds[[j]]
      }else{
        dtrain = rbind(dtrain, lfolds[[j]])
      }
    }
    modelpcr = pcr(Aff~., data=dtrain, ncomp = nbcomp)
    predpcr = predict(modelpcr, ncomp = nbcomp, dtest, type = "response")
    
    #print(summary(modelpcr)$coefficient)
    names(predpcr) = rownames(dtest)
    vpred = append(vpred, predpcr)
    vreal = append(vreal, dtest[,"Aff"])
    i = i + 1
  }
  corpred = cor(vreal, vpred)
  rmsepcp = vrmsep(vreal, vpred)
  R2cp = calR2(vreal, vpred)
  MAEcp = MAE(vreal, vpred)
  R02cp = R02(vreal, vpred)
  
  #print (paste(prout, "PerfPCRreg_CV", length(lfolds), ".png", sep = ""))
  
  pdf(paste(prout, "PerfPCRreg_CV", length(lfolds), ".pdf", sep = ""), 20, 20)
  plot(vreal, vpred, type = "n", cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 5)
  text(vreal, vpred, labels = names(vpred))
  abline(a = 0, b = 1, col = "red", lty = 4)
  
  plot(vreal, vpred, pch = 19, cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", lty = 4)
  
  # for cluster 
  dcluster = dcluster[names(vpred)]
  plot(vreal, vpred, type = "n", cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 5)
  text(vreal, vpred, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", lty = 4)
  dev.off()  
  
  tperf = cbind(vpred, vreal)
  
  #write.table(tperf, paste(prout, "perfPCRRegCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  print("Perfomances in CV")
  print(paste("R2=", R2cp, sep = ""))
  print(paste("R02=", R02cp, sep = ""))
  print(paste("MAE=", MAEcp, sep = ""))
  print(paste("Cor=", corpred, sep = ""))
  print(paste("RMSEP=", rmsepcp, sep = ""))
  print("")
  print("")
}
  
# PCR - real #
PCRTrainTest = function(dtrain, dtest, dcluster, nbcp){
  
  modelpcr = pcr(Aff~., data=dtrain, ncomp = nbcp)
  predpcrtest = predict(modelpcr, ncomp = nbcp, newdata = dtest)
  predpcrtrain = predict(modelpcr, ncomp = nbcp, newdata = dtrain)
  
  names(predpcrtrain) = rownames(dtrain)
  names(predpcrtest) = rownames(dtest)
  
  cortrain = cor(dtrain[,"Aff"], predpcrtrain)
  cortest = cor(dtest[,"Aff"], predpcrtest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predpcrtrain)
  rmseptest = vrmsep(dtest[,"Aff"], predpcrtest)
  
  R2train = calR2(dtrain[,"Aff"], predpcrtrain)
  R2test = calR2(dtest[,"Aff"], predpcrtest)
  
  R02train = R02(dtrain[,"Aff"], predpcrtrain)
  R02test = R02(dtest[,"Aff"], predpcrtest)
  
  MAEtrain = MAE(dtrain[,"Aff"], predpcrtrain)
  MAEtest = MAE(dtest[,"Aff"], predpcrtest)
  
  pdf(paste(prout, "PerfPCRreg_TrainTest.pdf", sep = ""), 20, 20)

  
  # print performances
  print("====PCR model on external test====")
  print(paste("NB components = ", nbcp, sep = ""))
  print(paste("Perf training (dim = ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("R2=", R2train))
  print(paste("R02=", R02train))
  print(paste("MAE=", MAEtrain))
  print(paste("Corval=", cortrain))
  print(paste("RMSEP=", rmseptrain))
  print("****")
  print(paste("Perf test (dim = ", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("R2=", R2test))
  print(paste("R02=", R02test))
  print(paste("MAE=", MAEtest))
  print(paste("Corval=", cortest))
  print(paste("RMSEP=", rmseptest))
  print("")
  print("")
  
  
  # train
  plot(dtrain[,"Aff"], predpcrtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predpcrtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], predpcrtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predpcrtrain), digits = 3)))
  text(dtrain[,"Aff"], predpcrtrain, labels = names(predpcrtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(predpcrtrain)]
  plot(dtrain[,"Aff"], predpcrtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predpcrtrain), digits = 3)))
  text(dtrain[,"Aff"], predpcrtrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # test
  plot(dtest[,"Aff"], predpcrtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], predpcrtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], predpcrtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predpcrtest), digits = 3)))
  text(dtest[,"Aff"], predpcrtest, labels = names(predpcrtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclusterTest = dcluster[names(predpcrtest)]
  plot(dtest[,"Aff"], predpcrtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predpcrtest), digits = 3)))
  text(dtest[,"Aff"], predpcrtest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()
  
  return(modelpcr)
  
}


##################
#    PLS MODELS  #
##################


# PLS in cross validation
PLSCV = function(lfolds, dcluster, prout){
  
  # performance in CV and number of cmp
  # return best number of cmp
  
  maxCp = dim(lfolds[[1]])[2]-1
  
  vcpRMSEP = NULL
  vcpR2 = NULL
  vcpcor = NULL
  vcpR02 = NULL
  vcpMAE = NULL
  for (cp in seq(1,maxCp)){
    i = 1
    imax = length(lfolds)
    vpred = NULL
    vreal = NULL
    while(i <= imax){
      dtrain = NULL
      dtest = NULL
      for (j in seq(1:imax)){
        if (j == i){
          dtest = lfolds[[j]]
        }else{
          dtrain = rbind(dtrain, lfolds[[j]])
        }
      }
      modelpls = plsr(Aff~., data=dtrain, ncomp = cp)
      predpls = predict(modelpls, ncomp = cp, dtest, type = "response")
      
      names(predpls) = rownames(dtest)
      vpred = append(vpred, predpls)
      vreal = append(vreal, dtest[,"Aff"])
      i = i + 1
    }

    corpred = cor(vreal, vpred)
    rmsepcp = vrmsep(vreal, vpred)
    valR2 = calR2(vreal, vpred)
    R02pred = R02(vreal, vpred)
    MAEpred = MAE(vreal, vpred)
    
    vcpR2 = append(vcpR2, valR2)
    vcpRMSEP = append(vcpRMSEP,rmsepcp)
    vcpcor = append(vcpcor, corpred)
    vcpR02 = append(vcpR02, R02pred)
    vcpMAE = append(vcpMAE, MAEpred)
    
  }
  
  pdf(paste(prout ,"PerfPLSRef_CV10.pdf", sep = ""), width = 20, height = 20)
  plot(seq(1,length(vcpR2)), vcpR2, type = "l", cex = 3, main = paste("R2 by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "R2")
  plot(seq(1,length(vcpRMSEP)), vcpRMSEP, type = "l", cex = 3, main = paste("RMSEP by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "RMSEP")
  plot(seq(1,length(vcpcor)), vcpcor, type = "l", cex = 3, main = paste("Cor by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "Cor")
  
  plot(vreal, vpred, type = "n", cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 5)
  text(vreal, vpred, labels = names(vpred))
  abline(a = 0, b = 1, col = "red", lty = 4)
  plot(vreal, vpred, pch = 19, cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", lty = 4)
  
  # by cluster
  dcluster = dcluster[names(vpred)]
  plot(vreal, vpred, type = "n", cex.lab = 2.5, main = paste("Correlation = ", round(cor(vreal, vpred), digits = 3)), cex = 5)
  text(vreal, vpred, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", lty = 4)
  dev.off()
  
  # optimal number of component
  outcp = cbind(seq(1, length(vcpcor)), vcpcor)

  #for (i in seq(1,length(vcpRMSEP))){
  #  dout = append(dout, eucdist(0, outcp[i,1], 0, outcp[i,2]))
  #}
  #outcp = cbind(outcp, dout)
  nbCPoptimun = outcp[which(outcp[,2] == max(outcp[,2])),1]
  print("====PLS model on CV====")
  print(paste("Optimal component: ", nbCPoptimun, sep = ""))
  print("Perfomances in CV")
  print(paste("R2=", vcpR2[nbCPoptimun], sep = ""))
  print(paste("R02=", vcpR02[nbCPoptimun], sep = ""))
  print(paste("MAE=", vcpMAE[nbCPoptimun], sep = ""))
  print(paste("Cor=", vcpcor[nbCPoptimun], sep = ""))
  print(paste("RMSEP=", vcpRMSEP[nbCPoptimun], sep = ""))
  print("")
  print("")
  
  return (nbCPoptimun)
  
}

# PLS - real #
PLSTrainTest = function(dtrain, dtest, dcluster, nbcp){
  
  modelpls = plsr(Aff~., data=dtrain, ncomp = nbcp)
  predplstest = predict(modelpls, ncomp = nbcp, newdata = dtest)
  predplstrain = predict(modelpls, ncomp = nbcp, newdata = dtrain)
  
  names(predplstrain) = rownames(dtrain)
  names(predplstest) = rownames(dtest)
  
  r2train = calR2(dtrain[,"Aff"], predplstrain)
  r2test = calR2(dtest[,"Aff"], predplstest)
  
  cortrain = cor(dtrain[,"Aff"], predplstrain)
  cortest = cor(dtest[,"Aff"], predplstest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predplstrain)
  rmseptest = vrmsep(dtest[,"Aff"], predplstest)
  
  R02train = R02(dtrain[,"Aff"], predplstrain)
  R02test = R02(dtest[,"Aff"], predplstest)
  
  MAEtrain = MAE(dtrain[,"Aff"], predplstrain)
  MAEtest = MAE(dtest[,"Aff"], predplstest)
  
  print("===== PLS model train-Test =====")
  #print(modelpls$coefficients)
  print(paste("NB components = ", nbcp, sep = ""))
  print(paste("Perf training (dim= ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("R2 train=", r2train))
  print(paste("R02 train=", R02train))
  print(paste("MAE train=", MAEtrain))
  print(paste("Cor train=", cortrain))
  print(paste("RMSEP train=", rmseptrain))
  print("")  
  print(paste("Perf test (dim=", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("R2 test=", r2test))
  print(paste("R02 test=", R02test))
  print(paste("MAE test=", MAEtest))
  print(paste("Cor test=", cortest))
  print(paste("RMSEP test=", rmseptest))
  print("")
  print("")
  
  
  # train
  pdf(paste(prout ,"PerfPLSreg_TrainTest.pdf", sep = ""), width = 20, height = 20)
  plot(dtrain[,"Aff"], predplstrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predplstrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], predplstrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predplstrain), digits = 3)))
  text(dtrain[,"Aff"], predplstrain, labels = names(predplstrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(predplstrain)]
  plot(dtrain[,"Aff"], predplstrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predplstrain), digits = 3)))
  text(dtrain[,"Aff"], predplstrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # test
  plot(dtest[,"Aff"], predplstest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], predplstest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], predplstest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predplstest), digits = 3)))
  text(dtest[,"Aff"], predplstest, labels = names(predplstest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclusterTest = dcluster[names(predplstest)]
  plot(dtest[,"Aff"], predplstest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predplstest), digits = 3)))
  text(dtest[,"Aff"], predplstest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()
  
  return(modelpls)
  
}


#################
#      SVM      #
#################

######################
# case of regression #
######################

SVMRegCV = function(lfolds, vgamma, vcost, dcluster, prout){
  
  print(paste("==== SVM in CV with ", length(lfolds), " Automatic optimization CV ====", sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  timportance = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modtune = SVMTune(dtrain, vgamma, vcost, 10)
    
    #print(dtest)
    y_real = append(y_real, dtest[,"Aff"])
    dtest = dtest[,-c(which(colnames(dtest) == "Aff"))]
    
    vpred = predict (modtune, dtest)
    names(vpred) = rownames(dtest)
    y_predict = append(y_predict, vpred)
    
    k = k + 1
  }
  
  # performances
  valr2 = calR2(y_real, y_predict)
  corval = cor(y_real, y_predict)
  RMSEP = vrmsep(y_real, y_predict)
  MAEval = MAE(y_real, y_predict)
  R02val = R02(y_real, y_predict)
  
  print("Perfomances in CV")
  print(paste("R2=", valr2, sep = ""))
  print(paste("R02=", R02val, sep = ""))
  print(paste("MAE=", MAEval, sep = ""))
  print(paste("Cor=", corval, sep = ""))
  print(paste("RMSEP=", RMSEP, sep = ""))
  print("")
  print("")
  
  
  pdf(paste(prout, "PerfSVMreg_CV", length(lfolds), ".pdf", sep = ""), 20, 20)
  plot(y_real, y_predict, pch = 20, main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster
  dcluster = dcluster[names(y_predict)]
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()  
  
  tperf = cbind(y_predict, y_real)
  write.table(tperf, paste(prout, "perfSVMRegCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  return(modtune$best.model)
}


SVMRegTrainTest = function(dtrain, dtest, vgamma, vcost, dcluster, prout){
  
  print(paste("==== SVM in train-test --- Automatic optimization CV-10====", sep = ""))
  
  # optimisation on CV-10
  modelsvm = SVMTune(dtrain, vgamma, vcost, 10)
  
  predsvmtest = predict(modelsvm, dtest[,-c(which(colnames(dtest) == "Aff"))])
  predsvmtrain = predict(modelsvm, dtrain[,-c(which(colnames(dtrain) == "Aff"))])
  
  names(predsvmtrain) = rownames(dtrain)
  names(predsvmtest) = rownames(dtest)
  
  r2train = calR2(dtrain[,"Aff"], predsvmtrain)
  r2test = calR2(dtest[,"Aff"], predsvmtest)
  
  cortrain = cor(dtrain[,"Aff"], predsvmtrain)
  cortest = cor(dtest[,"Aff"], predsvmtest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predsvmtrain)
  rmseptest = vrmsep(dtest[,"Aff"], predsvmtest)
  
  R02train = R02(dtrain[,"Aff"], predsvmtrain)
  R02test = R02(dtest[,"Aff"], predsvmtest)
  
  MAEtrain = MAE(dtrain[,"Aff"], predsvmtrain)
  MAEtest = MAE(dtest[,"Aff"], predsvmtest)
  
  print("===== SVM model train-Test =====")
  #print(modelpls$coefficients)
  print(paste("Perf training (dim= ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("R2 train=", r2train))
  print(paste("R02 train=", R02train))
  print(paste("MAE train=", MAEtrain))
  print(paste("Corval train=", cortrain))
  print(paste("RMSEP=", rmseptrain))
  print("")
  print("")
  
  
  print(paste("Perf test (dim=", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("R2 test=", r2test))
  print(paste("R02 test=", R02test))
  print(paste("MAE test=", MAEtest))
  print(paste("Corval test=", cortest, sep = ""))
  print(paste("RMSEP test=", rmseptest, sep = ""))
  print("")
  print("")
  
  
  # train
  pdf(paste(prout ,"PerfSVMreg_TrainTest.pdf", sep = ""), width = 20, height = 20)
  plot(dtrain[,"Aff"], predsvmtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predsvmtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], predsvmtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predsvmtrain), digits = 3)))
  text(dtrain[,"Aff"], predsvmtrain, labels = names(predsvmtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(predsvmtrain)]
  plot(dtrain[,"Aff"], predsvmtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predsvmtrain), digits = 3)))
  text(dtrain[,"Aff"], predsvmtrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # test
  plot(dtest[,"Aff"], predsvmtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], predsvmtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], predsvmtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predsvmtest), digits = 3)))
  text(dtest[,"Aff"], predsvmtest, labels = names(predsvmtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclusterTest = dcluster[names(predsvmtest)]
  plot(dtest[,"Aff"], predsvmtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predsvmtest), digits = 3)))
  text(dtest[,"Aff"], predsvmtest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()
  
}



SVMTune = function(dtrain, vgamma, vcost, nbCV){
  
  
  lfolds = samplingDataNgroup(dtrain, nbCV)
  lmodel = list()
  lR2best = NULL
  
  k = 1
  kmax = length(lfolds)
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      lcpd = rownames(lfolds[[m]])
      if (m == k){
        dtest = as.data.frame(lfolds[[m]])
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    #dtrain = as.data.frame(scale(dtrain))
    ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
    #ddestrain = scale(ddestrain)
    Aff = dtrain[,c("Aff")]
    
    
    #dtest = as.data.frame(scale(dtest))
    dtestAff = dtest[,"Aff"]
    ddesctest = dtest[,-c(which(colnames(dtest) == "Aff"))]
    #ddesctest = scale(ddesctest, center = attr(ddestrain, 'scaled:center'), scale = attr(ddestrain, 'scaled:scale'))
    
    modelsvm = tune(svm, train.x = ddestrain, train.y = Aff, ranges = list(gamma = vgamma, cost = vcost), tunecontrol = tune.control(sampling = "cross"), kernel = "polynomial")
    modelsvm = modelsvm$best.model
    
    vpred = predict(modelsvm, ddesctest)
    
    R2 = calR2(dtestAff, vpred)
    #print(R2)
    lR2best = append(lR2best, R2)
    lmodel[[k]] =  modelsvm
    k = k + 1 
    
  }
  return(lmodel[[which(lR2best == max(lR2best))]])
}  



##################
# NEURAL NETWORK #
##################


NNRegCV = function(lfolds, dcluster, vdecay, vsize, prout){
  
  print(paste("==== NN in CV with ", length(lfolds), " Automatic optimization CV10 decay and size ====", sep = ""))
  
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      lcpd = rownames(lfolds[[m]])
      if (m == k){
        dtest = as.data.frame(lfolds[[m]])
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    # grid optimisation 
    modelNN = NNRegOptimizeGrid(dtrain, "0", vdecay, vsize, 10)
    
    #dtrain = as.data.frame(scale(dtrain))
    ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
    Aff = dtrain[,c("Aff")]
    
    #print (as.vector(Aff))
    ddestrain = scale(ddestrain)

    #dtest = as.data.frame(scale(dtest))
    dtestAff = dtest[,"Aff"]
    ddesctest = dtest[,-c(which(colnames(dtest) == "Aff"))]
    ddesctest = scale(ddesctest, center = attr(ddestrain, 'scaled:center'), scale = attr(ddestrain, 'scaled:scale'))
    
    vpred = predict(modelNN, ddesctest)
    names(vpred) = rownames(dtest)
    y_predict = append(y_predict, vpred)
    y_real = append(y_real, dtestAff)
    k = k + 1 
  }
  
  # performances
  valr2 = calR2(y_real, y_predict)
  corval = cor(y_real, y_predict)
  RMSEP = vrmsep(y_real, y_predict)
  MAEval = MAE(y_real, y_predict)
  R02val = R02(y_real, y_predict)
  
  print("Perfomances in CV")
  print(paste("R2=", valr2, sep = ""))
  print(paste("R02=", R02val, sep = ""))
  print(paste("MAE=", MAEval, sep = ""))
  print(paste("Cor=", corval, sep = ""))
  print(paste("RMSEP=", RMSEP, sep = ""))
  print("")
  print("")
  
  
  pdf(paste(prout, "PerfNNCV", length(lfolds), ".pdf", sep = ""), 20, 20)
  plot(y_real, y_predict, pch = 20, main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster
  dcluster = dcluster[names(y_predict)]
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()  
  
  tperf = cbind(y_predict, y_real)
  write.table(tperf, paste(prout, "perfNNCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
}




NNReg = function(dtrain, dtest, dcluster,  vdecay, vsize, prout){
  
  print(paste("==== Neural network in train-test ===", sep = ""))
  
  # find best model
  modelNN = NNRegOptimizeGrid(dtrain, "0", vdecay, vsize, 10)
  
  # scale data
  dtrain = as.data.frame(dtrain)
  ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
  dAfftrain = dtrain[,c("Aff")]
  ddestrain = scale(ddestrain)
  
  dtest = as.data.frame(dtest)
  ddestest = dtest[,-c(which(colnames(dtest) == "Aff"))]
  dAfftest = dtest[,c("Aff")]
  ddestest = scale(ddestest, attr(ddestrain, 'scaled:center'), scale = attr(ddestrain, 'scaled:scale'))
  
  # predict  #
  ############
  predNNtest = predict(modelNN, ddestest)
  predNNtrain = predict(modelNN, ddestrain)
  
  names(predNNtrain) = rownames(dtrain)
  names(predNNtest) = rownames(dtest)
  
  r2train = calR2(dtrain[,"Aff"], predNNtrain)
  r2test = calR2(dtest[,"Aff"], predNNtest)
  
  cortrain = cor(dtrain[,"Aff"], predNNtrain)
  cortest = cor(dtest[,"Aff"], predNNtest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predNNtrain)
  rmseptest = vrmsep(dtest[,"Aff"], predNNtest)
  
  R02train = R02(dtrain[,"Aff"], predNNtrain)
  R02test = R02(dtest[,"Aff"], predNNtest)
  
  MAEtrain = MAE(dtrain[,"Aff"], predNNtrain)
  MAEtest = MAE(dtest[,"Aff"], predNNtest)
  
  print("===== NN model train-Test =====")
  #print(modelpls$coefficients)
  print(paste("Perf training (dim= ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("R2 train=", r2train))
  print(paste("R02 train=", R02train))
  print(paste("MAE train=", MAEtrain))
  print(paste("Corval train=", cortrain))
  print(paste("RMSEP=", rmseptrain))
  print("")
  print("")
  
  
  print(paste("Perf test (dim=", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("R2 test=", r2test))
  print(paste("R02 test=", R02test))
  print(paste("MAE test=", MAEtest))
  print(paste("Corval test=", cortest, sep = ""))
  print(paste("RMSEP test=", rmseptest, sep = ""))
  print("")
  print("")
  
  
  # train
  pdf(paste(prout ,"PerfNNreg_TrainTest.pdf", sep = ""), width = 20, height = 20)
  plot(dtrain[,"Aff"], predNNtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], predNNtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)))
  text(dtrain[,"Aff"], predNNtrain, labels = names(predNNtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(predNNtrain)]
  plot(dtrain[,"Aff"], predNNtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)))
  text(dtrain[,"Aff"], predNNtrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # test
  plot(dtest[,"Aff"], predNNtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], predNNtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)))
  text(dtest[,"Aff"], predNNtest, labels = names(predNNtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclusterTest = dcluster[names(predNNtest)]
  plot(dtest[,"Aff"], predNNtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)))
  text(dtest[,"Aff"], predNNtest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()
  
}





NNRegOptimizeGrid = function(dtrain, pfig, vdecay, vsize, nbCV){
  
  
  lfolds = samplingDataNgroup(dtrain, nbCV)
  lmodel = list()
  lR2best = NULL
  
  if(pfig != "0"){
    pdf(paste(pfig, "_", length(lfolds), ".pdf", sep = ""), 20, 20)
  }
    
  k = 1
  kmax = length(lfolds)
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      lcpd = rownames(lfolds[[m]])
      if (m == k){
        dtest = as.data.frame(lfolds[[m]])
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    #dtrain = as.data.frame(scale(dtrain))
    ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
    Aff = dtrain[,c("Aff")]
    
    #print (as.vector(Aff))
    ddestrain = scale(ddestrain)
    
    #dtest = as.data.frame(scale(dtest))
    dtestAff = dtest[,"Aff"]
    ddesctest = dtest[,-c(which(colnames(dtest) == "Aff"))]
    ddesctest = scale(ddesctest, center = attr(ddestrain, 'scaled:center'), scale = attr(ddestrain, 'scaled:scale'))
    
    gridOptimize = NULL
    for (d in vdecay){
      optsize = NULL
      for(i in vsize){
        #2000
        modelNN = nnet(ddestrain, Aff, size = i, linout = T, maxit = 2000, MaxNWts = i*(dim(ddestrain)[2]+1)+i+1, decay = d)
        vpred = predict (modelNN, ddesctest)
        valr2 = calR2(dtestAff, vpred)
        if (valr2 < 0){
          valr2 = 0
        }
        optsize = append(optsize, valr2)
      }
      gridOptimize = cbind(gridOptimize, optsize)
    }
    colnames(gridOptimize) = vdecay
    rownames(gridOptimize) = vsize
    
    if(pfig != "0"){
      plot(rownames(gridOptimize), gridOptimize[,1] , type = "l", col = "red", lwd = 3, ylim = c(min(gridOptimize),max(gridOptimize)), xlab = "Size", ylab = "R2 in CV 10", main = paste("fold: ", k, sep = ""))
      lines(rownames(gridOptimize), gridOptimize[,2], lwd = 3, col = "green")
      lines(rownames(gridOptimize), gridOptimize[,3], lwd = 3, col = "blue")
      lines(rownames(gridOptimize), gridOptimize[,4], lwd = 3, col = "pink")
      lines(rownames(gridOptimize), gridOptimize[,5], lwd = 3, col = "yellow")
      lines(rownames(gridOptimize), gridOptimize[,6], lwd = 3, col = "black")
      legend("right" ,col = c("red", "green", "blue", "pink", "yellow", "black"), legend = c("0.001", "0.01", "0.1", "1", "10", "100"), pch = 19)
    }
    
    idecaybest = ceiling(which(gridOptimize == max(gridOptimize))/ dim(gridOptimize)[1])
    isizebest = which(gridOptimize[,idecaybest] == max(gridOptimize))
    
    bestsize = as.double(rownames(gridOptimize)[isizebest])
    bestdecay = as.double(colnames(gridOptimize)[idecaybest])
    
    # best model
    modelNN = nnet(ddestrain, Aff, size = bestsize, linout = T,  maxit = 2000, MaxNWts = bestsize*(dim(ddestrain)[2]+1)+bestsize+1, decay = bestdecay)
    vpred = predict(modelNN, ddesctest)
    R2 = calR2(dtestAff, vpred)
    lR2best = append(lR2best, R2)
    lmodel[[k]] =  modelNN
    k = k + 1 
  }
  
  if (pfig != "0"){
    dev.off() 
  }
  return(lmodel[[which(lR2best == max(lR2best))]])
}


#################
# DEEP LEARNING #
#################

DLRegCV = function(lfolds, dcluster, prout){
  
  # set seed and load library only if function is activated
  library(mxnet)
  require(mxnet)
  
  mx.set.seed(1)
  
  print(paste("==== DL in CV with ", length(lfolds), " Automatic optimization CV ====", sep = ""))
  
  
  # data combination #
  ####################
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  timportance = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      lcpd = rownames(lfolds[[m]])
      if (m == k){
        dtest = as.data.frame(lfolds[[m]])
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    ddesctrest = dtest[,-c(which(colnames(dtest) == "Aff"))]
    dtestAff = dtest[,c("Aff")]
    
    dtrain = as.data.frame(dtrain)
    ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
    dAff = dtrain[,c("Aff")]
    print (dAff)
    
    
    # number of hidden node
    data = mx.symbol.Variable("data")
    fc1 = mx.symbol.FullyConnected(data, num_hidden=)
    
    # Linear regression for output layer
    lro = mx.symbol.LinearRegressionOutput(fc1)
    
    mlpmodel <- mx.model.FeedForward.create(lro, X = ddestrain
                       ,Y = dAff[,1]
                       ,ctx=mx.gpu()
                       ,num.round=50
                       ,array.batch.size=20
                       ,learning.rate = 2e-6 #same as step size
                       ,eval.metric= mx.metric.mae)
    
    vpred = predict (dlmodel, dtest)
    print(vpred)
    names(vpred) = rownames(dtest)
    
    y_predict = append(y_predict, vpred)
    y_real = append(y_real, dtest[,"Aff"])
    k = k + 1
  }
  
  # performances
  valr2 = calR2(y_real, y_predict)
  corval = cor(y_real, y_predict)
  RMSEP = vrmsep(y_real, y_predict)
  MAEval = MAE(y_real, y_predict)
  R02val = R02(y_real, y_predict)
  
  print("Perfomances in CV")
  print(paste("R2=", valr2, sep = ""))
  print(paste("R02=", R02val, sep = ""))
  print(paste("MAE=", MAEval, sep = ""))
  print(paste("Cor=", corval, sep = ""))
  print(paste("RMSEP=", RMSEP, sep = ""))
  print("")
  print("")
  
  
  pdf(paste(prout, "PerfNNCV", length(lfolds), ".pdf", sep = ""), 20, 20)
  plot(y_real, y_predict, pch = 20, main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster
  dcluster = dcluster[names(y_predict)]
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()  
  
  tperf = cbind(y_predict, y_real)
  write.table(tperf, paste(prout, "perfNNCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
}







DLReg = function(dtrain, dtest, dcluster, prout){
  
  print(paste("==== Deep Learning in train-test ====", sep = ""))
  
  # set seed and load library only if function is activated
  library(mxnet)
  require(mxnet)
  
  mx.set.seed(1)
  
  dtrain = as.matrix(dtrain)
  ddestrain = dtrain[,-c(which(colnames(dtrain) == "Aff"))]
  dtrainAff = dtrain[,c("Aff")]
  
  dtest = as.matrix(dtest)
  ddestest = dtest[,-c(which(colnames(dtest) == "Aff"))]
  dtestAff = dtest[,c("Aff")]
  
  
  data = mx.symbol.Variable("data")
  label <- mx.symbol.Variable("label")
  
  fc1 <- mx.symbol.FullyConnected(data, num_hidden=5, name="fc1")
  tanh1 <- mx.symbol.Activation(fc1, act_type="tanh", name="tanh1")
  fc2 <- mx.symbol.FullyConnected(tanh1, num_hidden=5, name="fc2")
  tanh2 = mx.symbol.Activation(fc2, act_type="tanh", name="tanh1")
  fc3 <- mx.symbol.FullyConnected(tanh1, num_hidden=1, name="fc3")
  lro2 <- mx.symbol.LinearRegressionOutput(data=fc3, label=label, name="lro2")
  
  mx.set.seed(0)
  
  tuneNN <- mx.model.FeedForward.create(lro2, X=ddestrain, y=dtrainAff,
                                        ctx=mx.cpu(), num.round=4000, array.batch.size=20,
                                       learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse, verbose = FALSE)
  
    
  predNNtest = predict(tuneNN, ddestest)[1,]
  predNNtrain = predict(tuneNN, ddestrain)[1,]
  
  
  names(predNNtrain) = rownames(ddestrain)
  names(predNNtest) = rownames(ddestest)
  
  r2train = calR2(dtrain[,"Aff"], predNNtrain)
  r2test = calR2(dtest[,"Aff"], predNNtest)
  
  cortrain = cor(dtrain[,"Aff"], predNNtrain)
  cortest = cor(dtest[,"Aff"], predNNtest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predNNtrain)
  rmseptest = vrmsep(dtest[,"Aff"], predNNtest)
  
  R02train = R02(dtrain[,"Aff"], predNNtrain)
  R02test = R02(dtest[,"Aff"], predNNtest)
  
  MAEtrain = MAE(dtrain[,"Aff"], predNNtrain)
  MAEtest = MAE(dtest[,"Aff"], predNNtest)
  
  print("===== DL model train-Test =====")
  #print(modelpls$coefficients)
  print(paste("Perf training (dim= ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("R2 train=", r2train))
  print(paste("R02 train=", R02train))
  print(paste("MAE train=", MAEtrain))
  print(paste("Corval train=", cortrain))
  print(paste("RMSEP=", rmseptrain))
  print("")
  print("")
  
  
  print(paste("Perf test (dim=", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("R2 test=", r2test))
  print(paste("R02 test=", R02test))
  print(paste("MAE test=", MAEtest))
  print(paste("Corval test=", cortest, sep = ""))
  print(paste("RMSEP test=", rmseptest, sep = ""))
  print("")
  print("")
  
  
  # train
  pdf(paste(prout ,"PerfDLreg_TrainTest.pdf", sep = ""), width = 20, height = 20)
  plot(dtrain[,"Aff"], predNNtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], predNNtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)))
  text(dtrain[,"Aff"], predNNtrain, labels = names(predNNtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(predNNtrain)]
  plot(dtrain[,"Aff"], predNNtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], predNNtrain), digits = 3)))
  text(dtrain[,"Aff"], predNNtrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # test
  plot(dtest[,"Aff"], predNNtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], predNNtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)))
  text(dtest[,"Aff"], predNNtest, labels = names(predNNtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclusterTest = dcluster[names(predNNtest)]
  plot(dtest[,"Aff"], predNNtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], predNNtest), digits = 3)))
  text(dtest[,"Aff"], predNNtest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  dev.off()
  
}


#################
# Random forest #
#################

RFGridRegCV = function(lntree, lmtry, lfolds, prout){
  
  gridOpt = data.frame ()
  i = 0
  for (ntree in lntree){
    i = i + 1
    j = 0
    for (mtry in lmtry){
      j = j + 1
      
      # data combination
      k = 1
      kmax = length(lfolds)
      y_predict = NULL
      y_real = NULL
      while(k <= kmax){
        dtrain = NULL
        dtest = NULL
        for (m in seq(1:kmax)){
          if (m == k){
            dtest = lfolds[[m]]
          }else{
            dtrain = rbind(dtrain, lfolds[[m]])
          }
        }
        
        modelRF = randomForest( Aff~., data = dtrain, mtry=mtry, ntree = ntree, type = "response",  importance=TRUE)
        vpred = predict (modelRF, dtest, type = "response")
        
        y_predict = append(y_predict, vpred)
        y_real = append(y_real, dtest[,"Aff"])
        k = k + 1
      }
      
      # R2 for grid
      
      valr2 = calR2(y_real, y_predict)
      #print(valr2)
      gridOpt[i,j] = valr2
      
      # R conversion 
      
      #rate = calculTaux2  (changeList(as.double(l_predict)-1),changeList(y_real))
      #mcc = MCC (rate[1], rate[2], rate[3], rate[4])
      #grid[i,j] = mcc
    }
  }
  colnames (gridOpt) = lmtry
  rownames (gridOpt) = lntree
  
  write.table (gridOpt, paste(prout, "RFreg.grid", sep = ""))
  
  print(paste("=== RF grid optimisation in CV=", length(lfolds), " ntree = ", rownames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[1]], " mtry=", colnames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[2]], sep = ""))
  return (list(rownames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[1]],colnames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[2]] ))
}

RFregCV = function(lfolds, ntree, mtry, dcluster, prout){
  
  print(paste("==RF in CV with ", length(lfolds), " folds ntree = ", ntree, " mtry = ", mtry, sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  timportance = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modelRF = randomForest( Aff~., data = dtrain, mtry=as.integer(mtry), ntree=as.integer(ntree), type = "response",  importance=TRUE)
    vpred = predict (modelRF, dtest, type = "response")
    
    timportance = cbind(timportance, modelRF$importance[,1])
    
    names(vpred) = rownames(dtest)
    y_predict = append(y_predict, vpred)
    y_real = append(y_real, dtest[,"Aff"])
    k = k + 1
  }
  
  # performances
  valr2 = calR2(y_real, y_predict)
  corval = cor(y_real, y_predict)
  RMSEP = vrmsep(y_real, y_predict)
  MAEval = MAE(y_real, y_predict)
  R02val = R02(y_real, y_predict)
  
  print("== Perfomances in CV ==")
  print(paste("R2=", valr2, sep = ""))
  print(paste("R02=", R02val, sep = ""))
  print(paste("MAE=", MAEval, sep = ""))
  print(paste("Cor=", corval, sep = ""))
  print(paste("RMSEP=", RMSEP, sep = ""))
  print("")
  print("")
  
  
  # importance descriptors
  Mimportance = apply(timportance, 1, mean)
  SDimportance = apply(timportance, 1, sd)
  
  dimportance = cbind(Mimportance, SDimportance)
  rownames(dimportance) = rownames(timportance)
  colnames(dimportance) = c("M", "SD")
  dimportance = dimportance[order(dimportance[,1], decreasing = TRUE),]
  
  
  pdf(paste(prout, "PerfRFreg_CV", length(lfolds), ".pdf", sep = ""), 20, 20)
  par( mar=c(10,4,4,4))
  plot(dimportance[,1], xaxt ="n", xlab="", pch = 19, ylab="M importance")
  axis(1, 1:length(dimportance[,1]), labels = rownames(dimportance), las = 2, cex.axis = 0.7, cex = 2.75)
  for (i in 1:(dim(dimportance)[1])){
    segments(i, dimportance[i,1] - dimportance[i,2], i, dimportance[i,1] + dimportance[i,2])
  }
  
  plot(y_real, y_predict, pch = 20, main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster
  dcluster = dcluster[names(y_predict)]
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = dcluster, col = dcluster)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  dev.off()  
  
  write.table(dimportance, paste(prout, "ImportanceDescFRRegCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  
  # plot importance descriptors
  dimportance = read.table(paste(prout, "ImportanceDescFRRegCV_", length(lfolds), ".txt", sep = ""))
  
  ORDER = order(dimportance[,1], decreasing = T)
  NAME = rownames(dimportance)
  
  dimportance = cbind(dimportance, NAME)
  dimportance = cbind (dimportance, ORDER)
  dimportance = dimportance[ORDER[seq(1,10)],]
  dimportance = as.data.frame(dimportance)
  
  
  p = ggplot(dimportance, aes(-ORDER, M, fill = 1)) + 
    geom_bar(stat = "identity", show.legend = FALSE) + 
    scale_x_continuous(breaks = -dimportance$ORDER, labels = dimportance$NAME)+
    theme(axis.text.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 15, hjust = 0.5, vjust =0.1))+
    labs(y = "", x = "") + 
    ylim (c(0, 0.5)) +
    coord_flip()
  
  ggsave(paste(prout, "ImportanceDescRF_CV10.png", sep = ""), width = 6,height = 6, dpi = 300)
  
  
  tperf = cbind(y_predict, y_real)
  write.table(tperf, paste(prout, "perfRFRegCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  # test - plot for publication
  dpred = cbind(y_real, y_predict)
  colnames(dpred) = c("Yreal", "Ypredict")
  dpred = as.data.frame(dpred)
  
  p = ggplot(dpred, aes(Yreal, Ypredict))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=-2.1, y=2.2, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    xlim (c(4, 12)) +
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  ggsave(paste(prout, "PerfRFregpoint_CV10.png", sep = ""), width = 6,height = 6, dpi = 300)
  
  
  #dpred = cbind(dtest[,"Aff"], vpredtest)
  Vcluster = dcluster
  d = cbind(dpred, Vcluster)
  #colnames(dpred) = c("NAME", "Yreal", "Ypredict", "cluster")
  #dpred = as.data.frame(dpred)
  #print(dpred)
  
  p = ggplot(dpred, aes(Yreal, Ypredict, label=rownames(dpred)))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=-2.1, y=2.2, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    geom_text(size = 2.6, aes(label= paste(rownames(dpred), "^(", Vcluster, ")", sep = "")), parse = TRUE, color="black", nudge_y = 0.06) + 
    labs(x = "pAff", y = "Predict pAff") + 
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    xlim (c(4, 12)) +
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  print(p)
  ggsave(paste(prout, "PerfRFregname_CV10.png", sep = ""), width = 8,height = 8, dpi = 300)
}

RFreg = function (dtrain, dtest, ntree, mtry, dcluster, prout){
  
  
  modelRF = randomForest( Aff~., data = dtrain, mtry=as.integer(mtry), ntree=as.integer(ntree), type = "response",  importance=TRUE)
  vpredtrain = predict (modelRF, dtrain, type = "response")
  vpredtest = predict (modelRF, dtest, type = "response")
  
  names(vpredtrain) = rownames(dtrain)
  names(vpredtest) = rownames(dtest)
  
  testw = cbind(dtest[,"Aff"], vpredtest)
  colnames(testw) = c("Yreal", "Ypredict")
  write.csv(testw, paste(prout, "perfTestRFRegPred.csv"))

  trainw = cbind(dtrain[,"Aff"], vpredtrain)
  colnames(trainw) = c("Yreal", "Ypredict")
  write.csv(trainw, paste(prout, "perfTrainRFRegPred.csv"))  
    
  r2train = calR2(dtrain[,"Aff"], vpredtrain)
  cortrain = cor(dtrain[,"Aff"], vpredtrain)
  RMSEPtrain = vrmsep(dtrain[,"Aff"], vpredtrain)
  R02train = R02(dtrain[,"Aff"], vpredtrain)
  MAEtrain = MAE(dtrain[,"Aff"], vpredtrain)
  
  r2test = calR2(dtest[,"Aff"], vpredtest)
  cortest = cor(dtest[,"Aff"], vpredtest)
  RMSEPtest = vrmsep(dtest[,"Aff"], vpredtest)
  R02test = R02(dtest[,"Aff"], vpredtest)
  MAEtest = MAE(dtest[,"Aff"], vpredtest)
  
  
  print("===Perf RF===")
  print(paste("Dim train: ", dim(dtrain)[1]," ", dim(dtrain)[2], sep = ""))
  print(paste("Dim test: ", dim(dtest)[1]," ", dim(dtest)[2], sep = ""))
  
  print("==Train==")
  print(paste("R2 train=", r2train, sep = ""))
  print(paste("R02 train=", R02train, sep = ""))
  print(paste("MAE train=", MAEtrain, sep = ""))
  print(paste("cor train=", cortrain, sep = ""))
  print(paste("RMSEP train=", RMSEPtrain, sep = ""))
  
  
  print("==Test==")
  print(paste("R2 test=", r2test, sep = ""))
  print(paste("R02 test=", R02test, sep = ""))
  print(paste("MAE test=", MAEtest, sep = ""))
  print(paste("cor test=", cortest, sep = ""))
  print(paste("RMSEP test=", RMSEPtest, sep = ""))
  print("")
  print("")
  
  
  pdf(paste(prout, "PerfRFreg_TrainTest.pdf", sep = ""), 20, 20)
  plot(modelRF)
  
  # train
  plot(dtrain[,"Aff"], vpredtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], vpredtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], vpredtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], vpredtrain), digits = 3)))
  text(dtrain[,"Aff"], vpredtrain, labels = names(vpredtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclusterTrain = dcluster[names(vpredtrain)]
  text(dtrain[,"Aff"], vpredtrain, labels = dclusterTrain, col = dclusterTrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  
  # test
  plot(dtest[,"Aff"], vpredtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], vpredtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], vpredtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], vpredtest), digits = 3)))
  text(dtest[,"Aff"], vpredtest, labels = names(vpredtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster-test
  dclusterTest = dcluster[names(vpredtest)]
  text(dtest[,"Aff"], vpredtest, labels = dclusterTest, col = dclusterTest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  dev.off() 
  
  
  # test - plot for publication
  dpred = cbind(dtest[,"Aff"], vpredtest)
  colnames(dpred) = c("Yreal", "Ypredict")
  dpred = as.data.frame(dpred)
  
  p = ggplot(dpred, aes(Yreal, Ypredict))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=5, y=11.5, label = paste("R2=",round(r2test,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    xlim (c(4, 12)) +
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  #print(p)
  ggsave(paste(prout, "PerfRFregpoint_Test.png", sep = ""), width = 6,height = 6, dpi = 300)
  
  
  #dpred = cbind(dtest[,"Aff"], vpredtest)
  Vcluster = dclusterTest
  #dpred = cbind(dpred, as.character(dclusterTest))
  dpred = cbind(dpred, Vcluster)
  #colnames(dpred) = c("NAME", "Yreal", "Ypredict", "cluster")
  #dpred = as.data.frame(dpred)
  #print(dpred)
  
  p = ggplot(dpred, aes(Yreal, Ypredict, label=rownames(dpred)))+
    geom_point(size=1.5, colour="black", shape=21) + 
    geom_text(x=5, y=11.5, label = paste("R2=",round(r2test,2), sep = ""), size = 8)+
    labs(x = "pAff", y = "Predicted pAff") +
    geom_text(size = 2.6, aes(label= paste(rownames(dpred), "^(", Vcluster, ")", sep = "")), parse = TRUE, color="black", nudge_y = 0.06) + 
    labs(x = "Real pAff", y = "Predict pAff") + 
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    xlim (c(4, 12)) +
    geom_segment(aes(x = 4, y = 4, xend = 12, yend = 12), linetype=2, size = 0.1) + 
    ylim (c(4, 12)) 
  #print(p)
  ggsave(paste(prout, "PerfRFregname_Test.png", sep = ""), width = 8,height = 8, dpi = 300)
  
}


#############
#   CART    #
#############


CARTRegCV = function(lfolds, dcluster, prout){
  
  print(paste("== CART in CV with ", length(lfolds), "==", sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  y_proba = NULL
  #print (paste(prout, "TreeCARTReg-CV", length(lfolds), ".pdf",sep = ""))
  pdf(paste(prout, "TreeCARTReg-CV", length(lfolds), ".pdf",sep = ""))
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modelCART = rpart( Aff~., data = dtrain, method = "anova")#, control = rpart.control(cp = 0.05))
    vpred = predict(modelCART, dtest[,-dim(dtest)[2]], type = "vector")
    names(vpred) = rownames(dtest)
    
    #print(vpred)
    # plot tree in pdf
    plotcp(modelCART)
    rpart.plot( modelCART , # middle graph
                box.palette="GnBu",
                branch.lty=3, shadow.col="gray", nn=TRUE)
    
    y_predict = append(y_predict, vpred)
    y_real = append(y_real, dtest[,"Aff"])
    
    k = k + 1
  }
  dev.off()
  
  # performances
  valr2 = calR2(y_real, y_predict)
  corval = cor(y_real, y_predict)
  RMSEP = vrmsep(y_real, y_predict)
  R02val = R02(y_real, y_predict)
  MAEval = MAE(y_real, y_predict)
  
  print("=== Perfomances in CV ===")
  print(paste("R2=", valr2, sep = ""))
  print(paste("R02=", R02val, sep = ""))
  print(paste("MAE=", MAEval, sep = ""))
  print(paste("Cor=", corval, sep = ""))
  print(paste("RMSEP=", RMSEP, sep = ""))
  
  
  pdf(paste(prout, "PerfCARTReg_CV", length(lfolds), ".pdf", sep = ""), 20, 20)
  plot(y_real, y_predict, pch = 20, main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster
  dcluster = dcluster[names(y_predict)]
  plot(y_real, y_predict, type = "n", main = paste("Correlation = ", round(cor(y_real, y_predict), digits = 3)))
  text(y_real, y_predict, labels = names(y_predict))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  dev.off()  
  
  dpred = cbind(y_predict, y_real)
  colnames(dpred) = c("Predict", "Real")
  write.table(dpred, file = paste(prout, "PerfCARTRegCV", length(lfolds), ".txt", sep = ""), sep = "\t")
  
}


CARTreg = function (dtrain, dtest, dcluster, prout){
  
  modelCART = rpart( Aff~., data = dtrain, method = "anova")#, control = rpart.control(cp = 0.05))
  
  pdf(paste(prout, "TreeCARTReg-Train.pdf",sep = ""))
  plotcp(modelCART)
  rpart.plot( modelCART , # middle graph
              box.palette="GnBu",
              branch.lty=3, shadow.col="gray", nn=TRUE)
  dev.off()
  
  vpredtest = predict(modelCART, dtest[,-dim(dtest)[2]], type = "vector") # remove affinity
  vpredtrain = predict(modelCART, dtrain[,-dim(dtrain)[2]], type = "vector") # remove affinity
  
  names(vpredtest) = rownames(dtest)
  names(vpredtrain) = rownames(dtrain)
  
  write.csv(vpredtest, paste(prout, "perfCARTPred.csv"))
  
  r2train = calR2(dtrain[,"Aff"], vpredtrain)
  cortrain = cor(dtrain[,"Aff"], vpredtrain)
  RMSEPtrain = vrmsep(dtrain[,"Aff"], vpredtrain)
  R02train = R02(dtrain[,"Aff"], vpredtrain)
  MAEtrain = MAE(dtrain[,"Aff"], vpredtrain)
  
  r2test = calR2(dtest[,"Aff"], vpredtest)
  cortest = cor(dtest[,"Aff"], vpredtest)
  RMSEPtest = vrmsep(dtest[,"Aff"], vpredtest)
  R02test = R02(dtest[,"Aff"], vpredtest)
  MAEtest = MAE(dtest[,"Aff"], vpredtest)
  
  print("===Perf CART===")
  print(paste("Dim train: ", dim(dtrain)[1]," ", dim(dtrain)[2], sep = ""))
  print(paste("Dim test: ", dim(dtest)[1]," ", dim(dtest)[2], sep = ""))
  
  print("==Train==")
  print(paste("R2 train=", r2train, sep = ""))
  print(paste("R02 train=", R02train, sep = ""))
  print(paste("MAE train=", MAEtrain, sep = ""))
  print(paste("cor train=", cortrain, sep = ""))
  print(paste("RMSEP train=", RMSEPtrain, sep = ""))
  
  
  print("==Test==")
  print(paste("R2 test=", r2test, sep = ""))
  print(paste("R02 test=", R02test, sep = ""))  
  print(paste("MAE test=", MAEtest, sep = ""))
  print(paste("cor test=", cortest, sep = ""))
  print(paste("RMSEP test=", RMSEPtest, sep = ""))
  
  pdf(paste(prout, "PerfCARTreg_TrainTest.pdf", sep = ""), 20, 20)
  
  # train
  plot(dtrain[,"Aff"], vpredtrain, pch = 20, main = paste("Correlation = ", round(cor(dtrain[,"Aff"], vpredtrain), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtrain[,"Aff"], vpredtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], vpredtrain), digits = 3)))
  text(dtrain[,"Aff"], vpredtrain, labels = names(vpredtrain))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - train
  dclustertrain = dcluster[names(vpredtrain)]
  plot(dtrain[,"Aff"], vpredtrain, type = "n", main = paste("Correlation = ", round(cor(dtrain[,"Aff"], vpredtrain), digits = 3)))
  text(dtrain[,"Aff"], vpredtrain, labels = dclustertrain, col = dclustertrain)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  
  # test
  plot(dtest[,"Aff"], vpredtest, pch = 20, main = paste("Correlation = ", round(cor(dtest[,"Aff"], vpredtest), digits = 3)), cex = 2)
  abline(a = 0, b = 1, col = "red", cex = 3)
  plot(dtest[,"Aff"], vpredtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], vpredtest), digits = 3)))
  text(dtest[,"Aff"], vpredtest, labels = names(vpredtest))
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  # cluster - test
  dclustertest = dcluster[names(vpredtest)]
  plot(dtest[,"Aff"], vpredtest, type = "n", main = paste("Correlation = ", round(cor(dtest[,"Aff"], vpredtest), digits = 3)))
  text(dtest[,"Aff"], vpredtest, labels = dclustertest, col = dclustertest)
  abline(a = 0, b = 1, col = "red", cex = 3)
  
  dev.off() 
  
}




######################################################################################################################################
######################################################################################################################################


#####################
#  classification   #
#####################

##########
#  SVM   #
##########

SVMClassCV = function(lfolds, vgamma, vcost, prout){
  
  print(paste("== SVM in CV with ", length(lfolds), " Automatic optimization by folds", sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  y_proba = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modtune = tune(svm, Aff~., data = dtrain, ranges = list(gamma = vgamma, cost = vcost), tunecontrol = tune.control(sampling = "fix"))
    
    vpred = predict (modtune$best.model, dtest, type = "class")
    y_proba = append(y_proba, vpred)
    
    vpred[which(vpred < 0.5)] = 0
    vpred[which(vpred >= 0.5)] = 1
    
    y_predict = append(y_predict, vpred)
    y_real = append(y_real, dtest[,"Aff"])
    k = k + 1
  }
  
  # performances
  lpref = classPerf(y_real, y_predict)
  acc = lpref[[1]]
  se = lpref[[2]]
  sp = lpref[[3]]
  mcc = lpref[[4]]
  
  png(paste(prout, "PerfSVMClassCV", length(lfolds), ".png", sep = ""), 800, 800)
  plot(y_real, y_proba, type = "n")
  text(y_real, y_proba, labels = names(y_predict), cex = 0.8)
  abline(a = 0.5, b = 0, col = "red", cex = 3)
  dev.off()
  
  dpred = cbind(y_proba, y_real)
  colnames(dpred) = c("Predict", "Real")
  write.table(dpred, file = paste(prout, "PerfRFClassCV", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  print("Perfomances in CV")
  print(paste("acc=", acc, sep = ""))
  print(paste("se=", se, sep = ""))
  print(paste("sp=", sp, sep = ""))
  print(paste("mcc=", mcc, sep = "")) 
  
  tperf = cbind(y_predict, y_real)
  write.table(tperf, paste(prout, "perfSVMRegCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
}




#####################
#  Random forest    #
#####################

RFGridClassCV = function(lntree, lmtry, lfolds, prout){
  
  gridOpt = data.frame ()
  i = 0
  for (ntree in lntree){
    i = i + 1
    j = 0
    for (mtry in lmtry){
      j = j + 1
      
      # data combination
      k = 1
      kmax = length(lfolds)
      y_predict = NULL
      y_real = NULL
      while(k <= kmax){
        dtrain = NULL
        dtest = NULL
        for (m in seq(1:kmax)){
          if (m == k){
            dtest = lfolds[[m]]
          }else{
            dtrain = rbind(dtrain, lfolds[[m]])
          }
        }
        
        modelRF = randomForest( Aff~., data = dtrain, mtry=mtry, ntree = ntree, type = "class",  importance=TRUE)
        vpred = predict (modelRF, dtest)
        vpred[which(vpred < 0.5)] = 0
        vpred[which(vpred >= 0.5)] = 1
        
        y_predict = append(y_predict, vpred)
        y_real = append(y_real, dtest[,"Aff"])
        k = k + 1
      }
      
      
      # R2 for grid
      #print(y_predict)
      l_perf = classPerf(y_real, y_predict)
      gridOpt[i,j] = l_perf[[4]]
      
      # R conversion 
    }
  }
  colnames (gridOpt) = lmtry
  rownames (gridOpt) = lntree
  
  write.table (gridOpt, paste(prout, "RFclassMCC.grid", sep = ""))
  print(which(gridOpt == max(gridOpt), arr.ind = TRUE))
  
  print(paste("=== RF grid optimisation in CV = ", length(lfolds), " ntree = ", rownames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[1]], " mtry=", colnames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[2]], sep = ""))
  return (list(rownames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[1]],colnames (gridOpt)[which(gridOpt==max(gridOpt), arr.ind=T)[2]] ))
}

RFClassCV = function(lfolds, ntree, mtry, prout){
  
  print(paste("== RF in CV with ", length(lfolds), " folds ntree = ", ntree, " mtry = ", mtry, sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  timportance = NULL
  y_proba = NULL
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modelRF = randomForest( Aff~., data = dtrain, mtry=as.integer(mtry), ntree=as.integer(ntree), type = "class",  importance=TRUE)
    vpred = predict (modelRF, dtest)
    vproba = vpred
    vpred[which(vpred < 0.5)] = 0
    vpred[which(vpred >= 0.5)] = 1
    
    timportance = cbind(timportance, modelRF$importance[,1])
    y_predict = append(y_predict, vpred)
    y_proba = append(y_proba, vproba)
    y_real = append(y_real, dtest[,"Aff"])
    k = k + 1
  }
  
  # performances
  lpref = classPerf(y_real, y_predict)
  acc = lpref[[1]]
  se = lpref[[2]]
  sp = lpref[[3]]
  mcc = lpref[[4]]
  
  png(paste(prout, "PerfRFClassCV", length(lfolds), ".png", sep = ""), 800, 800)
  plot(y_real, y_proba, type = "n")
  text(y_real, y_proba, labels = names(y_predict), cex = 0.8)
  abline(a = 0.5, b = 0, col = "red", cex = 3)
  dev.off()
  
  dpred = cbind(y_proba, y_real)
  colnames(dpred) = c("Predict", "Real")
  write.table(dpred, file = paste(prout, "PerfRFClassCV", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  print("Perfomances in CV")
  print(paste("acc=", acc, sep = ""))
  print(paste("se=", se, sep = ""))
  print(paste("sp=", sp, sep = ""))
  print(paste("mcc=", mcc, sep = ""))
  
  # importance descriptors
  Mimportance = apply(timportance, 1, mean)
  SDimportance = apply(timportance, 1, sd)
  
  dimportance = cbind(Mimportance, SDimportance)
  rownames(dimportance) = rownames(timportance)
  dimportance = dimportance[order(dimportance[,1], decreasing = TRUE),]
  
  write.table(dimportance, paste(prout, "ImportanceDescClassCV_", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  png(paste(prout, "ImportanceRFClassCV_", length(lfolds), ".png", sep = ""), 1000, 800)
  par( mar=c(10,4,4,4))
  plot(dimportance[,1], xaxt ="n", xlab="", pch = 19, ylab="M importance")
  axis(1, 1:length(dimportance[,1]), labels = rownames(dimportance), las = 2, cex.axis = 0.7, cex = 2.75)
  for (i in 1:(dim(dimportance)[1])){
    segments(i, dimportance[i,1] - dimportance[i,2], i, dimportance[i,1] + dimportance[i,2])
  }
  dev.off()
  return(y_predict)
}


RFClass = function (dtrain, dtest, ntree, mtry, prout){
  
  modelRF = randomForest( Aff~., data = dtrain, mtry=as.integer(mtry), ntree=as.integer(ntree), type = "class",  importance=TRUE)
  vpredtrain = predict (modelRF, dtrain, type = "class")
  vpredtest = predict (modelRF, dtest, type = "class")
  
  vpredtrainprob = vpredtrain
  vpredtestprob = vpredtest
  
  vpredtest[which(vpredtest < 0.5)] = 0
  vpredtest[which(vpredtest >= 0.5)] = 1
  vpredtrain[which(vpredtrain < 0.5)] = 0
  vpredtrain[which(vpredtrain >= 0.5)] = 1
  
  
  vperftrain = classPerf(dtrain[,c("Aff")], vpredtrain)
  vperftest = classPerf(dtest[,c("Aff")], vpredtest)
  
  print("===Perf RF===")
  print(paste("Dim train: ", dim(dtrain)[1]," ", dim(dtrain)[2], sep = ""))
  print(paste("Dim test: ", dim(dtest)[1]," ", dim(dtest)[2], sep = ""))
  
  print("==Train==")
  print(paste("acc=", vperftrain[[1]], sep = ""))
  print(paste("se=", vperftrain[[2]], sep = ""))
  print(paste("sp=", vperftrain[[3]], sep = ""))
  print(paste("mcc=", vperftrain[[4]], sep = ""))
    
  
  print("==Test==")
  print(paste("acc=", vperftest[[1]], sep = ""))
  print(paste("se=", vperftest[[2]], sep = ""))
  print(paste("sp=", vperftest[[3]], sep = ""))
  print(paste("mcc=", vperftest[[4]], sep = ""))
  

  png(paste(prout, "PerfTrainTest.png", sep = ""), 1600, 800)
  par(mfrow = c(1,2))
  plot(dtrain[,"Aff"], vpredtrainprob, type = "n")
  text(dtrain[,"Aff"], vpredtrainprob, labels = names(vpredtrainprob))
  abline(a = 0.5, b = 0, col = "red", cex = 3)
  
  plot(dtest[,"Aff"], vpredtest, type = "n")
  text(dtest[,"Aff"], vpredtestprob, labels = names(vpredtestprob))
  abline(a = 0.5, b = 0, col = "red", cex = 3)
  dev.off()
  write.csv(vpredtestprob, paste(prout,"classTest.csv", sep = ""))
  
}



##########
#  CART  #
##########

CARTClassCV = function(lfolds, prout){
  
  print(paste("== CART in CV with ", length(lfolds), "==", sep = ""))
  
  # data combination
  k = 1
  kmax = length(lfolds)
  y_predict = NULL
  y_real = NULL
  y_proba = NULL
  pdf(paste(prout, "TreeCARTClass-CV", length(lfolds), ".pdf",sep = ""))
  while(k <= kmax){
    dtrain = NULL
    dtest = NULL
    for (m in seq(1:kmax)){
      if (m == k){
        dtest = lfolds[[m]]
      }else{
        dtrain = rbind(dtrain, lfolds[[m]])
      }
    }
    
    modelCART = rpart( Aff~., data = dtrain, method = "class")
    vpred = predict(modelCART, dtest)
    vpred = vpred[,2]
    vproba = vpred
    vpred[which(vpred < 0.5)] = 0
    vpred[which(vpred >= 0.5)] = 1
    
    plotcp(modelCART)
    # plot tree in pdf
    rpart.plot( modelCART , # middle graph
               extra=104, box.palette="GnBu",
               branch.lty=3, shadow.col="gray", nn=TRUE)
    
    y_predict = append(y_predict, vpred)
    y_proba = append(y_proba, vproba)
    y_real = append(y_real, dtest[,"Aff"])
    
    k = k + 1
  }
  dev.off()
  
  # performances
  lpref = classPerf(y_real, y_predict)
  acc = lpref[[1]]
  se = lpref[[2]]
  sp = lpref[[3]]
  mcc = lpref[[4]]
  
  png(paste(prout, "PerfCARTClassCV", length(lfolds), ".png", sep = ""), 800, 800)
  plot(y_real, y_proba, type = "n")
  text(y_real, y_proba, labels = names(y_predict), cex = 0.8)
  abline(a = 0.5, b = 0, col = "red", cex = 3)
  dev.off()
  
  dpred = cbind(y_proba, y_real)
  colnames(dpred) = c("Predict", "Real")
  write.table(dpred, file = paste(prout, "PerfCARTClassCV", length(lfolds), ".txt", sep = ""), sep = "\t")
  
  print("Perfomances in CV")
  print(paste("acc=", acc, sep = ""))
  print(paste("se=", se, sep = ""))
  print(paste("sp=", sp, sep = ""))
  print(paste("mcc=", mcc, sep = ""))
  print("")
  print("")
  
}

