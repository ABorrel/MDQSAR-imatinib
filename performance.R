


#####################
# Perf for classes  #
#####################

# change case of d or nd
perftable = function (list_predict, list_real){
  
  nb_value = length (list_real)
  i = 1
  tp = 0
  fp = 0
  tn = 0
  fn = 0
  while(i <= nb_value ){
    if (list_predict[i]==1){
      if (list_predict[i] == list_real[i]){
        tp = tp + 1
      }else {
        fp = fp + 1
      }
    }else{
      if (list_predict[i] == list_real[i]){
        tn = tn + 1
      }else {
        fn = fn + 1
      }
    }
    i = i + 1
  }
  #print (paste ("TP : ", tp, sep = ""))
  #print (paste ("TN : ", tn, sep = ""))
  #print (paste ("FP : ", fp, sep = ""))
  #print (paste ("FN : ", fn, sep = ""))
  
  tableval = c(tp,tn,fp,fn)
  return (tableval)
}



accuracy = function (tp, tn, fp, fn){
  return ((tp + tn)/(tp + fp + tn +fn))
}

precision = function (tp, fp){
  return (tp/(tp + fp))
}

recall = function (tp, fn){
  return (tp/(tp + fn))
}

specificity = function (tn, fp){
  return (tn/(tn + fp))
}

sensibility = function (tp, fn){
  return (tp/(tp + fn))
}


BCR = function (tp, tn, fp, fn){
  return (0.5*(tp/(tp+fn) + tn/(tn+fp)))
  
}

MCC = function (tp, tn, fp, fn){
  numerator = tp*tn-fp*fn
  denumerator = (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
  return (numerator / sqrt(denumerator))
}


qualityPredict = function (predict, Y2){
  print (as.vector(predict)[[1]])
  print (as.vector(Y2)[[1]])
  v_predict = calculTaux (as.vector(predict)[[1]], as.vector(Y2)[[1]])
  print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
  #print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  return (v_predict)
}




qualityPredictList = function (test_vector, real_vector){
  v_predict = calculTaux (test_vector, real_vector)
  print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
  #print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  return (v_predict)
}


qualityShowModelSelection = function (list_des_model, coef, v_real_train, v_predict_train, v_real_test, v_predict_test, v_real_loo, v_predict_loo, l_out_CV){
  
  # loo
  criteria_loo = computedCriteria(v_real_loo, v_predict_loo)
  # CV
  criteria_CV = computedCriteriaCV(l_out_CV)
  # train
  criteria_train = computedCriteria(v_real_train, v_predict_train)
  # test
  criteria_test = computedCriteria(v_real_test, v_predict_test)
  
  # show
  print ("descriptor")
  print (list_des_model)
  print (as.vector(abs(coef[list_des_model])))
  print ("Acc_loo --- Acc_train --- Acc_test --- Acc_CV_train --- SD_CV_train --- Acc_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[1]], criteria_train[[1]], criteria_test[[1]], criteria_CV["acc_train"], criteria_CV["acc_train_SD"], criteria_CV["acc_test"], criteria_CV["acc_test_SD"], sep = "---"))
  print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test --- Se_CV_train --- SD_CV_train --- Se_CV_test --- SD_CV_test --- Sp_CV_train --- SD_CV_train --- Sp_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[2]], criteria_loo[[3]], criteria_train[[2]], criteria_train[[3]], criteria_test[[2]], criteria_test[[3]], criteria_CV["se_train"], criteria_CV["se_train_SD"], criteria_CV["se_test"],criteria_CV["se_test_SD"], criteria_CV["sp_train"], criteria_CV["sp_train_SD"], criteria_CV["sp_test"],criteria_CV["sp_test_SD"], sep = "---"))
  print ("MCC_loo --- MCC_train --- MCC_test --- MCC_CV_train --- SD_CV_train --- MCC_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[4]], criteria_train[[4]], criteria_test[[4]], criteria_CV["mcc_train"], criteria_CV["mcc_train_SD"], criteria_CV["mcc_test"], criteria_CV["mcc_test_SD"], sep = "---"))
  print ("**********************************************************************")
}



classPerf = function (v_real, v_predict){
  
  rate = perftable (v_predict, v_real)
  acc = accuracy(rate[1], rate[2], rate[3], rate[4])
  se = sensibility(rate[1], rate[4])
  sp = sensibility(rate[2], rate[3])
  mcc = MCC(rate[1], rate[2], rate[3], rate[4])
  
  return (list (acc, se, sp, mcc))
  
}


computedCriteriaCV = function (l_out_CV){
  
  CV_train = l_out_CV[[1]]
  CV_test = l_out_CV[[2]]
  
  v_acc_train = NULL
  v_acc_test = NULL
  v_se_train = NULL
  v_se_test = NULL
  v_sp_train = NULL
  v_sp_test = NULL
  v_mcc_train = NULL
  v_mcc_test = NULL
  
  
  for (i in seq (1,length (CV_train))){
    v_acc_train = append (v_acc_train, CV_train[[i]][[1]])
    v_acc_test = append (v_acc_test, CV_test[[i]][[1]])
    v_se_train = append (v_se_train, CV_train[[i]][[2]])
    v_se_test = append (v_se_test, CV_test[[i]][[2]])
    v_sp_train = append (v_sp_train, CV_train[[i]][[3]])
    v_sp_test = append (v_sp_test, CV_test[[i]][[3]])
    v_mcc_train = append (v_mcc_train, CV_train[[i]][[4]])
    v_mcc_test = append (v_mcc_test, CV_test[[i]][[4]])
  }
  
  v_out = c(mean (v_acc_train), sd (v_acc_train), mean (v_acc_test), sd (v_acc_test), mean (v_se_train), sd (v_se_train), mean (v_se_test), sd (v_se_test), mean (v_sp_train), sd (v_sp_train), mean (v_sp_test), sd (v_sp_test),mean (v_mcc_train), sd (v_mcc_train), mean (v_mcc_test), sd (v_mcc_test) )
  names (v_out) = c("acc_train", "acc_train_SD","acc_test", "acc_test_SD","se_train", "se_train_SD", "se_test", "se_test_SD", "sp_train", "sp_train_SD", "sp_test", "sp_test_SD", "mcc_train", "mcc_train_SD", "mcc_test", "mcc_test_SD")
  return (v_out)
}



cumulTaux = function (taux1, taux2){
  
  tp = taux1[1] + taux2[1]
  tn = taux1[2] + taux2[2]
  fp = taux2[3]
  fn = taux2[4]
  
  print (paste ("accuracy : ", accuracy(tp, tn, fp, fn), sep = ""))
  print (paste ("precision : ",precision(tp, fp), sep = ""))
  #print (paste ("recall : ", recall(tp, fn), sep = ""))
  print (paste ("sensibility : ", sensibility(tp, fn), sep = ""))
  print (paste ("specificity : ", sensibility(tn, fp), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (tp, tn, fp, fn), sep = ""))	
  print (paste ("BER (balanced error rate) : ", 1 - BCR (tp, tn, fp, fn), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (tp, tn, fp, fn), sep = ""))	
  
}


# for ROC curve -> calcul vecteur prediction with probability (just for druggability)
generateVect = function(proba_out_predict, threshold){
  
  proba_class1 = proba_out_predict[,1]
  
  vect_out = NULL
  
  for (proba in proba_class1){
    if (proba > threshold){
      vect_out = c(vect_out, "d")
    }else{
      vect_out = c(vect_out, "nd")
    }
  }
  return (vect_out)
}


#########################
#    PERF regression    #
#########################


vrmsep = function(dreal, dpredict){
  
  #dpredict = dpredict[rownames(dreal),]
  
  i = 1
  imax = length(dreal)
  
  valout = 0
  while(i <= imax){
    valout = valout + ((dreal[i] - dpredict[i])^2)
    i = i + 1
  }
  return(sqrt(valout))
  
}



calR2 = function(dreal, dpredict){
  
  dreal = as.vector(dreal)
  dpredict = as.vector(dpredict)
  
  #print("Nb val in perf:")
  #print(length(dreal))
  
  dperf = cbind(dreal, dpredict)
  dperf = na.omit(dperf)
  
  #print("Nb val predict:")
  #print(dim(dperf))
  
  M = mean(dperf[,1])
  SCEy = 0
  SCEtot = 0
  for (i in seq(1, dim(dperf)[1])){
    #print (i)
    SCEy = SCEy + (dperf[i, 1] - dperf[i, 2])*(dperf[i, 1] - dperf[i, 2])
    SCEtot = SCEtot + (dperf[i, 1] - M)*(dperf[i, 1] - M)
  }
  
  r2 = 1 - (SCEy/SCEtot)
  return (as.double(r2))
}


MAE = function(dreal, dpredict){
  
  #dpredict = dpredict[rownames(dreal),]
  
  i = 1
  imax = length(dreal)
  
  valout = 0
  while(i <= imax){
    valout = valout + (abs(dreal[i] - dpredict[i]))
    i = i + 1
  }
  return(valout/imax)
}

R02 = function(dreal, dpredict){
  
  dreal = as.vector(dreal)
  dpredict = as.vector(dpredict)
  
  #print("Nb val in perf:")
  #print(length(dreal))
  
  dperf = cbind(dreal, dpredict)
  dperf = na.omit(dperf)
  
  #print("Nb val predict:")
  #print(dim(dperf))
  
  Mreal = mean(dperf[,1])
  Mpredict = mean(dperf[,2])
  
  #print(paste("Mpred - ",Mpredict))
  
  A = 0
  B = 0
  k = 0
  yypred = 0
  Sumpredict = 0
  # first loop for k
  for (i in seq(1, dim(dperf)[1])){
    #print (i)
    yypred = yypred + (dperf[i,1]*dperf[i,2])
    Sumpredict = Sumpredict + (dperf[i,2]^2)
  }
  #print(yypred)
  #print(Sumpredict)
  k = yypred/Sumpredict
  #print(paste("k - ", k))
  
  for (i in seq(1, dim(dperf)[1])){
    #print (i)
    tempA = ((dperf[i,2]-(k*dperf[i,2]))^2)
    tempB = ((dperf[i,2]-Mpredict)^2)
    #print(paste(tempA, tempB))
    
    A = A + ((dperf[i,2]-(k*dperf[i,2]))^2)
    B = B + ((dperf[i,2]-Mpredict)^2)
  }
  #print(k)
  #print(paste("A -", A))
  #print(paste("B -",B))
  
  r02 = as.double(A/B)
  return (1 - r02)
}
