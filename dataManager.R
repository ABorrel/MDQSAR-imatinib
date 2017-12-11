#!/usr/bin/env Rscript
source("PCAplot.R")
library("ggplot2")


###############################
# divise the dataset in folds #
###############################

samplingDataNgroupClass = function (t_din, i_nb_group, s_nameclass){
  
  # divise two classes
  v_class = as.factor(t_din[,s_nameclass])
  t_dc0 = t_din [which(v_class == 0),]
  t_dc1 = t_din [which(v_class == 1),]
  
  
  # sample data
  v_sampledc0 = sample (dim (t_dc0)[1])
  v_sampledc1 = sample (dim (t_dc1)[1])
  
  # ind limit
  i_limitc0 = as.integer (dim(t_dc0)[1] / i_nb_group)
  i_limitc1 = as.integer (dim(t_dc1)[1] / i_nb_group)
  
  #print (i_limitc0)
  #print (i_limitc1)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = rbind (t_dc0[v_sampledc0[1:i_limitc0],], t_dc1[v_sampledc1[1:i_limitc1],])
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = rbind (t_dc0[v_sampledc0[((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))],], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(length (v_sampledc1))]),])
    }
    else{
      t_group = rbind (t_dc0[(v_sampledc0[(i_limitc0 * (g-1) + 1):(i_limitc0 * g)]),], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(i_limitc1 * g)]),])
    }
    # append list
    output[[g]] = t_group
  }
  return (output)
}

sampligDataFractionCluster = function(t_din, fract, pcluster){
  
  dclust = read.csv(pcluster, sep=",", header = TRUE)
  print (t_din)
  lclust = unique(dclust[,2])
  
  dtrain = NULL
  dtest = NULL
  for (clust in lclust){
    lcpIDclust = dclust[which(dclust[,2] == clust),1]
    descclust = t_din[lcpIDclust,]
    
    # sample data
    v_sample = sample (dim (descclust)[1])
    
    # ind limit
    i_limitc = round((dim (descclust)[1]) * fract)
    
    dtrain = rbind(dtrain, descclust[v_sample[(i_limitc + 1):(length (v_sample))],])
    dtest = rbind(dtest, descclust[v_sample[1:i_limitc],])
  }
  # sample data
  return (list(dtrain, dtest))
  
}


samplingDataFraction = function (t_din, fract){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = round((dim (t_din)[1]) * fract)
  
  dtrain = t_din[v_sample[(i_limitc + 1):(length (v_sample))],]
  dtest = t_din[v_sample[1:i_limitc],]
  
  return (list(dtrain, dtest))
  
}


samplingDataNgroup = function (t_din, i_nb_group){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = as.integer (dim(t_din)[1] / i_nb_group)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = t_din[v_sample[1:i_limitc],]
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = t_din[v_sample[((i_limitc * (g-1) + 1):(length (v_sample)))],]
    }
    else{
      t_group = t_din[(v_sample[(i_limitc * (g-1) + 1):(i_limitc * g)]),]
    }
    # append list
    output[[g]] = t_group
  }
  
  return (output)
}


#################################
#   Control dataset integrity   #
#################################

controlDatasets = function(ldataset, prin){
  
  nbsplit = length(ldataset)
  
  pdf(paste(prin,".pdf", sep = ""), width = 10, height = 10)
  colorrainbow = rainbow(nbsplit)
  i = 1
  colorpoint = NULL
  dPCA = NULL
  dglobal = NULL
  for (d in ldataset){
    dglobal = rbind(dglobal,d)
    h = ggplot(d, aes(x=Aff)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
      labs(x = "pAff", y = "Frequencies") + 
      ggtitle(paste("Fold ", i, "-Dim=", dim(d)[1], sep = ""))+
      geom_density(alpha=.2, fill="#FAFFA5")
    print(h)
    
    #points for PCA
    colorpoint = append(colorpoint, rep(colorrainbow[i], (dim(d)[1])))
    dPCA = rbind(dPCA, d[,-dim(d)[2]]) # remove affinity coulum
    i = i + 1
  }
  
  h = ggplot(dglobal, aes(x=Aff)) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "pAff", y = "Frequencies") + 
    ggtitle(paste("Fold ", i, "-Dim=", dim(dglobal)[1], sep = "")) +
    geom_density(alpha=.2, fill="#FAFFA5")
  print (h)
  
  # plot PCA
  dplot = generatePCAcoords(dPCA)
  var_cap = dplot[[2]]
  data_plot = dplot[[1]]
  #par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = colorpoint, xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 1.5, cex.main = 2, cex.axis = 1.5, cex = 2)
  abline(h=0,v=0)
  dev.off()
  
}
