#!/usr/bin/env Rscript

# By Alexandre BORREL
# 02-2016

library (MASS)
require(plotrix)
require(lattice)
library(scatterplot3d)

source("elimcor_sansY.R")


#######################
#  GRAPHICS MANAGERS  #
#######################

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}



factorACP = function (coor_point, vector_arrows){

	factor = 1
	orgin = vector_arrows
	while (max (vector_arrows[,1]) < max (coor_point[,1]) && max (vector_arrows[,2]) < max (coor_point[,2]) && min (vector_arrows[,1]) > min (coor_point[,1]) && min (vector_arrows[,2]) > min (coor_point[,2]) ){
		factor = factor + 1
		vector_arrows[,1] = vector_arrows[,1] + orgin[,1]
		vector_arrows[,2] = vector_arrows[,2] + orgin[,2]

	}
	return (factor-1)

}

# MAYBE USEFULL
#colorACPByTypeOfDescriptors = function (){
#	d.col = read.csv ("temp_color", sep = "\t", header = TRUE)
#	d.col[which (d.col == "black")] = "grey"
	#print (as.matrix(d.col))
#	return (d.col)
#}


colorDesc = function(ldesc){
  mcolor = read.table("colorParameter", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  rownames(mcolor) = mcolor[,1]
  #mcolor = mcolor[,-1]
  lout = mcolor[ldesc,2]
  return(lout)
}


calculLimListData = function (data1, data2){

	x1 = min (data1)
	x2 = min (data2)
	X1 = max(data1)
	X2 = max (data2)
	l = c (min (x1,x2), max (X1,X2))
	return (l)
}


MDSElimcor = function (din, out_elimcor, path_file, type_dist){

	groupe_elimcor = out_elimcor$groupes
	descriptor_selected = out_elimcor$possetap

	#print(descriptor_selected)
	name_descriptor = colnames(din)[descriptor_selected]
	din = din[,descriptor_selected]
	
	din = na.omit (din)
	if (type_dist == "corr"){
	  MC = cor(din)
	  dist1 = abs(1-MC)  
	}

	color_desc = "black" #colorDesc(name_descriptor)
	#print (color_desc)

	png (paste (path_file, "_", type_dist, ".png", sep = ""), 3500, 3500, res = 180, bg = "white")
	par( mar=c(6,6,6,6))
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	#c = as.vector (col.desc[,colnames (data)])
	plot (fit$points[,1], fit$points[,2], main="MDS Descriptor", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1, 1))
	#text (fit$points[,1], fit$points[,2]+0.05, labels = num_desc,  cex = 1.8, col = color_desc, font = 12)
	#text (fit$points[descriptor_selected,1], fit$points[descriptor_selected,2]+0.02, labels = "*",  cex = 4, col = color_desc)
	text (fit$points[,1], fit$points[,2], labels = name_descriptor,  cex = 2.6, col = color_desc)

	dev.off()
}



# distribution
histData = function (data1, data2, path_pdf){

	nb_descriptor = dim (data1)[2]
	pdf (path_pdf)
	for (i_descriptor in 1:nb_descriptor){
		l = list (data1[,i_descriptor], data2[,i_descriptor])
		multhist(l, col = c(2,1), cex.names = 1.5, freq = TRUE, cex.axis = 1.5)
		title(main=colnames (data1)[i_descriptor], cex.main = 1.5, ylab = "Frequency", cex.lab = 1.5, sub=paste("Means:", round(mean(data1[,i_descriptor]),3),round(mean(data2[,i_descriptor]),3),"// SD:", round(sd(data1[,i_descriptor]),3), round(sd(data2[,i_descriptor]),3) ," // Med:",round(median(data1[,i_descriptor]),3), round(median(data2[,i_descriptor]),3)))
		legend("topright", col=c(2,1), legend = c("Class1", "Class2"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 1.5)
	}
	dev.off()
}

histDataOne = function (data1, path_pdf){

	nb_descriptor = dim (data1)[2]
	pdf (path_pdf)
	for (i_descriptor in 1:nb_descriptor){
		hist(data1[,i_descriptor], col = "red", cex.names = 1.5, freq = TRUE, cex.axis = 1.5, main = "", xlab = "", ylab = "", breaks = 100)
		title(main = colnames (data1)[i_descriptor], cex.main = 1.5, ylab = "Frequency", cex.lab = 1.5, sub=paste("Means:", round(mean(data1[,i_descriptor]),3), "// SD:", round(sd(data1[,i_descriptor]),3), " // Med:",round(median(data1[,i_descriptor]),3), " // Min:", round(min(data1[,i_descriptor]),3), " // Max:", round(max(data1[,i_descriptor]),3)))
	}
	dev.off()
}

histOneData = function (d, p_filin){
	nb_descriptor = dim (d)[2]
	for (i_descriptor in 1:nb_descriptor){
		svg (paste (p_filin, colnames (d)[i_descriptor],".svg", sep = ""), 20, 15)
		hist(d[,i_descriptor], col = "grey", cex.names = 1.5,cex.main = 6, freq = TRUE, cex.axis = 3.5, breaks = 10, main = colnames (d)[i_descriptor], ylab = "", xlab = "")
		dev.off()
	}

}


boxplotData = function (data1,data2, p_filout, format){
	nb_descriptor = dim (data1)[2]
	if (format == "pdf"){
		pdf (p_filout)
		for (i_descriptor in 1:nb_descriptor){
			l = list (data1[,i_descriptor], data2[,i_descriptor])
			boxplot(l, col = c(2,8), cex.names = 1.5, freq = TRUE, cex.axis = 1.5)
			title(main=colnames (data1)[i_descriptor], cex.main = 1.5, ylab = "Frequency", cex.lab = 1.5, sub=paste("Means:", round(mean(data1[,i_descriptor]),3),round(mean(data2[,i_descriptor]),3),"// SD:", round(sd(data1[,i_descriptor]),3), round(sd(data2[,i_descriptor]),3) ," // Med:",round(median(data1[,i_descriptor]),3), round(median(data2[,i_descriptor]),3)))
			#legend("topright", col=c(2,1), legend = c("Class1", "Class2"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 1.5)
		}
	dev.off ()
	}
	if (format == "svg"){
		for (i_descriptor in 1:nb_descriptor){
			svg (paste (p_filout, colnames (data1)[i_descriptor],".svg", sep = ""), 7, 20)
			par(mar=c(8,8,8,8))
			l = list (data1[,i_descriptor], data2[,i_descriptor])
			boxplot(l, col = c(2,8), cex.names = 2.5, freq = TRUE, cex.axis = 1.5)
			title(main="", cex.main = 1.5, ylab = colnames (data1)[i_descriptor], cex.lab = 4.5,cex.axis = 2, sub=paste("Means:", round(mean(data1[,i_descriptor]),3),round(mean(data2[,i_descriptor]),3),"// SD:", round(sd(data1[,i_descriptor]),3), round(sd(data2[,i_descriptor]),3) ," // Med:",round(median(data1[,i_descriptor]),3), round(median(data2[,i_descriptor]),3)))
			#legend("topright", col=c(2,1), legend = c("Class1", "Class2"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 1.5)
			dev.off ()
		}
	}
}




# check color vector
CheckColorVector = function (l_descriptor, col.des){
	for (d in l_descriptor){
		if (is.integer0(which(d == names(col.des))) == TRUE){
			out = rep (1,length (l_descriptor))
			names (out) = l_descriptor
			return (out)
		}
	}
	return (col.des)
}



gifGeneration = function(pin, d){
  
  png(file=paste(pin, "3DPlot%03d.png", sep = ""), width=800, height=800)
  for (i in seq(0, 350 ,10)){
    print(cloud(d[,1]~d[,2]*d[,3], col = "black", pch = 16, screen = list(z = i, x = -80)))
  }
  dev.off()
  
  # convert pngs to one gif using ImageMagick
  system(paste("convert -delay 50 ", pin, "*.png ", pin, ".gif", sep = ""))
  
}






####################
# COMPARISON TEST  #
####################

# choice paramtric or not with bartlet and shapiro

conditionANOVA = function (d_in, class_name){
	# 1 -> parametric // 0 -> no paramtric

	l_value_class = unique (d_in[,class_name])
	for (val_class in l_value_class){

		d_temp = d_in[which (d_in[,class_name] == val_class), ]
		d_temp = d_in [,-(which (colnames (d_in) == class_name))]

		pval = shapiro.test (d_temp)$p.value
		# H0 -> normalite
		if (pval < 0.005 ){
			return (0)
		}
	}

	p_val_bartlett = bartlett.test(d_in[,1]~d_in[,which (colnames (d_in) == class_name)])$p.value
	if (p_val_bartlett < 0.005){
		return (0)

	}
	return (1)
}


# choice parametric or non parametric test comparison



conditionTtest = function (data1, data2){

	sd1 = sd (data1)
	sd2 = sd (data2)

	if (is.na(sd1)){
		sd1 = 0
	}
	if (is.na(sd2)){
		sd2 = 0
	}

	if(sd1 == 0 || sd2 == 0){
		return (2)
	}

	if (length (data1) < 30 || length (data2) < 30 ){
		#normalisation
		pval1 = shapiro.test (data1)$p.value
		pval2 = shapiro.test (data2)$p.value

		if (pval1 < 0.005 || pval2 < 0.005){
			return (0)
		}
	}
	pval_var = var.test (data1, data2)$p.value
	if (pval_var < 0.05){
		return (0)
	}
	return (1)
}

# run test comparison
comparisonTest = function (vector_value1, vector_value2, type){

	if (type == "parametric"){
		result = t.test (vector_value1, vector_value2)

	}else{
		result = wilcox.test (vector_value1, vector_value2)

	}
	return (result$p.value)
}


# test ANOVA
ANOVA = function (d_in, class_name){


	condition = conditionANOVA (d_in, class_name)
	print (condition)
	# 1 -> parametric // 2 -> non parametric

	if (condition == 1){
		r_anova = aov(d_in[,1]~d_in[,which (colnames (d_in) == class_name)])
		p_val = summary (r_anova)[[1]][,5][1]
		signif = signifPvalue (p_val)

		l_pval = c (format(p_val, scientific = TRUE, digit = 2), signif)
		return (l_pval)
	}
	else {

		r_krutal_wallis_test = kruskal.test(d_in[,1]~d_in[,which (colnames (d_in) == class_name)])
		p_val =r_krutal_wallis_test$p.value

		signif = signifPvalue (p_val)
		l_pval = c (format(p_val, scientific = TRUE, digit = 2), signif)
		return (l_pval)

	}

}



###########################################
# retrieve more significative descriptors #
###########################################


moreSignif = function (data1, data2){
	p_val = 100
	des_out = NULL
	nb_descriptor = dim (data2)[2]
	for (i_des in seq (1,nb_descriptor)){
		if (conditionTtest(data1[,i_des], data2[,i_des]) == 1){
			pval_temp = comparisonTest(data1[,i_des], data2[,i_des], "parametric")
		}else if (conditionTtest(data1[,i_des], data2[,i_des]) == 0){
			pval_temp = comparisonTest(data1[,i_des], data2[,i_des], "non-parametric")
		}else{
			pval_temp = 100
		}
		if (pval_temp < p_val){
			des_out = i_des
			p_val = pval_temp
		}
	}
	return (colnames (data1)[des_out])
}



#################
# DATA MANAGERS #
#################

# separate data with class group
separeData = function (data, descriptor_class){

	data1 = data [which(data[,descriptor_class] == 0),]
	data2 = data [which(data[,descriptor_class] == 1),]

	return (list (data1, data2))
}


openData = function (pfilin, valcor, prout){
	desc = read.csv (pfilin, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
	#print("ddddd")
	#print (desc)
	#print(dim(desc))
	#print (rownames(desc))
	rownames(desc) = desc[,1]
	
	desc = desc[,-1]
	#print (desc)

	# deleted col with NA
	lcoldel = NULL
	print (dim(desc))
	for (icol in seq(1, dim(desc)[2])){
	  #print (desc[,icol])
	  if (sum(is.na(as.vector(desc[,icol]))) > 10){
	    lcoldel = append(lcoldel, icol)
	  }
	}
	
	if(is.null(lcoldel)){
	  desc = na.omit(desc)
	}else{
	  desc = desc[,-lcoldel]	  
	  desc = na.omit(desc)  
	}
	
	# dell when sd = 0
	sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)

	#print (sd_desc)
	#print ("--------")
	sd_0 = which (sd_desc == 0.0)

	#print ("------------")
	#print (mode(sd_0))
	#print (length (sd_0))
	#print ("------------")
	if (length(sd_0) != 0){
		#print (as.factor (sd_0))
		#desc = desc[,-sd_0]
		desc=subset(desc,select=-sd_0)
	}
	if (valcor != 0){
	  out_elimcor = elimcor_sansY (desc, valcor)
		descriptor = out_elimcor$possetap

		MDSElimcor (desc, out_elimcor, paste (prout, "MDSDesc_", valcor, sep = ""), "corr")
		descriptor = colnames (desc) [descriptor]
		desc = desc[,descriptor]
		#print (dim(desc))
	}
	
	# again with SD null
	sd_desc =apply (desc[,1:(dim(desc)[2])], 2, sd)
	sd_0 = which (sd_desc == 0)
	if (length(sd_0) != 0){
	  desc=subset(desc,select=-sd_0)
	}
	return (list((desc),colnames (desc)))
}

delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


# change 0 or 1 by d and nd
changeList = function (list_element){
	Y = as.factor(list_element)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	return (Y2)
}


changeProbaList = function (list_element){
	Y = list_element
	Y2 = rep(1,length(Y))
	Y2[which(Y<0.5)]=0
	return (Y2)
}

# del colum
delCol = function (data, list_name_col){
	for (name_col in list_name_col){
		data = data[,-which(colnames(data)==name_col)]
	}
	return (data)
}


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


delnohomogeniousdistribution = function(din, cutoff = 80){

  countMax = dim(din)[1]*cutoff/100  
  
  i = 1
  imax = dim(din)[2]
  while(i <= imax){
    #print (i)
    #print (din[,i])
    qt = hist(din[,i], breaks = 10, plot = FALSE)$counts
      
    for (qtc in qt){
      if (qtc >= countMax){
        din = din[,-i]
        imax = imax - 1
        i = i - 1
        break()
      }
    }
    i = i + 1
  }
  return(din)
}

###############################
# divise the dataset in folds #
###############################

samplingDataNgroup = function (t_din, i_nb_group, s_nameclass){

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


#######################
# order data rownames #
#######################


orderByType = function (d.col, d){
	l_temp = NULL

	for (desc in colnames (d.col)){
		if (is.integer0 (which(rownames (d) == desc))== FALSE){
			l_temp = append (l_temp, desc)
		}
	}
	if (length (l_temp) != length (rownames (d))){
		print ("TO DO -> return list desc")
		return (l_temp)
	}
	d = d[l_temp,]
	#names (d) = l_temp
	return (d)
}


orderByTypeOneCol = function (d.col, d){
	l_temp = NULL

	for (desc in colnames (d.col)){
		if (is.integer0 (which(rownames (d) == desc))== FALSE){
			l_temp = append (l_temp, desc)
		}
	}
	if (length (l_temp) != length (rownames (d))){
		print ("TO DO -> return list desc")
		return (l_temp)
	}
	d = d[l_temp,]
	names (d) = l_temp
	return (d)
}



listDescriptorShipshape = function (d.col, l_des){
	l_temp = NULL

	for (desc in names (d.col)){
		if (desc %in% l_des == TRUE){
			l_temp = append (l_temp, desc)
		}
	}
	return (l_temp)
}


######################
# Normalization LDA  #
######################


normalizationCoef = function (coef, data_train){

	d_class1 = data_train[which(data_train[,"drugg"]==0),]
	d_class2 = data_train[which(data_train[,"drugg"]==1),]

	m_class1 = mean (d_class1[,1])
	m_class2 = mean (d_class2[,1])

	v_c1 =  sum((d_class1[,1]-m_class1) * (d_class1[,1]-m_class1))
	v_c2 =  sum((d_class2[,1]-m_class2) * (d_class2[,1]-m_class2))


	v_out = sqrt((v_c1 + v_c2) / (dim(data_train)[1] - 2))

	return (coef*v_out)

}


normalizationScalingLDA = function (scalingLDA, d){

	l_out = NULL
	l_des = names (scalingLDA)

	for (desc in l_des){
		l_out = append (l_out, normalizationCoef (scalingLDA[desc], d[,c(desc,"drugg")]))
	}
	names (l_out) = names (scalingLDA)
	return (l_out)
}



###############
#   Signif    #
###############



signifPvalue = function (a){
	if (a < 0.001){
		return ("***")
	}
	else if (a < 0.01){
		return ("**")
	}
	else if (a < 0.05){
		return ("*")
	}
	else {
		return ("-")
	}
}



##############
# plot tree  #
##############

plot.nnet<-function(mod.in,nid=T,all.out=T,all.in=T,bias=T,wts.only=F,rel.rsc=5,
                    circle.cex=5,node.labs=T,var.labs=T,x.lab=NULL,y.lab=NULL,
                    line.stag=NULL,struct=NULL,cex.val=1,alpha.val=1,
                    circle.col='lightblue',pos.col='black',neg.col='grey',
                    bord.col='black', ...){
  
  require(scales)
  
  #sanity checks
  if('mlp' %in% class(mod.in)) warning('Bias layer not applicable for rsnns object')
  if('numeric' %in% class(mod.in)){
    if(is.null(struct)) stop('Three-element vector required for struct')
    if(length(mod.in) != ((struct[1]*struct[2]+struct[2]*struct[3])+(struct[3]+struct[2])))
      stop('Incorrect length of weight matrix for given network structure')
  }
  if('train' %in% class(mod.in)){
    if('nnet' %in% class(mod.in$finalModel)){
      mod.in<-mod.in$finalModel
      warning('Using best nnet model from train output')
    }
    else stop('Only nnet method can be used with train object')
  }
  
  #gets weights for neural network, output is list
  #if rescaled argument is true, weights are returned but rescaled based on abs value
  nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct){
    
    require(scales)
    require(reshape)
    
    if('numeric' %in% class(mod.in)){
      struct.out<-struct
      wts<-mod.in
    }
    
    #neuralnet package
    if('nn' %in% class(mod.in)){
      struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
      struct.out<-struct.out[-length(struct.out)]
      struct.out<-c(
        length(mod.in$model.list$variables),
        struct.out,
        length(mod.in$model.list$response)
      )    		
      wts<-unlist(mod.in$weights[[1]])   
    }
    
    #nnet package
    if('nnet' %in% class(mod.in)){
      struct.out<-mod.in$n
      wts<-mod.in$wts
    }
    
    #RSNNS package
    if('mlp' %in% class(mod.in)){
      struct.out<-c(mod.in$nInputs,mod.in$archParams$size,mod.in$nOutputs)
      hid.num<-length(struct.out)-2
      wts<-mod.in$snnsObject$getCompleteWeightMatrix()
      
      #get all input-hidden and hidden-hidden wts
      inps<-wts[grep('Input',row.names(wts)),grep('Hidden_2',colnames(wts)),drop=F]
      inps<-melt(rbind(rep(NA,ncol(inps)),inps))$value
      uni.hids<-paste0('Hidden_',1+seq(1,hid.num))
      for(i in 1:length(uni.hids)){
        if(is.na(uni.hids[i+1])) break
        tmp<-wts[grep(uni.hids[i],rownames(wts)),grep(uni.hids[i+1],colnames(wts)),drop=F]
        inps<-c(inps,melt(rbind(rep(NA,ncol(tmp)),tmp))$value)
      }
      
      #get connections from last hidden to output layers
      outs<-wts[grep(paste0('Hidden_',hid.num+1),row.names(wts)),grep('Output',colnames(wts)),drop=F]
      outs<-rbind(rep(NA,ncol(outs)),outs)
      
      #weight vector for all
      wts<-c(inps,melt(outs)$value)
      assign('bias',F,envir=environment(nnet.vals))
    }
    
    if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))
    
    #convert wts to list with appropriate names 
    hid.struct<-struct.out[-c(length(struct.out))]
    row.nms<-NULL
    for(i in 1:length(hid.struct)){
      if(is.na(hid.struct[i+1])) break
      row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
    }
    row.nms<-c(
      row.nms,
      rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
    )
    out.ls<-data.frame(wts,row.nms)
    out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
    out.ls<-split(out.ls$wts,f=out.ls$row.nms)
    
    assign('struct',struct.out,envir=environment(nnet.vals))
    
    out.ls
    
  }
  
  wts<-nnet.vals(mod.in,nid=F)
  
  if(wts.only) return(wts)
  
  #circle colors for input, if desired, must be two-vector list, first vector is for input layer
  if(is.list(circle.col)){
    circle.col.inp<-circle.col[[1]]
    circle.col<-circle.col[[2]]
  }
  else circle.col.inp<-circle.col
  
  #initiate plotting
  x.range<-c(0,100)
  y.range<-c(0,100)
  #these are all proportions from 0-1
  if(is.null(line.stag)) line.stag<-0.011*circle.cex/2
  layer.x<-seq(0.17,0.9,length=length(struct))
  bias.x<-layer.x[-length(layer.x)]+diff(layer.x)/2
  bias.y<-0.95
  circle.cex<-circle.cex
  
  #get variable names from mod.in object
  #change to user input if supplied
  print (attributes(mod.in$nnet))
  
  if('numeric' %in% class(mod.in)){
    x.names<-paste0(rep('X',struct[1]),seq(1:struct[1]))
    y.names<-paste0(rep('Y',struct[3]),seq(1:struct[3]))
  }
  if('mlp' %in% class(mod.in)){
    all.names<-mod.in$snnsObject$getUnitDefinitions()
    x.names<-all.names[grep('Input',all.names$unitName),'unitName']
    y.names<-all.names[grep('Output',all.names$unitName),'unitName']
  }
  if('nn' %in% class(mod.in)){
    x.names<-mod.in$model.list$variables
    y.names<-mod.in$model.list$respons
  }
  if('xNames' %in% names(mod.in)){
    x.names<-mod.in$xNames
    y.names<-attr(terms(mod.in),'factor')
    y.names<-row.names(y.names)[!row.names(y.names) %in% x.names]
  }
  if(!'xNames' %in% names(mod.in) & 'nnet' %in% class(mod.in)){
    if(is.null(mod.in$call$formula)){
      x.names<-colnames(eval(mod.in$call$x))
      y.names<-colnames(eval(mod.in$call$y))
    }
    else{
      forms<-eval(mod.in$call$formula)
      x.names<-mod.in$coefnames
      facts<-attr(terms(mod.in),'factors')
      y.check<-mod.in$fitted
      if(ncol(y.check)>1) y.names<-colnames(y.check)
      else y.names<-as.character(forms)[2]
    } 
  }
  #change variables names to user sub 
  if(!is.null(x.lab)){
    if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
    else x.names<-x.lab
  }
  if(!is.null(y.lab)){
    if(length(y.names) != length(y.lab)) stop('y.lab length not equal to number of output variables')
    else y.names<-y.lab
  }
  
  #initiate plot
  plot(x.range,y.range,type='n',axes=F,ylab='',xlab='',...)
  
  #function for getting y locations for input, hidden, output layers
  #input is integer value from 'struct'
  get.ys<-function(lyr){
    spacing<-diff(c(0*diff(y.range),0.9*diff(y.range)))/max(struct)
    seq(0.5*(diff(y.range)+spacing*(lyr-1)),0.5*(diff(y.range)-spacing*(lyr-1)),
        length=lyr)
  }
  
  #function for plotting nodes
  #'layer' specifies which layer, integer from 'struct'
  #'x.loc' indicates x location for layer, integer from 'layer.x'
  #'layer.name' is string indicating text to put in node
  layer.points<-function(layer,x.loc,layer.name,cex=cex.val){
    x<-rep(x.loc*diff(x.range),layer)
    y<-get.ys(layer)
    points(x,y,pch=21,cex=circle.cex,col=bord.col,bg=in.col)
    if(node.labs) text(x,y,paste(layer.name,1:layer,sep=''),cex=cex.val)
    if(layer.name=='I' & var.labs) text(x-line.stag*diff(x.range),y,x.names,pos=2,cex=cex.val)      
    if(layer.name=='O' & var.labs) text(x+line.stag*diff(x.range),y,y.names,pos=4,cex=cex.val)
  }
  
  #function for plotting bias points
  #'bias.x' is vector of values for x locations
  #'bias.y' is vector for y location
  #'layer.name' is  string indicating text to put in node
  bias.points<-function(bias.x,bias.y,layer.name,cex,...){
    for(val in 1:length(bias.x)){
      points(
        diff(x.range)*bias.x[val],
        bias.y*diff(y.range),
        pch=21,col=bord.col,bg=in.col,cex=circle.cex
      )
      if(node.labs)
        text(
          diff(x.range)*bias.x[val],
          bias.y*diff(y.range),
          paste(layer.name,val,sep=''),
          cex=cex.val
        )
    }
  }
  
  #function creates lines colored by direction and width as proportion of magnitude
  #use 'all.in' argument if you want to plot connection lines for only a single input node
  layer.lines<-function(mod.in,h.layer,layer1=1,layer2=2,out.layer=F,nid,rel.rsc,all.in,pos.col,
                        neg.col,...){
    
    x0<-rep(layer.x[layer1]*diff(x.range)+line.stag*diff(x.range),struct[layer1])
    x1<-rep(layer.x[layer2]*diff(x.range)-line.stag*diff(x.range),struct[layer1])
    
    if(out.layer==T){
      
      y0<-get.ys(struct[layer1])
      y1<-rep(get.ys(struct[layer2])[h.layer],struct[layer1])
      src.str<-paste('out',h.layer)
      
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-wts[grep(src.str,names(wts))][[1]][-1]
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-wts.rs[grep(src.str,names(wts.rs))][[1]][-1]
      
      cols<-rep(pos.col,struct[layer1])
      cols[wts<0]<-neg.col
      
      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)
      
    }
    
    else{
      
      if(is.logical(all.in)) all.in<-h.layer
      else all.in<-which(x.names==all.in)
      
      y0<-rep(get.ys(struct[layer1])[all.in],struct[2])
      y1<-get.ys(struct[layer2])
      src.str<-paste('hidden',layer1)
      
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-unlist(lapply(wts[grep(src.str,names(wts))],function(x) x[all.in+1]))
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-unlist(lapply(wts.rs[grep(src.str,names(wts.rs))],function(x) x[all.in+1]))
      
      cols<-rep(pos.col,struct[layer2])
      cols[wts<0]<-neg.col
      
      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)
      
    }
    
  }
  
  bias.lines<-function(bias.x,mod.in,nid,rel.rsc,all.out,pos.col,neg.col,...){
    
    if(is.logical(all.out)) all.out<-1:struct[length(struct)]
    else all.out<-which(y.names==all.out)
    
    for(val in 1:length(bias.x)){
      
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      
      if(val != length(bias.x)){
        wts<-wts[grep('out',names(wts),invert=T)]
        wts.rs<-wts.rs[grep('out',names(wts.rs),invert=T)]
        sel.val<-grep(val,substr(names(wts.rs),8,8))
        wts<-wts[sel.val]
        wts.rs<-wts.rs[sel.val]
      }
      
      else{
        wts<-wts[grep('out',names(wts))]
        wts.rs<-wts.rs[grep('out',names(wts.rs))]
      }
      
      cols<-rep(pos.col,length(wts))
      cols[unlist(lapply(wts,function(x) x[1]))<0]<-neg.col
      wts.rs<-unlist(lapply(wts.rs,function(x) x[1]))
      
      if(nid==F){
        wts.rs<-rep(1,struct[val+1])
        cols<-rep('black',struct[val+1])
      }
      
      if(val != length(bias.x)){
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1]),
          lwd=wts.rs,
          col=cols
        )
      }
      
      else{
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1])[all.out],
          lwd=wts.rs[all.out],
          col=cols[all.out]
        )
      }
      
    }
  }
  
  #use functions to plot connections between layers
  #bias lines
  if(bias) bias.lines(bias.x,mod.in,nid=nid,rel.rsc=rel.rsc,all.out=all.out,pos.col=alpha(pos.col,alpha.val),
                      neg.col=alpha(neg.col,alpha.val))
  
  #layer lines, makes use of arguments to plot all or for individual layers
  #starts with input-hidden
  #uses 'all.in' argument to plot connection lines for all input nodes or a single node
  if(is.logical(all.in)){  
    mapply(
      function(x) layer.lines(mod.in,x,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[1]
    )
  }
  else{
    node.in<-which(x.names==all.in)
    layer.lines(mod.in,node.in,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,all.in=all.in,
                pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
  }
  #connections between hidden layers
  lays<-split(c(1,rep(2:(length(struct)-1),each=2),length(struct)),
              f=rep(1:(length(struct)-1),each=2))
  lays<-lays[-c(1,(length(struct)-1))]
  for(lay in lays){
    for(node in 1:struct[lay[1]]){
      layer.lines(mod.in,node,layer1=lay[1],layer2=lay[2],nid=nid,rel.rsc=rel.rsc,all.in=T,
                  pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
    }
  }
  #lines for hidden-output
  #uses 'all.out' argument to plot connection lines for all output nodes or a single node
  if(is.logical(all.out))
    mapply(
      function(x) layer.lines(mod.in,x,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[length(struct)]
    )
  else{
    node.in<-which(y.names==all.out)
    layer.lines(mod.in,node.in,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                pos.col=pos.col,neg.col=neg.col,all.out=all.out)
  }
  
  #use functions to plot nodes
  for(i in 1:length(struct)){
    in.col<-circle.col
    layer.name<-'H'
    if(i==1) { layer.name<-'I'; in.col<-circle.col.inp}
    if(i==length(struct)) layer.name<-'O'
    layer.points(struct[i],layer.x[i],layer.name)
  }
  
  if(bias) bias.points(bias.x,bias.y,'B')
  
}

