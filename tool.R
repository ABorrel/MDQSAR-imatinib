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


openData = function (pfilin, valcor, prout, vexclude){
	desc = read.csv (pfilin, header = TRUE, sep = "\t")
	#print("ddddd")
	#print (desc)
	#print(dim(desc))
	#print (rownames(desc))
	#rownames(desc) = desc[,1]
	#desc = desc[,-1]
	#print (desc)

	# deleted line with NA
	#rownames (desc) = seq (1, dim(desc)[1])
	desc = na.omit(desc)
  #print (dim(desc))
  cexclude = desc[,vexclude]
  #print(dim(desc))
  desc = desc[,-vexclude]
  
	# dell when sd = 0

	sd_desc =apply (desc[,1:(dim(desc)[2])], 2, sd)

	#print (sd_desc)
	#print ("--------")
	sd_0 = which (sd_desc == 0)

	#print (sd_0)

	#print ("------------")
	#print (mode(sd_0))
	#print (length (sd_0))
	#print ("------------")
	if (length(sd_0) != 0){
		#print (as.factor (sd_0))
		#desc = desc[,-sd_0]
		desc=subset(desc,select=-sd_0)
		cexclude = subset(cexclude,select=-sd_0)
		#print(dim(desc_new))
	}
	if (valcor != 0){
		out_elimcor = elimcor_sansY (desc, valcor)
		descriptor = out_elimcor$possetap

		MDSElimcor (desc, out_elimcor, paste (prout, "MDSDesc_", valcor, sep = ""), "corr")
		descriptor = colnames (desc) [descriptor]
		desc = desc[,descriptor]
		#print (dim(desc))
	}
  desc = cbind(cexclude, desc)
  #print(dim(desc))
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
