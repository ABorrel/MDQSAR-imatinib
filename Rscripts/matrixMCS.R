#!/usr/bin/env Rscript


require(graphics)

generateLegend = function (length_tot, cut){

  list_out = seq(0,length_tot, cut)
  if (list_out[length(list_out)] != length_tot){
    list_out = append(list_out, length_tot)
  }
  return (list_out)
}

generatePosition = function (list_legend, value_ecart){

  list_out = c()
  for (element in list_legend){
    list_out = append(list_out,element * value_ecart)
  }
  return (list_out)
}


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


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






cardAffinityText = function(matrixIN, daff,dtext ,name_file){

  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]

  dim_x = nb_line
  dim_y = nb_col

  if (nb_col == 1){
    return ()
  }
  if (nb_line == 1){
    return ()
  }

  if (nb_col < 30){
    dim_x = 30
  }

  if (nb_line < 30){
    dim_y = 30
  }

  bk = c(0,0.20,0.40,0.60,0.80,1)

  png (file = paste (name_file, ".png", sep = ""), dim_x * 55, dim_y * 55)
  par(mar=c(20,20,3,10))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", breaks = bk, col = c("#FFFFFF", "#FFBFBF","#FF8080", "#FF4040", "#FF0000"))
  grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  box()
  # place les petites barres
  axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)

  axis(4,seq(0,1,(1/(nb_col-1))), labels = FALSE)

  # place les positions en fonction du cut
  ecart1 = 1/(nb_line-1)
  ecart2 = 1/(nb_col-1)
  list_L1 = generateLegend (nb_line,1)
  list_L2 = generateLegend (nb_col,1)

  nbcol = dim(matrixIN)[2]
  nbline = dim(matrixIN)[1]
  for (i in seq(0,nbline-1)){
    for (j in seq(0, nbcol-1)){
      text((1/(nbline-1))*i,(1/(nbcol-1))*j, labels = round(dtext[i+1,j+1], digits = 1), cex = 1.5)
    }
  }

  # place les legendes
  posX = generatePosition(list_L1, ecart1)
  posY = generatePosition(list_L2, ecart2)
  axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 2, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),rownames (matrixIN), cex.axis = 2, las = 2)

  print(daff)
  axis(4,seq(0,1,(1/(nb_col-1))),daff[rownames(matrixIN),2], cex.axis = 2, las = 2)


  #legend ("right", fill = c("black", "darkred", "red", "darkmagenta", "darkorchid","deepskyblue", "cyan", "white"), legend = c("0-2", "2-3", "3-4", "4-5", "5-6","6-7", "7-8", "> 10"), bg = addTrans ("#FFFFFF", 120) )
  dev.off()
}

cardAffinityTextVar = function(matrixIN, daff, dtext ,name_file){
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  dim_x = nb_line
  dim_y = nb_col
  
  if (nb_col == 1){
    return ()
  }
  if (nb_line == 1){
    return ()
  }
  
  if (nb_col < 30){
    dim_x = 30
  }
  
  if (nb_line < 30){
    dim_y = 30
  }
  
  bk = seq(min(matrixIN), max(matrixIN), length.out = 6)
  
  png (file = paste (name_file, ".png", sep = ""), dim_x * 55, dim_y * 55)
  par(mar=c(20,20,3,10))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", breaks = bk, col = c("#FFFFFF", "#FFBFBF","#FF8080", "#FF4040", "#FF0000"))
  grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  box()
  # place les petites barres
  axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  axis(4,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  # place les positions en fonction du cut
  ecart1 = 1/(nb_line-1)
  ecart2 = 1/(nb_col-1)
  list_L1 = generateLegend (nb_line,1)
  list_L2 = generateLegend (nb_col,1)
  
  nbcol = dim(matrixIN)[2]
  nbline = dim(matrixIN)[1]
  for (i in seq(0,nbline-1)){
    for (j in seq(0, nbcol-1)){
      text((1/(nbline-1))*i,(1/(nbcol-1))*j, labels = dtext[i+1,j+1], cex = 1.5)
    }
  }
  
  # place les legendes
  posX = generatePosition(list_L1, ecart1)
  posY = generatePosition(list_L2, ecart2)
  axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 2, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),rownames (matrixIN), cex.axis = 2, las = 2)
  
  axis(4,seq(0,1,(1/(nb_col-1))),daff[rownames(matrixIN),2], cex.axis = 2, las = 2)
  #legend ("right", fill = c("black", "darkred", "red", "darkmagenta", "darkorchid","deepskyblue", "cyan", "white"), legend = c("0-2", "2-3", "3-4", "4-5", "5-6","6-7", "7-8", "> 10"), bg = addTrans ("#FFFFFF", 120) )
  dev.off()
}



###########
#  MAIN   #
###########

args = commandArgs(TRUE)
pmatrix = args[1]
paffinity = args[2]
ptext = args[3]
probMatrix = as.integer(args[4])

# for test

#pmatrix = "/home/aborrel/imitanib/results/analysis/MCS/CHEMBL941_Tanimoto"
#paffinity = "/home/aborrel/imitanib/results/analysis/MCS/CHEMBL941_Aff"
#ptext = "/home/aborrel/imitanib/results/analysis/MCS/CHEMBL941_MAX"

#paffinity = "C://Users/Alexandre\ Borrel/Desktop/LSR/OT-55_1PTW/affinity"
#pmatrix = "C://Users/Alexandre\ Borrel/Desktop/LSR/OT-55_1PTW/matriceMCSTanimoto"
#ptext = "C://Users/Alexandre\ Borrel/Desktop/LSR/OT-55_1PTW/matriceMCSNbAtomDiff"
#pLSR = "C://Users/Alexandre\ Borrel/Desktop/LSR/OT-55_1PTW/listLSRsmiles"

d = read.csv (pmatrix, header = T, sep = "\t")
daff = read.table(paffinity, header = T, sep = "\t")
rownames(daff) = daff[,1]
dtext = read.table(ptext, header = T, sep = "\t")


if (probMatrix == 1){
  cardAffinityText(d, daff, dtext, pmatrix)
}else if (probMatrix == 0){
  cardAffinityTextVar(d, daff, dtext, pmatrix)
}else if (probMatrix == 3){
  # change color with difference of affinity
  ddiff = data.frame()
  for(i in seq(1,dim(daff)[1])){
    for (j in seq(1,dim(daff)[1])){
      ddiff[i,j] = abs(daff[i,2] - daff[j,2])
    }
  }
  
  colnames(ddiff) = colnames(d)
  rownames(ddiff) = rownames(d)
  
  # card with texte
  cardAffinityTextVar(ddiff, daff, dtext, pmatrix)
}