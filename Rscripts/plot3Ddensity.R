#!/usr/bin/env Rscript

library(plotly)
library(interp)
library(ggplot2)
library(RColorBrewer)

source("PCAplot.R")
source ("~/development/Rglobal/source/dataManager.R")


##########
#  MAIN  #
##########

args <- commandArgs(TRUE)
pdesc = args[1]
paff = args[2]
prout = args[3]

pdesc = "/home/borrela2/imatinib/results/analysis/densityMD/masterM_Ki"
#pdesc = "~/imatinib/results/analysis/QSARs/Lig2D_Ki/globalSet.csv"
paff = "/home/borrela2/imatinib/results/CHEMBL/AffAllcurated"
prout = "/home/borrela2/imatinib/results/analysis/densityMD/"

din = openData(pdesc, 0.9, prout, NbmaxNA = 10000)
din = din[[1]]

daff = read.csv(paff, sep = "\t")
rownames(daff) = daff[,1]

coordPCA = generatePCAcoords(din)

# create vector for all frame
vaff = NULL
for (chemIDframe in rownames(coordPCA[[1]])){
  chemID = strsplit(chemIDframe, "_")[[1]][1]
  for(chemIDaff in rownames(daff)){
    if(chemID == chemIDaff){
      vaff = append(vaff, daff[chemID,"Aff"])
      break()
    }
  }
  print (length(vaff))
}

f <- list(
  family = "Old Standard TT, serif",
  size = 34,
  color = "black"
)


X = coordPCA[[1]][,1]
Y = coordPCA[[1]][,2]

di <- interp(X, Y, vaff)#,
             #xo=seq(min(X), max(X), length=75),
             #yo=seq(min(Y), max(Y), length=75))



dat_interp <- data.frame(x=di$x, y=di$y, z=di$z)

#di <- data.frame(x=coordPCA[[1]][,1], y=coordPCA[[1]][,2], daff[,1])

p <- plot_ly(x = di$x, y = di$y, z=t(di$z), type = "surface", surfacecolor= t(di$z),
             colors=colorRampPalette(brewer.pal(11,"Spectral"))(500))


p <- p %>% layout(p,              # all of layout's properties: /r/reference/#layout
                  scene = list(
                    xaxis = list(title = "PC1", titlefont = f, range= list(-10.5,10), showticklabels = F),
                    yaxis = list(title = "PC2", titlefont = f, range = list(-4.3,6), showticklabels = F),
                    zaxis = list(title = "pKi", titlefont = f, showticklabels = F),
                    camera = list(eye = list(x = (-1.25/2), y = (-2.25/2), z = (3/2)))
                  )
)

