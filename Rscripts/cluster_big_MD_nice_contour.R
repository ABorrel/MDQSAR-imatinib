list.of.packages <- c("magrittr","plotly","gplot","akima","plyr",
                      "rgeos", "RColorBrewer","cluster","fpc","mclust","fastcluster","ggbiplot",
                      "boot","bootstrap","vegan","caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# you will need your own plotly account to get this to work

Sys.setenv("plotly_username" = "jrash")
Sys.setenv("plotly_api_key" = "XXXXXXX")
library(plotly)


rm(list=ls())
setwd("C:/Users/Vestige/Google Drive/FourchesLab/3I60/big_MD/")


# Load data ---------------------

load("clusterBigMDmore.RDATA")

#-----experiment with heatmap

# Make the contour data

PC1 = matTransYeo.pca$x[,1]
PC2 = matTransYeo.pca$x[,2]


di <- interp(PC1, PC2, mat_pki$pki,
             xo=seq(min(PC1), max(PC1), length=100),
             yo=seq(min(PC2), max(PC2), length=100))


dat_interp <- data.frame(x=di$x, y=di$y, z=di$z)

library(extrafont)
font_import(c("courier"))

# set up fonts

f <- list(
  family = "Old Standard TT, serif",
  size = 34,
  color = "black"
)
f2 <- list(
  family = "arialbd",
  size = 18,
  color = "grey"
)
f3 <- list(
  family = "Old Standard TT, serif",
  size = 24,
  color = "grey"
)
f4 <- list(
  family = "Old Standard TT, serif",
  size = 24,
  color = "black"
)

a <- list(xref= 'paper',
          yref= 'paper',
          x= 0,
          xanchor= 'right',
          y= 1,
          yanchor= 'bottom',
          text= 'X axis label',
          showarrow= F
)


# use pki on z-axis

p <- plot_ly(x = di$x, y = di$y, z=t(di$z), type = "surface", surfacecolor= t(di$z),
             colors=colorRampPalette(brewer.pal(11,"Spectral"))(100), colorbar = list(title = "pKi", titlefont = f4, tickfont=f4))

p <- p %>% layout(p,              # all of layout's properties: /r/reference/#layout
                  scene = list(
                    xaxis = list(title = "PC1", titlefont = f, range= list(-10.5,10), showticklabels = F),
                    yaxis = list(title = "PC2", titlefont = f, range = list(-4.3,6), showticklabels = F),
                    zaxis = list(title = "pKi", titlefont = f, showticklabels = F),
                    camera = list(eye = list(x = (-1.25/2), y = (-2.25/2), z = (3/2)))
                  )
)


print(p)

#export image
plotly_IMAGE(p, format = "png", out_file = "contour_3d_spectral_view1.png", width =1100, height = 1000)
