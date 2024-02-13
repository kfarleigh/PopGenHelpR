## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Load the package
library(PopGenHelpR)

## ----load data----------------------------------------------------------------
data("Fst_dat")
data("Het_dat")
data("Q_dat")
data("HornedLizard_Pop")
data("HornedLizard_VCF")

## ----Heterozygosity, echo=TRUE, eval=FALSE------------------------------------
#  Obs_Het <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Ho")

## ----Differentiation, echo=TRUE, eval=FALSE-----------------------------------
#  Fst <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Fst")

## ----Private.alleles, echo=TRUE, eval=FALSE-----------------------------------
#  PA <- Private.alleles(data = HornedLizard_VCF, pops = HornedLizard_Pop)

## ----Ancestry barchart, echo=TRUE, eval=FALSE---------------------------------
#  # First, we seperate the list elements into two seperate objects. The q-matrix (Qmat) and the locality information for each individual (Loc).
#  Qmat <- Q_dat[[1]]
#  Loc <- Q_dat[[2]]
#  
#  # Now we will generate both population and individual plots by setting plot.type to 'all'. If you wanted, you could only generate individual or population plots by setting plot.type to "individual" and "population", respectively.
#  Test_all <- Ancestry_barchart(anc.mat = Qmat, pops = Loc, K = 5,
#  plot.type = 'all', col = c('#d73027', '#f46d43', '#e0f3f8', '#74add1', '#313695'))
#  
#  Test_all$`Individual Ancestry Plot`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Anc_barchart1.png")

## ----eval=FALSE---------------------------------------------------------------
#  Test_all$`Population Ancestry Plot`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Anc_barchart2.png")

## ----Piechart map, echo=TRUE, eval=FALSE, fig.align='center'------------------
#  # First, we seperate the list elements into two seperate objects. The q-matrix (Qmat) and the locality information for each individual (Loc).
#  Qmat <- Q_dat[[1]]
#  Loc <- Q_dat[[2]]
#  
#  # Now we will generate both population and individual plots by setting plot.type to 'all'. If you wanted, you could only generate individual or population plots by setting plot.type to "individual" and "population", respectively.
#  Test_all_piemap <- Piechart_map(anc.mat = Qmat, pops = Loc, K = 5,plot.type = 'all', col = c('#d73027', '#f46d43', '#e0f3f8', '#74add1', '#313695'),
#                                  Lat_buffer = 1, Long_buffer = 1)
#  
#  Test_all_piemap$`Individual Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Ind_PieMap.png")

## ----Piechart map2, echo=TRUE, eval=FALSE, fig.align='center'-----------------
#  Test_all_piemap$`Population Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Pop_PieMap.png")

## ----Differentiation_old, echo = TRUE, eval=FALSE-----------------------------
#  # Isolate our fst matrix and locality information
#  Fst <- Fst_dat[[1]]
#  Loc <- Fst_dat[[2]]
#  # Plot the heatmap, the statistic argument is used to label the plot.
#  Fstat_hmap <- Pairwise_heatmap(dat = Fst, statistic = 'FST')
#  # Look at the plot
#  Fstat_hmap

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE-----------
knitr::include_graphics("./img/Differentiation-1.png")

## ----Differentiation map, echo = TRUE, eval=FALSE-----------------------------
#  # Closest Neighbor
#  Fst_map <- Dif_Stats_Map(dat = Fst, pops = Loc, neighbors = 1,
#                       col = c('#FFFF00','#FFC000','#FFA500','#e31a1c','#800026'),
#                       Lat_buffer = 1, Long_buffer = 1)
#  
#  Fst_map2 <- Dif_Stats_Map(dat = Fst, pops = Loc, neighbors = c('East_West', 'East_South', 'South_West'),
#                       col = c('#FFFF00','#FFC000','#FFA500','#e31a1c','#800026'),
#                       Lat_buffer = 1, Long_buffer = 1)
#  
#  Fst_map$`Differentiation Map`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE-----------
knitr::include_graphics("./img/Differentiationmap-1.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  Fst_map2$`Differentiation Map`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE-----------
knitr::include_graphics("./img/Differentiationmap-2.png")

## ----Diversity, echo = TRUE, eval=FALSE---------------------------------------
#  # Similar to our heat map, we use the statistic argument to label our figures and any output raster you create.
#  Div_map <- Div_Stats_Map(dat = Het_dat, plot.type = 'point',
#  statistic = "Heterozygosity", col = c('blue', 'orange', 'red'), Lat_buffer = 1, Long_buffer = 1, prefix = 'Test_het')
#  
#  Div_map$`Heterozygosity Map`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE-----------
knitr::include_graphics("./img/Diversity-1.png")

