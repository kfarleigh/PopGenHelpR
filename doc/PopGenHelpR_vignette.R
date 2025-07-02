## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Load the package
library(PopGenHelpR)

## ----out.width= "750px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/PGH_workflow.png")

## ----VCF filtering, echo=TRUE, eval=FALSE-------------------------------------
# # vcftools
# vcftools --vcf myfile.vcf --max-alleles 2 --recode --recode-INFO-all --out my_biallelic_file.vcf
# 
# # bcftools
# bcftools view -m2 -M2 -v snps myfile.vcf > my_biallelic_file.vcf

## ----load data----------------------------------------------------------------
data("Fst_dat")
data("Het_dat")
data("Q_dat")
data("HornedLizard_Pop")
data("HornedLizard_VCF")

## ----Heterozygosity, echo=TRUE, eval=FALSE------------------------------------
# Obs_Het <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Ho")

## ----Differentiation, echo=TRUE, eval=FALSE-----------------------------------
# Fst <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Fst")

## ----Private.alleles, echo=TRUE, eval=FALSE-----------------------------------
# PA <- Private.alleles(data = HornedLizard_VCF, pops = HornedLizard_Pop)

## ----Ancestry barchart, echo=TRUE, eval=FALSE---------------------------------
# # First, we separate the list elements into two separate objects. The q-matrix (Qmat) and the locality information for each individual (Loc).
# Qmat <- Q_dat[[1]]
# Loc <- Q_dat[[2]]
# 
# # Now we will generate both population and individual plots by setting plot.type to 'all'. If you wanted, you could only generate individual or population plots by setting plot.type to "individual" and "population", respectively.
# Test_all <- Ancestry_barchart(anc.mat = Qmat, pops = Loc, K = 5,
# plot.type = 'all', col = c('#d73027', '#f46d43', '#e0f3f8', '#74add1', '#313695'))
# 
# Test_all$`Individual Ancestry Plot`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Anc_barchart1.png")

## ----eval=FALSE---------------------------------------------------------------
# Test_all$`Population Ancestry Plot`

## ----out.width= "600px", out.height= "350px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Anc_barchart2.png")

## ----Piechart map, echo=TRUE, eval=FALSE, fig.align='center'------------------
# # First, we seperate the list elements into two seperate objects. The q-matrix (Qmat) and the locality information for each individual (Loc).
# Qmat <- Q_dat[[1]]
# Loc <- Q_dat[[2]]
# 
# # Now we will generate both population and individual plots by setting plot.type to 'all'. If you wanted, you could only generate individual or population plots by setting plot.type to "individual" and "population", respectively.
# Test_all_piemap <- Piechart_map(anc.mat = Qmat, pops = Loc, K = 5,plot.type = 'all', col = c('#d73027', '#f46d43', '#e0f3f8', '#74add1', '#313695'),
#                                 Lat_buffer = 1, Long_buffer = 1)
# 
# Test_all_piemap$`Individual Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Ind_PieMap.png")

## ----Piechart map2, echo=TRUE, eval=FALSE, fig.align='center'-----------------
# Test_all_piemap$`Population Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Pop_PieMap.png")

## ----Pairwise Heatmap, echo=TRUE, eval=FALSE----------------------------------
# PW_hmap <- Pairwise_heatmap(Fst_dat[[1]], statistic = "Fst", col = c("#0000FF", "#FF0000"))

## ----out.width= "600px", out.height= "600px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/PWhmap.png")

## ----Network Map, echo=TRUE, eval=FALSE---------------------------------------
# NW_map <- Network_map(Fst_dat[[1]], pops = Fst_dat[[2]], neighbors = 2, statistic = "Fst")
# NW_map$Map

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/NW_map1.png")

## ----Network Map2, echo=TRUE, eval=FALSE--------------------------------------
# NW_map2 <- Network_map(Fst_dat[[1]], pops = Fst_dat[[2]], neighbors = c("East_West", "East_South"), statistic = "Fst")
# NW_map2$Map

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/NW_map2.png")

## ----Heterozygosity2, echo=TRUE, eval=FALSE-----------------------------------
# Het_map <- Point_map(Het_dat, statistic = "Heterozygosity")
# Het_map$`Heterozygosity Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Het_map1.png")

## ----Heterozygosity3, echo=TRUE, eval=FALSE-----------------------------------
# Het_map2 <- Point_map(Het_dat, statistic = "Heterozygosity", out.col = "#000000")
# Het_map2$`Heterozygosity Map`

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Het_map2.png")

## ----Coordinates, echo=TRUE, eval=FALSE---------------------------------------
# Sample_map <- Plot_coordinates(HornedLizard_Pop)
# Sample_map

## ----out.width= "500px", out.height= "750px", echo=FALSE, eval=TRUE, fig.align='center'----
knitr::include_graphics("./img/Samp_map1.png")

