## code to prepare `Het_dat` dataset goes here

usethis::use_data(Het_dat, overwrite = TRUE)

Het_dat <- read.csv('Test_Heterozygosity.csv')

use_data(Het_dat)
