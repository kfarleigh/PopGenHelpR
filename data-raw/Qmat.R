## code to prepare `Qmat` dataset goes here

usethis::use_data(Qmat, overwrite = TRUE)

Qmat <- read.csv('Test_Qmat.csv')
Qmat <- Qmat[1:30,]
Loc_dat <- read.csv('Test_Popmap.csv')

Q_dat <- list(Qmat, Loc_dat)

use_data(Q_dat)
