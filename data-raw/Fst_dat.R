## code to prepare `Fst_dat` dataset goes here

usethis::use_data(Fst_dat, overwrite = TRUE)

# Read in example Fst output
Fst_dat <- read.csv('Test_Fst.csv', row.names = 1)

# Read in example locality information
Loc_dat <- read.csv('Platyrhinos_Locdat.csv')
# Select one row from each cluster
Loc_dat <- Loc_dat[c(1,23,56),]
# Rename columns to adhere to formatting requirements
colnames(Loc_dat) <- c('Sample', 'Population', 'Long', 'Lat')
# Store in a list
Fst_dat <- list(Fst_dat, Loc_dat)

use_data(Fst_dat)
