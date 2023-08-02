## code to prepare `HornedLizard` dataset goes here

usethis::use_data(HornedLizard, overwrite = TRUE)

HornedLizard_VCF <- vcfR::read.vcfR('P_platyrhinos_SingleSNPperLoc_100.recode.vcf')
HornedLizard_Pop <- read.delim('P_platyrhinos_SingleSNPperLoc_100.012.indv', header = FALSE)

use_data(HornedLizard_VCF, HornedLizard_Pop)
