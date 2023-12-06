# PopGenHelpR 1.3.0
- December 6th, 2023
** To Do for KF
- Split `Div_Stats` into `Heterozygosity` and `Private_alleles`
- Split `Dif_Stats` into `Fst` and `NeisD`
- Rewrite all functions to allow 012/geno files
- Rewrite all functions to auto detect vcf vs 012/geno
- Rewrite error to warning for vcf/popmap sorting 
- Rename mapping functions
- Split `Plot_Ancestry` into `Piechart_Map` and `Ancestry_barchart`
- Add `Nucleotide_diversity`
- Add expected Heterozygosity (He) to `Heterozygosity`
- Add proportion of heterozygous loci (PHt) to `Heterozygosity`
- Add heterozygosity standardized by expected heterozygosity (StHe) to `Heterozygosity`
- Add heterozygosity standardized by observed heterozygosity (StHo) to `Heterozygosity`
- Add internal relatedness (IR) to `Heterozygosity`
- Add homozygosity by locus (HL) to `Heterozygosity`
- Update documentation
- Update vignette 

# PopGenHelpR 1.2.2
- October 2nd, 2023
* `Dif_stats_Map`, `Div_Stats_Map`, and `Plot_ancestry` have been updated to use a base color of #f4f4f4 instead of grey99 which was throwing an error for some users. 
* The piesize argument has been added to `Plot_ancestry`, the original value of 0.35 was found to be too high, especially in cases where users were mapping a smaller geographic area. 

# PopGenHelpR 1.2.1
- August 14th, 2023
* `Div_Stats` and `Dif_stats` have been updated to accept a vcf file or vcfR object as input. 
* `Div_Stats` and `Dif_stats` have been updated to accept a csv file or data frame for population assignment. 
* `Plot_ancestry` has been updated to generate structure-like plots using ggplot2 instead of base R and to handle character and numeric values for individual and population names. Note that individual and populations must be of the same type (i.e., both numeric or both characters). 
* The dependency `rnaturalearth` is no longer used. We now use `spData` for mapping data. 
* The vignette has been updated to accommodate the changes noted above. 

# PopGenHelpR 1.1.1 
- July 17th, 2023
* The horned lizard data was added so that examples can be run by users.
* The write argument has been added to `Div_Stats` and `Dif_stats` so that files are not automatically written to the working directory. 

# PopGenHelpR 1.0.1 
- July 17th, 2023
*  PopGenHelpR has been updated to rnaturalearthhires in Suggests field of the DESCRIPTION file and now use it conditonally. 


# PopGenHelpR 1.0.0
*  First development of PopGenHelpR, publication on Github, and submission to the CRAN (02/06/2023)
