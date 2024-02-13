# PopGenHelpR 1.3.0
- December 6th, 2023
** To Do for KF
- Add tests

* The `Plot_coordinates` function has been added to make sample maps from coordinates. 
* The `Point_map` function has replaced `Div_stats_map` and a statistic argument has been added to highlight the functions utility and allow users to name the map legend.
* `Dif_stats_map` has been changed to `Network_map` and a statistic argument has been added to highlight the functions utility and allow users to name the map legend.
* `Plot_ancestry` has been split into `Piechart_Map` and `Ancestry_barchart` so that it is easier for users to determine which function is most appropriate for their analysis. 
* `Differentiation` has been added to estimate Fst, Nei's D, and Jost's D. Please see the documentation for more details. 
* `Heterozygosity` has been added to estimate 7 different measures of heterozygosity. Please see the documentation for more details. 
* `Private.alleles` has been added to calculate the number of private alleles in each population. 
* The `Dif_Stats` function has been deprecated, please used the `Differentiation` function to calculate pairwise differentiation between populations (Fst, Nei's D, Jost's D) or individuals (Nei's D). 
* The `Div_Stats` function has been deprecated, please use the `Heterozygosity` function if you wish to estimate heterozygosity or the `Private.alleles` function if you wish to calculate the number of private alleles per population. Please use the `Point_Map` function if you wish to visualize the results on a map or plot.

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
