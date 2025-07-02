# PopGenHelpR 1.4.0
- July 2nd, 2025
## `Differentiation` and `Heterozygosity`
- Functions can now detect whether "|" or "\" separate genotypes. This was causing functions to throw NAs.
- Minor bug was fixed for Jost's D calculation. A parentheses was removed in the calculation of the `D.bar` object. 

## `Ancestry_barchart`
- Code has been simplified, individual and population information is now coerced to character before visualization.
- It is now an error if the individual names in the population data do not match the names in the individual data.
- Plots are automatically sorted by cluster so that individuals and populations are plotted with others in the cluster, if there is no specified order.
- Users can specify the individual and population order using the `ind.order` and `pop.order` arguments, respectively. 
- Users can specify the legend position using the `legend_pos` argument. 

## Mapping functions 
The following updates laregly apply for the `Plot_coordinates`, `Piechart_map`, `Point_map`, and `Network_map` functions. We note where any updates are function specific. 

- The maps are now built using the R package geodata to avoid the mapping NULL layer error.
- Users can now include administrative borders like states using the `country_code` argument. Country codes must match those used in geodata.
- Users can color points by group using the `group` argument and `group_col` argument to specify colors for each group (only for `Plot_coordinates`).
- Users can specify the legend position using the `legend_pos` argument. 
- Users can add a scale bar or north arrow to plots using the `scale_bar` and `north_arrow` arguments. Users can also modify the position and style of the north arrow using `north_arrow_position` and `north_arrow_style` arguments, respectively. 
- Users can add a shapefile to the map with the `shapefile` argument. Users can also control where the shapefile is placed using `shape_file_position`, the fill color (`shapefile_col`), and the outline color (`shapefile_outline_col`). Shapefiles can be character strings or spatVector objects. 
- Users can add a raster (discrete or continuous) and control it's placement in a similar fashion to shapefiles, with the `raster`, `raster_col`, and `raster_plot_position` arguments. Custom breaks can be set with `raster_breaks` and  users can indicate if rasters are discrete using `discrete_raster`. **Note that sometimes `terra` interprets continuous rasters as discrete, and users should use the `discrete_raster` argument to accommodate this.** Rasters can be characters or spatRaster objects. 
- *Rasters cannot be added to the `Piechart_map output because it would require additional dependencies to accomodate two fills on the same ggplot2 object. Please email Keaka Farleigh (see package maintainer) if this poses a problem.*

## `Pairwise_heatmap`
- Users can now specify color gradients of their preferred length (as long as they have that many colors).

# PopGenHelpR 1.3.2
- July 31st, 2024
* Code to calculate Fst in the `Differentiation` function has been updated to be consistent with the rest of the function (using !is.na).

# PopGenHelpR 1.3.1
- July 26th, 2024
* The previously deprecated `Dif_stats`, `Dif_stats_map`, `Div_stats`, `Div_stats_map`, and `Plot_Ancestry` functions have been removed.

# PopGenHelpR 1.3.0
- May 7th, 2024
* The `Longitude_col` and `Latitude_col` arguments have been added to mapping functions; `Plot_coordinates`, `Ancestry_Piemap`, `Piechart_map`, `Network_map`, and `Point_map`.
* `Ind_order` and `Pop_order` arguments have been added to `Ancestry_Barchart` to allow users to specify the order of individuals and populations in the barchart.
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
