# Plot a map of ancestry pie charts.

Plot a map of ancestry pie charts.

## Usage

``` r
Piechart_map(
  anc.mat,
  pops,
  K,
  plot.type = "all",
  col,
  piesize = 0.35,
  Lat_buffer,
  Long_buffer,
  Latitude_col = NULL,
  Longitude_col = NULL,
  country_code = NULL,
  shapefile = NULL,
  legend_pos = "none",
  scale_bar = FALSE,
  north_arrow = FALSE,
  north_arrow_style = ggspatial::north_arrow_nautical(),
  north_arrow_position = NULL,
  shapefile_plot_position = NULL,
  shapefile_col = NULL,
  shapefile_outline_col = NULL,
  shp_outwidth = 1
)
```

## Arguments

- anc.mat:

  Data frame or character string that supplies the input data. If it is
  a character string, the file should be a csv. The first column should
  be the names of each sample/population, followed by the estimated
  contribution of each cluster to that individual/pop.

- pops:

  Data frame or character string that supplies the input data. If it is
  a character string, the file should be a csv. The columns should be
  named Sample, containing the sample IDs; Population indicating the
  population assignment of the individual, population and sample names
  must be the same type (i.e., both numeric or both characters); Long,
  indicating the longitude of the sample; Lat, indicating the latitude
  of the sample. Alternatively, see the Longitude_col and Latitude_col
  arguments.

- K:

  Numeric.The number of genetic clusters in your data set, please
  contact the package authors if you need help doing this.

- plot.type:

  Character string. Options are all, individual, and population. All is
  default and recommended, this will plot a piechart map for both the
  individuals and populations.

- col:

  Character vector indicating the colors you wish to use for plotting.

- piesize:

  Numeric. The radius of the pie chart for ancestry mapping.

- Lat_buffer:

  Numeric. A buffer to customize visualization.

- Long_buffer:

  Numeric. A buffer to customize visualization.

- Latitude_col:

  Numeric. The number of the column indicating the latitude for each
  sample. If this is not null, PopGenHelpR will use this column instead
  of looking for the Lat column.

- Longitude_col:

  Numeric. The number of the column indicating the longitude for each
  sample. If this is not null, PopGenHelpR will use this column instead
  of looking for the Long column.

- country_code:

  Character. A country code or vector of country codes from the R
  package [geodata](https://cran.r-project.org/package=geodata)
  specifying the country that you want to plot administrative borders
  for (e.g, US states). You can determine the correct codes using
  geodata's `country_codes` function.

- shapefile:

  Character. A file name, vector of file names of a shapefile(s) to plot
  on the map, or a spatvector object that is compatible with the R
  package terra. This should be used in conjunction with the
  shapefile_plot_position argument.

- legend_pos:

  Character. The desired position of the legend. The default is "none",
  which removes the legend. Other options include "left", "right", "top"
  or "bottom". Please see the ggplot2 documentation for all of the
  legend placement options.

- scale_bar:

  Boolean. Whether or not to add a scale bar. Note that maps with large
  areas or those that use unprojected spatial data (i.e., WGS 84) will
  generate a warning that the scale bar varies.

- north_arrow:

  Boolean. Whether or not to add a north arrow.

- north_arrow_style:

  Character. Which style of north arrow to add. See
  [ggspatial](https://cran.r-project.org/package=ggspatial)
  documentation for more details.

- north_arrow_position:

  Character. The position of the north arrow. See
  [ggspatial](https://cran.r-project.org/package=ggspatial)
  documentation for more details.

- shapefile_plot_position:

  Numeric. A number indicating which position to plot the shapefile in.
  The options are 1, which plots the shapefile on top of the base world
  map (under points and administrative boundaries), 2 which plots the
  shapefile on top of administrative boundaries (but under points), and
  3, which plots the shapefile on top of everything. Note that an error
  will occur if you do not provide a `country_code` argument.

- shapefile_col:

  Character. A color or color vector indicating the color to fill the
  shapefile(s) with. Similar to `group_col`, shapefiles will be colored
  alphabetically.

- shapefile_outline_col:

  Character. A color indicating the outline color of the shapefile.

- shp_outwidth:

  Numeric. The width of the shapefile outline.

## Value

A list containing your plots and the data frames used to generate the
plots.

## Author

Keaka Farleigh

## Examples

``` r
# \donttest{
data(Q_dat)
Qmat <- Q_dat[[1]]
rownames(Qmat) <- Qmat[,1]
Loc <- Q_dat[[2]]
Test_all <- Piechart_map(anc.mat = Qmat, pops = Loc, K = 5,
plot.type = 'all', col = c('#d73027', '#fc8d59', '#e0f3f8', '#91bfdb', '#4575b4'), piesize = 0.35,
Lat_buffer = 1, Long_buffer = 1)# }
#> Warning: There was 1 warning in `dplyr::summarise()`.
#> ℹ In argument: `dplyr::across(1:(K + 2), mean, na.rm = TRUE)`.
#> ℹ In group 1: `Pop = "1"`.
#> Caused by warning:
#> ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
#> Supply arguments directly to `.fns` through an anonymous function instead.
#> 
#>   # Previously
#>   across(a:b, mean, na.rm = TRUE)
#> 
#>   # Now
#>   across(a:b, \(x) mean(x, na.rm = TRUE))
#> ℹ The deprecated feature was likely used in the PopGenHelpR package.
#>   Please report the issue at <https://github.com/kfarleigh/PopGenHelpR/issues>.
```
