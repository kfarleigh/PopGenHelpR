# A function to map statistics as colored points on a map.

A function to map statistics as colored points on a map.

## Usage

``` r
Point_map(
  dat,
  statistic,
  size = 3,
  breaks = NULL,
  col,
  out.col = NULL,
  Lat_buffer = 1,
  Long_buffer = 1,
  Latitude_col = NULL,
  Longitude_col = NULL,
  country_code = NULL,
  shapefile = NULL,
  raster = NULL,
  legend_pos = "none",
  scale_bar = FALSE,
  north_arrow = FALSE,
  north_arrow_style = ggspatial::north_arrow_nautical(),
  north_arrow_position = NULL,
  shapefile_plot_position = NULL,
  raster_plot_position = NULL,
  shapefile_col = NULL,
  shapefile_outline_col = NULL,
  shp_outwidth = 1,
  raster_col = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
  interpolate_raster = NULL,
  raster_breaks = NULL,
  discrete_raster = NULL
)
```

## Arguments

- dat:

  Data frame or character string that supplies the input data. If it is
  a character string, the file should be a csv. The first column should
  be the statistic to be plotted. The coordinates of each row should be
  indicated by columns named Longitude and Latitude. Alternatively, see
  the Longitude_col and Latitude_col arguments.

- statistic:

  Character string. The statistic to be plotted.

- size:

  Numeric. The size of the points to plot.

- breaks:

  Numeric. The breaks used to generate the color ramp when plotting.
  Users should supply 3 values if custom breaks are desired.

- col:

  Character vector indicating the colors you wish to use for plotting,
  three colors are allowed (low, mid, high). The first color will be the
  low color, the second the middle, the third the high.

- out.col:

  Character. A color for outlining points on the map. There will be no
  visible outline if left as NULL.

- Lat_buffer:

  Numeric. A buffer to customize visualization.

- Long_buffer:

  Numeric. A buffer to customize visualization.

- Latitude_col:

  Numeric. The number of the column indicating the latitude for each
  sample. If this is not null, PopGenHelpR will use this column instead
  of looking for the Latitude column.

- Longitude_col:

  Numeric. The number of the column indicating the longitude for each
  sample. If this is not null, PopGenHelpR will use this column instead
  of looking for the Longitude column.

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

- raster:

  Character.A file name or a spatraster object that is compatible with
  the terra R package. This should be used in conjunction with the
  raster_plot_position argument.

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
  3, which plots the shapefile on top of everything.

- raster_plot_position:

  Numeric. A number indicating which position to plot the shapefile in.
  The options are 1, which plots the raster on top of the base world map
  (under points and administrative boundaries), 2 which plots the raster
  on top of administrative boundaries (but under points), and 3, which
  plots the raster on top of everything.

- shapefile_col:

  Character. A color or color vector indicating the color to fill the
  shapefile(s) with. Similar to `group_col`, shapefiles will be colored
  alphabetically.

- shapefile_outline_col:

  Character. A color indicating the outline color of the shapefile.

- shp_outwidth:

  Numeric. The width of the shapefile outline.

- raster_col:

  Character. A character vector indicating the colors used to visualize
  the raster. The function will seperate your raster data into the same
  number of bins as there are colors. If you provide 5 colors, for
  example, there will be 5 bins.

- interpolate_raster:

  Boolean. Whether or not to interpolate the raster. The default is to
  interpolate the raster.

- raster_breaks:

  Numeric or Character vector. Values to be used as breaks for the
  raster surface.

- discrete_raster:

  Boolean. Indicating whether or not the raster being supplied is
  discrete.

## Value

A list containing maps and the data frames used to generate them.

## Author

Keaka Farleigh

## Examples

``` r
# \donttest{
data(Het_dat)
Test <- Point_map(Het_dat, statistic = "Heterozygosity")# }
```
