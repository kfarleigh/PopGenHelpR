# A function to plot a heatmap from a symmetric matrix.

A function to plot a heatmap from a symmetric matrix.

## Usage

``` r
Pairwise_heatmap(
  dat,
  statistic,
  col = c("#abd9e9", "#2c7bb6", "#ffffbf", "#fdae61", "#d7191c"),
  breaks = NULL
)
```

## Arguments

- dat:

  Data frame or character string that supplies the input data. If it is
  a character string, the file should be a csv. If it is a csv, the 1st
  row should contain the individual/population names. The columns should
  also be named in this fashion.

- statistic:

  Character indicating the statistic represented in the matrix, this
  will be used to label the plot.

- col:

  Character vector indicating the colors to be used in plotting. The
  vector should contain two colors, the first will be the low value, the
  second will be the high value.

- breaks:

  Numeric. The breaks used to generate the color ramp when plotting. The
  number of breaks should match the number of colors.

## Value

A heatmap plot

## Examples

``` r
# \donttest{
#' data(Fst_dat)
Fst <- Fst_dat[[1]]
Fstat_plot <- Pairwise_heatmap(dat = Fst, statistic = 'FST')# }
```
