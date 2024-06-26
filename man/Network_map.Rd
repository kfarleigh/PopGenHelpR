% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Network_map.R
\name{Network_map}
\alias{Network_map}
\title{A function to map statistics (i.e., genetic differentiation) between points as a network on a map.}
\usage{
Network_map(
  dat,
  pops,
  neighbors,
  col,
  statistic = NULL,
  breaks = NULL,
  Lat_buffer = 1,
  Long_buffer = 1,
  Latitude_col = NULL,
  Longitude_col = NULL
)
}
\arguments{
\item{dat}{Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. If it is a csv, the 1st row should contain the individual/population names. The columns should also be named in this fashion.}

\item{pops}{Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The columns should be named Sample, containing the sample IDs; Population indicating the population assignment of the individual; Long, indicating the longitude of the sample; Lat, indicating the latitude of the sample. Alternatively, see the Longitude_col and Latitude_col arguments.}

\item{neighbors}{Numeric or character. The number of neighbors to plot connections with, or the specific relationship that you want to visualize. Names should match those in the population assignment file and be seperated by an underscore. If I want to visualize the relationship between East and West, for example, I would set neighbors = "East_West".}

\item{col}{Character vector indicating the colors you wish to use for plotting.}

\item{statistic}{Character indicating the statistic being plotted. This will be used to title the legend. The legend title will be blank if left as NULL.}

\item{breaks}{Numeric. The breaks used to generate the color ramp when plotting. Users should supply 3 values if custom breaks are desired.}

\item{Lat_buffer}{Numeric. A buffer to customize visualization.}

\item{Long_buffer}{Numeric. A buffer to customize visualization.}

\item{Latitude_col}{Numeric. The number of the column indicating the latitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Lat column.}

\item{Longitude_col}{Numeric. The number of the column indicating the longitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Long column.}
}
\value{
A list containing the map and the matrix used to plot the map.
}
\description{
A function to map statistics (i.e., genetic differentiation) between points as a network on a map.
}
\examples{
\donttest{
data(Fst_dat)
Fst <- Fst_dat[[1]]
Loc <- Fst_dat[[2]]
Test <- Network_map(dat = Fst, pops = Loc,
neighbors = 2,col = c('#4575b4', '#91bfdb', '#e0f3f8','#fd8d3c','#fc4e2a'),
statistic = "Fst", Lat_buffer = 1, Long_buffer = 1)}
}
\author{
Keaka Farleigh
}
