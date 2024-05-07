#' A function to plot coordinates on a map.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The coordinates of each row should be indicated by columns named Longitude and Latitude. Alternatively, see the Latitude_col and Longitude_col arugments.
#' @param col Character vector indicating the colors you wish to use for plotting, two colors are allowed. The first color will be the fill color, the second is the outline color. For example, if I want red points with a black outline I would set col to col = c("#FF0000", "#000000").
#' @param size Numeric. The size of the points to plot.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#' @param Latitude_col Numeric. The number of the column indicating the latitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Latitude column.
#' @param Longitude_col Numeric. The number of the column indicating the longitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Longitude column.
#'
#' @return A ggplot object.
#'
#'
#' @importFrom rlang .data
#'
#' @author Keaka Farleigh
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' Test <- Plot_coordinates(HornedLizard_Pop)}
Plot_coordinates <- function(dat, col = c("#A9A9A9", "#000000"), size = 3, Lat_buffer = 1, Long_buffer = 1, Latitude_col = NULL, Longitude_col = NULL){
  Long <- Lat <- x <- y <- z <- alpha <- world <- NULL
  ################### Get the data for mapping
  # Get map data
  map <- spData::world["continent"]
  states <- spData::us_states["REGION"]
  ### Make a base map for the countries of interest
  base_map <- ggplot2::ggplot() + ggplot2::geom_sf(data = map, fill = "#f4f4f4") +
    ggplot2::geom_sf(data = states, fill = ggplot2::alpha("#f4f4f4", 0))
  # Read in files
  if(is.data.frame(dat) == TRUE){
    Div_mat <- dat
  }
  else if(is.character(dat) == TRUE){
    Div_mat <- utils::read.csv(dat)
  }
  else{
    stop("Please supply a dataframe or .csv file name for analysis")
  }

  if(!is.null(Latitude_col)){
    colnames(Div_mat)[Latitude_col] <- "Latitude"
  }

  if(!is.null(Longitude_col)){
    colnames(Div_mat)[Longitude_col] <- "Longitude"
  }

  ### Get coordinate ranges for our data
  Lat_Min <- min(Div_mat$Latitude) - Lat_buffer
  Lat_Max <- max(Div_mat$Latitude) + Lat_buffer
  Long_Min <- min(Div_mat$Longitude) - Long_buffer
  Long_Max <- max(Div_mat$Longitude) + Long_buffer

 # Set colors
 Fill_col <- col[1]
 Outline_col <- col[2]

    ### Heterozygosity Map
    # Map it with colored points
    map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = .data$Longitude, y = .data$Latitude), shape = 21, size = 3, fill = Fill_col, color = Outline_col) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = NULL) +
      ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

  Output <- map
  return(Output)

}
