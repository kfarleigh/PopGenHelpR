#' A function to map statistics as colored points on a map.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the statistic to be plotted. The coordinates of each row should be indicated by columns named Longitude and Latitude.
#' @param statistic Character string. The statistic to be plotted.
#' @param size Numeric. The size of the points to plot.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. Users should supply 3 values if custom breaks are desired.
#' @param col Character vector indicating the colors you wish to use for plotting, three colors are allowed (low, mid, high). The first color will be the low color, the second the middle, the third the high.
#' @param out.col Character. A color for outlining points on the map. There will be no visible outline if left as NULL.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#'
#' @return A list containing maps and the data frames used to generate them.
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
#' data(Het_dat)
#' Test <- Point_map(Het_dat, statistic = "Heterozygosity")}
Point_map <- function(dat, statistic, size = 3, breaks = NULL, col, out.col = NULL, Lat_buffer = 1, Long_buffer = 1){
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

  ### Get coordinate ranges for our data
  Lat_Min <- min(Div_mat$Latitude) - Lat_buffer
  Lat_Max <- max(Div_mat$Latitude) + Lat_buffer
  Long_Min <- min(Div_mat$Longitude) - Long_buffer
  Long_Max <- max(Div_mat$Longitude) + Long_buffer

  # Set colors
  # Set breaks
  if(missing(col)){
    col <- c('#4575b4','#fee090','#d73027')
  }
  if(is.null(breaks) == TRUE){
    Breaks <- summary(Div_mat[,1])
    Breaks <- as.numeric(Breaks[1:5])
  }
  else{
    Breaks <- breaks
    Breaks <- as.numeric(Breaks)
  }
  Breaks <- round(Breaks,2)
  size <- size

  if(is.null(out.col)){
    out.col <- "#f4f4f4"
  } else{
    out.col <- out.col
  }

  ### Heterozygosity Map
  # Map it with colored points
  Div_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
    ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = .data$Longitude, y = .data$Latitude, fill = Div_mat[,1]), shape = 21, size = size, color = out.col) +
    ggplot2::scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = mean(Div_mat[,1]), breaks = Breaks, name = statistic) +
    ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") +
    ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

  Output <- list(Div_map, Div_mat)
  names(Output) <- c(paste(statistic," Map", sep = ""), "Plotting Dataframe")
  return(Output)

}
