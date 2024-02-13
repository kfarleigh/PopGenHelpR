#'\bold{WARNING! This function has been deprecated and is no longer supported. Please use the Point_map function instead.}
#'A function to map diversity statistics.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the statistic to be plotted. The coordinates of each row should be indicated by columns named Longitude and Latitude.
#' @param statistic Character string. The statistic to be plotted.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. Users should supply 3 values if custom breaks are desired.
#' @param col Character vector indicating the colors you wish to use for plotting, three colors are allowed (low, mid, high). The first color will be the low color, the second the middle, the third the high.
#' @param prefix  Character string that will be appended to file output.
#' @param Long_col Numeric. A number indicating which column contains the longitude information.
#' @param Lat_col Numeric. A number indicating which column contains the latitude information.
#' @param write Boolean. Whether or not to write the output to a file in the current working directory.
#'
#' @return A list containing maps and the data frames used to generate them.
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(Het_dat)
#' Test <- Point_map(Het_dat, statistic = "Heterozygosity")}
Point_map <- function(dat, statistic, breaks = NULL, col, Lat_buffer = 1, Long_buffer = 1){
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
  Locs <- Div_mat[,4:5]

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

    ### Heterozygosity Map
    # Map it with colored points
    Div_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = .data$Longitude, y = .data$Latitude, color = Div_mat[,1]), shape = 19, size = 3) +
      ggplot2::scale_color_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = mean(Div_mat[,1]), breaks = Breaks, name = statistic) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") +
      ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

  Output <- list(Div_map, Div_mat)
    names(Output) <- c(paste(statistic," Map", sep = ""), "Plotting Dataframe")
    return(Output)

}
