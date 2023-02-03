#' A function to map diversity statistics.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the statistic to be plotted and named the same as the statistic argument. The second column is Population indicating which population each row belongs to. The third column is the standard deviation, the fourth column is Long indicating the longitude, and the fifth column is Lat, indicating the latitude.
#' @param plot.type Character string. Options are all, point, or interpolated. All is recommended and will generate a map with points colored according to heterozygosity as well as a rater of interpolated heterozygosity values.
#' @param statistic Character string. The statistic to be plotted.
#' @param countries Character vector indicating the country borders that you wish to plot on a map. Can be any country that is valid in the ne_states function in the rnaturalearth package.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. Users should supply 3 values if custom breaks are desired.
#' @param col Character vector indicating the colors you wish to use for plotting, three colors are allowed (low, mid, high). The first color will be the low color, the second the middle, the third the high.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#' @param prefix  Character string that will be appended to file output.
#'
#' @return A list containing maps and the data frames used to generate them.
#' @export
#'
#' @examples
#' \dontrun{
#' data(Het_dat)
#' Test_het <- Div_Stats_Map(dat = Het_dat, plot.type = 'all',
#' statistic = "Heterozygosity", countries = c('united states of america', 'mexico'),
#' Lat_buffer = 1, Long_buffer = 1, prefix = 'Test_het')}
Div_Stats_Map <- function(dat, plot.type = 'all', statistic, countries, breaks = NULL, col, Lat_buffer = 1, Long_buffer = 1, prefix = NULL){
  Long <- Lat <- x <- y <- z <- NULL
  ################### Get the data for mapping
  # Get map data
  if(is.null(prefix) == TRUE){
    stop("Please supply a prefix for output files")
  }
  countries <- countries
  Country_borders <- list()
  for(i in 1:length(countries)) {
    ncon <- 1:length(countries)
    Country_borders[[i]] <- rnaturalearth::ne_states(country = countries[i], returnclass = 'sf')
    names(Country_borders)[[i]] <- paste("Country", ncon[i], sep = "_")
  }
  world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
  ### Make a base map for the countries of interest
  base_map <-  ggplot2::ggplot() + ggplot2::geom_sf(data = world, fill = 'grey99')
  for(i in 1:length(countries)) {
    base_map <- base_map + ggplot2::geom_sf(data = Country_borders[[i]], fill = "grey99")
  }

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
  Lat_Min <- min(Div_mat$Lat) - Lat_buffer
  Lat_Max <- max(Div_mat$Lat) + Lat_buffer
  Long_Min <- min(Div_mat$Long) - Long_buffer
  Long_Max <- max(Div_mat$Long) + Long_buffer

  # Set colors
  # Set breaks
  if(missing(col)){
    col <- c('#4575b4','#ffffbf','#d73027')
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
  if(plot.type == 'all') {
    ### Heterozygosity Map
    # Map it with colored points
    Div_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = Long, y = Lat, color = Div_mat[,1]), shape = 19, size = 3) +
      ggplot2::scale_color_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = mean(Div_mat[,1]), breaks = Breaks, name = statistic) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") +
      ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    # Interpolate the values
    sp::coordinates(Locs) <- c('Long', 'Lat')
    grd <- as.data.frame(sp::spsample(Locs,'regular', n=100000))
    names(grd) <- c("X", "Y")
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object
    P.idw <- gstat::idw(Div_mat[,1] ~ 1, Locs, newdata=grd, idp=2.0)

    r   <- raster::raster(P.idw)
    R_pts <- raster::rasterToPoints(r, spatial = TRUE)
    R_df <- data.frame(R_pts)
    colnames(R_df)[1] <- 'z'

    ras_map <-  ggplot2::ggplot() + ggplot2::geom_sf(data = world, fill = 'grey99') +
      ggplot2::geom_raster(data = R_df, ggplot2::aes(x = x, y= y, fill = z)) +
      ggplot2::scale_fill_gradientn(colors = col, name = statistic)
    for(i in 1:length(countries)) {
      ras_map <- ras_map + ggplot2::geom_sf(data = Country_borders[[i]], fill = "NA")
    }

    Div_inter_map <- ras_map +
      ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black") +
      ggplot2::xlab("Longitude") + ggplot2::ylab ("Latitude") +
      ggplot2::theme(panel.background = ggplot2::element_rect(color = "lightgray"), panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank())
    # Write out the interpolated raster for use in GIS software
    raster::crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    raster::writeRaster(r, filename = paste(prefix, '_', 'Interpolated', statistic, '.tiff', sep = ''), format = 'GTiff', overwrite = TRUE)

    Output <- list(Div_map, Div_inter_map, Div_mat)
    names(Output) <- c(paste(statistic," Map", sep = ""), paste("Interpolated ", statistic," Map", sep = ""), "Plotting Dataframe")
    return(Output)

  }
  else if(plot.type == 'point'){
    ### Heterozygosity Map
    # Map it with colored points
    Div_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = Long, y = Lat, color = Div_mat[,1]), shape = 19, size = 3) +
      ggplot2::scale_color_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = mean(Div_mat[,1]), breaks = Breaks, name = statistic) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") +
      ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    Output <- list(Div_map, Div_mat)
    names(Output) <- c(paste(statistic," Map", sep = ""), "Plotting Dataframe")
    return(Output)
  }
  else if(plot.type == 'interpolated'){

    sp::coordinates(Locs) <- c('Long', 'Lat')
    grd <- as.data.frame(sp::spsample(Locs,'regular', n=100000))
    names(grd) <- c("X", "Y")
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object
    P.idw <- gstat::idw(Div_mat[,1] ~ 1, Locs, newdata=grd, idp=2.0)

    r   <- raster::raster(P.idw)
    R_pts <- raster::rasterToPoints(r, spatial = TRUE)
    R_df <- data.frame(R_pts)
    colnames(R_df)[1] <- 'z'

    ras_map <-  ggplot2::ggplot() + ggplot2::geom_sf(data = world, fill = 'grey99') +
      ggplot2::geom_raster(data = R_df, ggplot2::aes(x = x, y= y, fill = z)) +
      ggplot2::scale_fill_gradientn(colors = col, name = statistic)
    for(i in 1:length(countries)) {
      ras_map <- ras_map + ggplot2::geom_sf(data = Country_borders[[i]], fill = "NA")
    }

    Div_inter_map <- ras_map +
      ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_point(data = Div_mat, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black") +
      ggplot2::xlab("Longitude") + ggplot2::ylab ("Latitude") +
      ggplot2::theme(panel.background = ggplot2::element_rect(color = "lightgray"), panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank())
    # Write out the interpolated raster for use in GIS software
    raster::crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    raster::writeRaster(r, filename = paste(prefix, '_', 'Interpolated', statistic, '.tiff', sep = ''), format = 'GTiff', overwrite = TRUE)

    Output <- list(Div_inter_map, Div_mat)
    names(Output) <- c(paste("Interpolated ", statistic," Map", sep = ""), "Plotting Dataframe")
    return(Output)
  }

}
