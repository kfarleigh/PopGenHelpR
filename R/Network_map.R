#' A function to map statistics (i.e., genetic differentiation) between points as a network on a map.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. If it is a csv, the 1st row should contain the individual/population names. The columns should also be named in this fashion.
#' @param pops Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The columns should be named Sample, containing the sample IDs; Population indicating the population assignment of the individual; Long, indicating the longitude of the sample; Lat, indicating the latitude of the sample. Alternatively, see the Longitude_col and Latitude_col arguments.
#' @param neighbors Numeric or character. The number of neighbors to plot connections with, or the specific relationship that you want to visualize. Names should match those in the population assignment file and be seperated by an underscore. If I want to visualize the relationship between East and West, for example, I would set neighbors = "East_West".
#' @param col Character vector indicating the colors you wish to use for plotting.
#' @param statistic Character indicating the statistic being plotted. This will be used to title the legend. The legend title will be blank if left as NULL.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. The number of breaks should match the number of colors.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#' @param Latitude_col Numeric. The number of the column indicating the latitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Lat column.
#' @param Longitude_col Numeric. The number of the column indicating the longitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Long column.
#' @param shapefile Character. A file name, vector of file names of a shapefile(s) to plot on the map, or a spatvector object that is compatible with the R package terra. This should be used in conjunction with the shapefile_plot_position argument.
#' @param raster Character.A file name or a spatraster object that is compatible with the terra R package. This should be used in conjunction with the raster_plot_position argument.
#' @param country_code Character. A country code or vector of country codes from the R package [geodata](https://cran.r-project.org/package=geodata) specifying the country that you want to plot administrative borders for (e.g, US states). You can determine the correct codes using geodata's `country_codes` function.
#' @param legend_pos Character. The desired position of the legend. The default is "none", which removes the legend. Other options include "left", "right", "top" or "bottom". Please see the ggplot2 documentation for all of the legend placement options.
#' @param scale_bar Boolean. Whether or not to add a scale bar. Note that maps with large areas or those that use unprojected spatial data (i.e., WGS 84) will generate a warning that the scale bar varies.
#' @param north_arrow Boolean. Whether or not to add a north arrow.
#' @param north_arrow_style Character. Which style of north arrow to add. See [ggspatial](https://cran.r-project.org/package=ggspatial) documentation for more details.
#' @param north_arrow_position Character. The position of the north arrow. See [ggspatial](https://cran.r-project.org/package=ggspatial) documentation for more details.
#' @param shapefile_plot_position Numeric. A number indicating which position to plot the shapefile in. The options are 1, which plots the shapefile on top of the base world map (under points and administrative boundaries), 2 which plots the shapefile on top of administrative boundaries (but under points), and 3, which plots the shapefile on top of everything.
#' @param raster_plot_position Numeric. A number indicating which position to plot the shapefile in. The options are 1, which plots the raster on top of the base world map (under points and administrative boundaries), 2 which plots the raster on top of administrative boundaries (but under points), and 3, which plots the raster on top of everything.
#' @param shapefile_col Character. A color or color vector indicating the color to fill the shapefile(s) with. Shapefiles will be colored alphabetically.
#' @param shapefile_outline_col Character. A color indicating the outline color of the shapefile.
#' @param shp_outwidth Numeric. The width of the shapefile outline.
#' @param raster_col Character. A character vector indicating the colors used to visualize the raster. The function will seperate your raster data into the same number of bins as there are colors. If you provide 5 colors, for example, there will be 5 bins.
#' @param interpolate_raster Boolean. Whether or not to interpolate the raster. The default is to interpolate the raster.
#' @param raster_breaks Numeric or Character vector. Values to be used as breaks for the raster surface.
#' @param discrete_raster Boolean. Indicating whether or not the raster being supplied is discrete.
#'
#' @return A list containing the map and the matrix used to plot the map.
#'
#' @author Keaka Farleigh
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(Fst_dat)
#' Fst <- Fst_dat[[1]]
#' Loc <- Fst_dat[[2]]
#' Test <- Network_map(dat = Fst, pops = Loc,
#' neighbors = 2,col = c('#4575b4', '#91bfdb', '#e0f3f8','#fd8d3c','#fc4e2a'),
#' statistic = "Fst", Lat_buffer = 1, Long_buffer = 1)}
Network_map <- function(dat, pops, neighbors, col, statistic = NULL, breaks = NULL, Lat_buffer = 1, Long_buffer = 1, Latitude_col = NULL, Longitude_col = NULL,
                        country_code = NULL, shapefile = NULL,
                        raster = NULL, legend_pos = "none", scale_bar = FALSE,
                        north_arrow = FALSE, north_arrow_style = ggspatial::north_arrow_nautical(),
                        north_arrow_position = NULL, shapefile_plot_position = NULL, raster_plot_position = NULL,
                        shapefile_col = NULL, shapefile_outline_col = NULL, shp_outwidth = 1,
                        raster_col = c('#2c7bb6','#abd9e9','#ffffbf','#fdae61', '#d7191c'),
                        interpolate_raster = NULL, raster_breaks = NULL, discrete_raster = NULL){
  X1 <- X2 <- X3 <- X4 <- Dif <- Long <- Lat <- alpha <- x <- y <- Outline_col <- value <- NULL

  # Read in files
  if(is.data.frame(dat) | is.matrix(dat)){
    Dif_mat <- dat
  }
  else if(is.character(dat) == TRUE){
    Dif_mat <- utils::read.csv(dat, row.names = 1)
  }
  else{
    stop("Please supply a dataframe or .csv file name for analysis")
  }
  if(is.data.frame(pops) == TRUE){
    Pops <- pops
  }
  else if(is.character(pops) == TRUE){
    Pops <- utils::read.csv(pops)
  }
  else{
    stop("Please supply a population assignment file as a dataframe or .csv file name")
  }

  # Make sure that the diagonal of the diferentiation matrix is 0
  diag(Dif_mat) <- 0


  if(!is.null(Latitude_col)){
    colnames(Pops)[Latitude_col] <- "Lat"
  }

  if(!is.null(Longitude_col)){
    colnames(Pops)[Longitude_col] <- "Long"
  }

  ################### Get the data for mapping
  ### Get coordinate ranges for our data
  Lat_Min <- min(Pops$Lat) - Lat_buffer
  Lat_Max <- max(Pops$Lat) + Lat_buffer
  Long_Min <- min(Pops$Long) - Long_buffer
  Long_Max <- max(Pops$Long) + Long_buffer

  # Set an extent to crop our data
  extent <- terra::ext((Long_Min-3), (Long_Max+3), (Lat_Min-3), (Lat_Max+3))

  # Get the world map
  world <- geodata::world(path = tempdir())

  # Crop the world to the extent
  world <- terra::crop(world, extent)

  # Get the country data if specified
  if(!is.null(country_code)){
    country_toplot <- geodata::gadm(country = country_code, path = tempdir())

    # Crop the country to extent
    country_toplot <- terra::crop(country_toplot, extent)

    # Convert to sf object for plotting
    world_sf <- sf::st_as_sf(world)
    countries_sf <- sf::st_as_sf(country_toplot)

    base_map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") +
      ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0))

  } else{

    # Convert to sf object for plotting
    world_sf <- sf::st_as_sf(world)

    base_map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4")
  }

  # Pull coordinates

  if(!is.null(Latitude_col) | !is.null(Longitude_col)){
    Coords <- Pops[,c(Longitude_col, Latitude_col)]
  } else{
    Coords <- Pops[,c(which(colnames(Pops) == "Long"), which(colnames(Pops) == "Lat"))]
  }

  # Get coordinates for each population
  Mapping <- Pops[!duplicated(Pops[,2]),]

  # Get the combination of all points
  All_pts <- data.frame(t(apply(utils::combn(seq_len(nrow(Mapping)), 2), 2, function(x) c(Mapping[x,]$Long, Mapping[x,]$Lat))))

  # Get the midpoint between each point pair
  Midpts<- data.frame(rowMeans(All_pts[,1:2]), rowMeans(All_pts[,3:4]))
  colnames(Midpts) <- c('x', 'y')

  # Get the names for each comparison
  Pop_combos <- utils::combn(Mapping$Population, 2, paste, collapse = "_")
  rownames(All_pts) <- Pop_combos
  Midpts$Combo <- Pop_combos

  # Reorganize diagonal matrix into long format
  Dif_Mapping <- data.frame(
    value=Dif_mat[lower.tri(Dif_mat,diag=T)],
    col=rep(1:ncol(Dif_mat),ncol(Dif_mat):1),
    row=unlist(sapply(1:nrow(Dif_mat),function(x) x:nrow(Dif_mat))))
  Dif_Mapping <- Dif_Mapping[-c(which(Dif_Mapping[,2] == Dif_Mapping[,3] | Dif_Mapping[,3] == Dif_Mapping[,2])),]
  Dif_Mapping$labels <- Pop_combos

  # Make sure that Dif_Mapping is in the same order as All_pts
  if(any(table(Dif_Mapping$labels == rownames(All_pts)) == FALSE)){
    stop("The matrices are not in the correct order, please contact the package authors if you see this error")
  }
  # Append the Dif statistic to the all pts matrix
  All_pts$Dif <- Dif_Mapping$value

  # Set colors
  if(missing(col)){
    col <- c('#4575b4', '#91bfdb', '#e0f3f8','#fd8d3c','#fc4e2a')
  } else{
    col <- col
  }

  # Set breaks
  if(is.null(breaks)){
    Breaks <- summary(All_pts$Dif)
    Breaks <- as.numeric(Breaks[1:5])
  }
  else{
    Breaks <- breaks
    Breaks <- as.numeric(Breaks)
  }
  # Get the neighbors
  if(missing(neighbors)){
    Coords <- as.matrix(Mapping[,3:4])
    Neigh <- spdep::knearneigh(Coords, k = 1, longlat = TRUE, use_kd_tree = FALSE)
    nNeigh <- as.data.frame(Neigh$nn)
    row.names(nNeigh) <- Mapping$Population
    nNeigh$Label <- Mapping$Population[nNeigh$V1]
    nNeigh$Comp1 <- paste(row.names(nNeigh), nNeigh$Label, sep = '_')
    nNeigh$Comp2 <- paste(nNeigh$Label, row.names(nNeigh), sep = '_')
    All_pts_toplot1 <- All_pts[which(rownames(All_pts) %in% nNeigh$Comp1),]
    All_pts_toplot2 <- All_pts[which(rownames(All_pts) %in% nNeigh$Comp2),]
    All_pts_toplot <- rbind(All_pts_toplot1, All_pts_toplot2) }
  else if(is.numeric(neighbors)){
    Coords <- as.matrix(Mapping[,3:4])
    Neigh <- spdep::knearneigh(Coords, k = neighbors, longlat = TRUE, use_kd_tree = FALSE)
    nNeigh <- as.data.frame(Neigh$nn)
    row.names(nNeigh) <- Mapping$Population
    combs1 <- c()
    combs2 <- c()
    for (i in seq(1:ncol(nNeigh))){
      nNeigh[,neighbors+i] <- Mapping$Population[nNeigh[,i]]
      combs1 <- c(combs1, paste(row.names(nNeigh), nNeigh[,neighbors+i], sep = '_'))
      combs2 <- c(combs2, paste(nNeigh[,neighbors+i], row.names(nNeigh), sep = '_'))}

    combs_tot <- c(combs1,combs2)
    combs_tot <- unique(combs_tot)
    All_pts_toplot <- All_pts[which(rownames(All_pts) %in% combs_tot),]
    All_pts_toplot$Comp <- rownames(All_pts_toplot)
  }
  else if(is.character(neighbors)){
    combs1 <- neighbors
    combs2 <- c()
    for (i in seq(1:length(neighbors))) {
      nsplit <- strsplit(neighbors[i], "_")
      combs2 <- c(combs2, paste(rev(nsplit[[1]]), collapse = "_"))
    }
    combs_tot <- c(combs1,combs2)
    All_pts_toplot <- All_pts[which(row.names(All_pts) %in% combs_tot),]
    All_pts_toplot$Comp <- rownames(All_pts_toplot)
    All_pts_toplot <- All_pts_toplot[!duplicated(All_pts_toplot$Comp),]
  }
  Breaks <- round(Breaks,2)
  All_pts_toplot$Dif <- round(All_pts_toplot$Dif,2)


  # Map it
  if(is.null(statistic)){
    map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4), color = 'black', linewidth= 1.05) +
      ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4, color = Dif), linewidth = 1) +
      ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = ggplot2::element_blank()) +
      ggplot2::geom_point(data = Mapping, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black") +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')
  } else{
    map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4), color = 'black', linewidth= 1.05) +
      ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4, color = Dif), linewidth = 1) +
      ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
      ggplot2::geom_point(data = Mapping, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black") +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')
  }

  # Create a layer for plotting with shapefiles and/or rasters

  net_layer_out <- ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4), color = 'black', linewidth= 1.05)
  net_layer_col <- ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4, color = Dif), linewidth = 1)
  net_layer_pt <- ggplot2::geom_point(data = Mapping, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black")

  ### Shapefile handling
  ### If there is or isn't a shapefile

  if(!is.null(shapefile)){
    # Get a list of the file names, create a list to store them
    shapefiles <- shapefile
    all_shp_mapped <- list()
    if(!is.null(shapefile_col)){
      # Get the shapefile colors if they are provided
      shp_fill <- shapefile_col
    } else{
      # Make transparent if there is no shapefile_col provided
      shp_fill <- rep(NA, length(shapefile))
    }

    if(!is.null(shapefile_outline_col)){
      # Set the shapefile outline color if one is provided
      shp_out <- shapefile_outline_col
    } else{
      # Otherwise, make it the Outline_col (probably black)
      shp_out <- rep(Outline_col, length(shapefiles))
    }

    # Read in all of the shapefiles as Spatvectors, convert to sf object, and plot.
    for(i in 1:length(shapefile)){
      # Read in as SpatVector
      if(is.character(shapefile[i])){
        tmp_vect <- terra::vect(shapefile[i])
      } else{
        tmp_vect <- shapefile[i]
      }
      # Convert to sf object
      tmp_sf <- sf::st_as_sf(tmp_vect)

      # Plot
      tmp_map <- ggplot2::geom_sf(data = tmp_sf, fill = ggplot2::alpha(shp_fill[i], 1), color = shp_out[i], linewidth = shp_outwidth)

      # Store as a list element
      all_shp_mapped[[i]] <- tmp_map

      # Remove tmp objects
      remove(tmp_vect, tmp_sf, tmp_map)
    }

    ### Plot according to the shapefile_plot_position argument
    # The first if and else if are required because the position of the shapefile depends on whether or not the users wants administrative boundaries.
    if(shapefile_plot_position == 1 && !is.null(country_code)){
      map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + all_shp_mapped +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (shapefile_plot_position == 1 && is.null(country_code)){
      map <- base_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 2){
      map <- base_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 3){
      map <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        all_shp_mapped +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    }
  }

  ### Raster handling

  ################### Adding a raster to the map
  if(!is.null(raster)){
    # Get a list of the file names, create a list to store them
    rasters <- raster

    if(is.character(raster)){
      raster <- terra::rast(raster)
    } else{
      raster <- raster
    }

    # The lines below are code to project the raster, I am wary to do this because it could take forever/crash users computers
    # crs <- "+proj=longlat +datum=WGS84"

    # raster_proj <- terra::project(raster, crs)

    # We will clip the raster to our extent so that it is smaller and doesn't crash computers/illustrator of advanced users
    raster_mask <- terra::mask(raster, mask = extent)

    # Create a data frame
    rast_df <- terra::as.data.frame(raster_mask, xy=TRUE)

    # Rename the 3rd column so the plotting works
    colnames(rast_df)[3] <- "value"

    # Set the colors for plotting
    if(!is.null(raster_col)){
      # Get the raster colors if they are provided
      ras_fill <- raster_col
    } else{
      # Use a default palette
      ras_fill <- c('#2c7bb6','#abd9e9','#ffffbf','#fdae61', '#d7191c')
    }

    # Shold raster be interpolated?
    if(is.null(interpolate_raster)){
      interpolate <- TRUE
    } else if(isTRUE(interpolate_raster)){
      interpolate <- TRUE
    } else{
      interpolate <- FALSE
    }

    # Get the raster breaks, or set our own breaks using the summary function
    if(!is.null(raster_breaks)){
      ras_breaks <- as.numeric(raster_breaks)

      # Check to make sure that the raster_breaks is equal to the ras_fill length, throw and error if not.
      if(length(ras_breaks) != length(ras_fill)){
        stop("The number of breaks does not match the number of colors.
             Please make sure that the raster breaks argument has the same amount of breaks as colors if you supplied them.
             Alternatively, there should be 5 breaks to use the default color palette.")
      }
    } else{
      ras_breaks <- c(summary(rast_df[,3])[1],summary(rast_df[,3])[2], summary(rast_df[,3])[4], summary(rast_df[,3])[5], summary(rast_df[,3])[6])

    }

    if(isTRUE(discrete_raster)){
      # Reclassify the raster using the ras_breaks object
      rast_rc <- terra::classify(raster_mask, ras_breaks, right = TRUE)

      rast_df <- terra::as.data.frame(rast_rc, xy = TRUE)
      colnames(rast_df)[3] <- "value"
      rast_df$value <- as.factor(rast_df$value)
    }

    # Create a layer for the raster
    raster_map <- ggplot2::geom_raster(data = rast_df, ggplot2::aes(x = x, y = y, fill = value), interpolate = interpolate)

    ### Plot according to the the raster_plot_position argument
    # The first if and else if are required because the position of the shapefile depends on whether or not the users wants administrative boundaries.
    if(raster_plot_position == 1 && !is.null(country_code) && isTRUE(discrete_raster)){

      map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") +
        raster_map +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(raster_plot_position == 1 && !is.null(country_code) && !isTRUE(discrete_raster)) {

      map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") +
        raster_map +
        ggplot2::scale_fill_gradientn(breaks = as.numeric(ras_breaks), colors = ras_fill) +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (raster_plot_position == 1 && is.null(country_code) && isTRUE(discrete_raster)){

      map <- base_map + raster_map +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (raster_plot_position == 1 && is.null(country_code) && !isTRUE(discrete_raster)){

      map <- base_map + raster_map +
        ggplot2::scale_fill_gradientn(breaks = as.numeric(ras_breaks), colors = ras_fill) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(raster_plot_position == 2 && isTRUE(discrete_raster)){
      map <- base_map +
        raster_map +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(raster_plot_position == 2 && !isTRUE(discrete_raster)){
      map <- base_map +
        raster_map +
        ggplot2::scale_fill_gradientn(breaks = as.numeric(ras_breaks), colors = ras_fill) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(raster_plot_position == 3 && isTRUE(discrete_raster)){
      map <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        raster_map +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(raster_plot_position == 3 && !isTRUE(discrete_raster)){
      map <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        raster_map +
        ggplot2::scale_fill_gradientn(breaks = as.numeric(ras_breaks), colors = ras_fill) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    }
  }

  # Create a map if both a shapefile and raster are provided. We will overlay the shapefile onto the raster and only show the outline.
  all_shp_mapped <- list()
  if(!is.null(shapefile) && !is.null(raster)){

    # Recreate the all_shp_mapped object, but make everything transparent
    for(i in 1:length(shapefile)){
      # Read in as SpatVector
      if(is.character(shapefile[i])){
        tmp_vect <- terra::vect(shapefile[i])
      } else{
        tmp_vect <- shapefile[i]
      }
      # Convert to sf object
      tmp_sf <- sf::st_as_sf(tmp_vect)

      # Plot
      tmp_map <- ggplot2::geom_sf(data = tmp_sf, fill = ggplot2::alpha(shp_fill[i], 0), color = shp_out[i], linewidth = shp_outwidth)

      # Store as a list element
      all_shp_mapped[[i]] <- tmp_map

      # Remove tmp objects
      remove(tmp_vect, tmp_sf, tmp_map)
    }

    #### Add in relation to base map
    if(shapefile_plot_position == 1 && !is.null(country_code) && isTRUE(discrete_raster)){
      map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + raster_map + all_shp_mapped +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 1 && !is.null(country_code) && !isTRUE(discrete_raster)){

      map <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + raster_map + all_shp_mapped +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::scale_fill_gradientn(colors = ras_fill, breaks = as.numeric(ras_breaks)) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (shapefile_plot_position == 1 && is.null(country_code) && isTRUE(discrete_raster)){
      map <- base_map + raster_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (shapefile_plot_position == 1 && is.null(country_code) && !isTRUE(discrete_raster)){
      map <- base_map + raster_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::scale_fill_gradientn(colors = ras_fill, breaks = as.numeric(ras_breaks)) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 2 && isTRUE(discrete_raster)){
      map <- base_map + raster_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 2 && !isTRUE(discrete_raster)){
      map <- base_map + raster_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        ggplot2::scale_fill_gradientn(colors = ras_fill, breaks = as.numeric(ras_breaks)) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 3 && isTRUE(discrete_raster)){
      map <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        raster_map + all_shp_mapped +
        ggplot2::scale_fill_manual(values = ras_fill) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 3 && !isTRUE(discrete_raster)){
      map <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        net_layer_out + net_layer_col + net_layer_pt +
        ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = statistic) +
        raster_map + all_shp_mapped +
        ggplot2::scale_fill_gradientn(colors = ras_fill, breaks = as.numeric(ras_breaks)) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    }
  }

  ################### Add north arrow and scale bar
  if(isTRUE(scale_bar)){
    map <- map + ggspatial::annotation_scale()
  } else{
    map <- map
  }

  if(isTRUE(north_arrow)){
    if(is.null(north_arrow_position)){
      map <- map + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
    } else{
      map <- map + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
    }
  } else{
    map <- map
  }


  Output_difmap <- list(All_pts_toplot[,c("Dif","Comp")], map)
  names(Output_difmap) <- c("Statistic Matrix", "Map")
  return(Output_difmap)
}
