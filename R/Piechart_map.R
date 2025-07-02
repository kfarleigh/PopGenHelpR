#' Plot a map of ancestry pie charts.
#'
#' @param anc.mat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the names of each sample/population, followed by the estimated contribution of each cluster to that individual/pop.
#' @param pops Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The columns should be named Sample, containing the sample IDs; Population indicating the population assignment of the individual, population and sample names must be the same type (i.e., both numeric or both characters); Long, indicating the longitude of the sample; Lat, indicating the latitude of the sample. Alternatively, see the Longitude_col and Latitude_col arguments.
#' @param K Numeric.The number of genetic clusters in your data set, please contact the package authors if you need help doing this.
#' @param plot.type Character string. Options are all, individual, and population. All is default and recommended, this will plot a piechart map for both the individuals and populations.
#' @param col Character vector indicating the colors you wish to use for plotting.
#' @param piesize Numeric. The radius of the pie chart for ancestry mapping.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#' @param Latitude_col Numeric. The number of the column indicating the latitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Lat column.
#' @param Longitude_col Numeric. The number of the column indicating the longitude for each sample. If this is not null, PopGenHelpR will use this column instead of looking for the Long column.
#' @param shapefile Character. A file name, vector of file names of a shapefile(s) to plot on the map, or a spatvector object that is compatible with the R package terra. This should be used in conjunction with the shapefile_plot_position argument.
#' @param country_code Character. A country code or vector of country codes from the R package [geodata](https://cran.r-project.org/package=geodata) specifying the country that you want to plot administrative borders for (e.g, US states). You can determine the correct codes using geodata's `country_codes` function.
#' @param legend_pos Character. The desired position of the legend. The default is "none", which removes the legend. Other options include "left", "right", "top" or "bottom". Please see the ggplot2 documentation for all of the legend placement options.
#' @param scale_bar Boolean. Whether or not to add a scale bar. Note that maps with large areas or those that use unprojected spatial data (i.e., WGS 84) will generate a warning that the scale bar varies.
#' @param north_arrow Boolean. Whether or not to add a north arrow.
#' @param north_arrow_style Character. Which style of north arrow to add. See [ggspatial](https://cran.r-project.org/package=ggspatial) documentation for more details.
#' @param north_arrow_position Character. The position of the north arrow. See [ggspatial](https://cran.r-project.org/package=ggspatial) documentation for more details.
#' @param shapefile_plot_position Numeric. A number indicating which position to plot the shapefile in. The options are 1, which plots the shapefile on top of the base world map (under points and administrative boundaries), 2 which plots the shapefile on top of administrative boundaries (but under points), and 3, which plots the shapefile on top of everything.
#' @param shapefile_col Character. A color or color vector indicating the color to fill the shapefile(s) with. Similar to `group_col`, shapefiles will be colored alphabetically.
#' @param shapefile_outline_col Character. A color indicating the outline color of the shapefile.
#' @param shp_outwidth Numeric. The width of the shapefile outline.
#' @return A list containing your plots and the data frames used to generate the plots.
#' @importFrom magrittr %>%
#'
#' @author Keaka Farleigh
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(Q_dat)
#' Qmat <- Q_dat[[1]]
#' rownames(Qmat) <- Qmat[,1]
#' Loc <- Q_dat[[2]]
#' Test_all <- Piechart_map(anc.mat = Qmat, pops = Loc, K = 5,
#' plot.type = 'all', col = c('#d73027', '#fc8d59', '#e0f3f8', '#91bfdb', '#4575b4'), piesize = 0.35,
#' Lat_buffer = 1, Long_buffer = 1)}
Piechart_map <- function(anc.mat, pops, K, plot.type = 'all', col, piesize = 0.35, Lat_buffer, Long_buffer,
                         Latitude_col = NULL, Longitude_col = NULL, country_code = NULL, shapefile = NULL,
                         legend_pos = "none", scale_bar = FALSE,
                         north_arrow = FALSE, north_arrow_style = ggspatial::north_arrow_nautical(),
                         north_arrow_position = NULL, shapefile_plot_position = NULL,
                         shapefile_col = NULL, shapefile_outline_col = NULL, shp_outwidth = 1){

  Pop <- coeff <- Sample <- value <- variable <- aes <- Long <- Lat <- alpha <- ID<- NULL
  # Read in ancestry matrix and pop file
  if(missing(anc.mat)){
    stop("Please supply an ancestry matrix file for plotting, if you have questions
         on how to generate an ancestry matrix please email the package author")
  }
  else if(is.data.frame(anc.mat) == TRUE){
    Ind_anc <- anc.mat
  }
  else if(is.character(anc.mat) == TRUE){
    Ind_anc <- utils::read.csv(anc.mat)
  }
  if(is.data.frame(pops) == TRUE){
    Pops <- pops
  }
  else if(is.character(pops) == TRUE){
    Pops <- utils::read.csv(pops)
  }
  if(missing(col)){
    stop("Please supply a vector of colors for plotting")
  }
  else if(length(col) < K){
    stop("Please supply at least as many colors as your K value")
  }
  else{
    col <- col
  }

  if(!is.null(Latitude_col)){
    colnames(Pops)[Latitude_col] <- "Lat"
  }

  if(!is.null(Longitude_col)){
    colnames(Pops)[Longitude_col] <- "Long"
  }

  # Pull coordinates

  if(!is.null(Latitude_col) | !is.null(Longitude_col)){
    Coords <- Pops[,c(Longitude_col, Latitude_col)]
  } else{
    Coords <- Pops[,3:4]
  }

  # Set population data names
  colnames(Pops)[1:2] <- c("Sample", "Population")

  # Convert the individual and population information to characters for plotting

  Ind_anc[,1] <- as.character(Ind_anc[,1])
  Pops[,1] <- as.character(Pops[,1])
  Pops[,2] <- as.character(Pops[,2])

  ### Get coordinate ranges for our data
  Lat_Min <- min(Coords$Lat) - Lat_buffer
  Lat_Max <- max(Coords$Lat) + Lat_buffer
  Long_Min <- min(Coords$Long) - Long_buffer
  Long_Max <- max(Coords$Long) + Long_buffer

  ################### Get the data for mapping      !!! This needs to be fixed
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


  ######################## Piechart layers          !!! Need to figure out how to just make piechart layers !!!

  if(plot.type == 'all') {
      # Individual plots
      colnames(Ind_anc) <- c("Sample", paste0(rep("cluster", K), 1:K))
      Ind_anc2 <- Ind_anc
      ord <- Ind_anc2[,1]
      Ind_anc2$ID <- Ind_anc2[,1]
      Ind_anc2$ID <- factor(Ind_anc2$ID, levels = ord)



      # Population plots
      Pop_anc <- Ind_anc[,-1]
      Pop_anc$Pop <- Pops$Population
      Pop_anc[,c(ncol(Pop_anc) +1,ncol(Pop_anc) +2)] <- Pops[,3:4]
      Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise(dplyr::across(1:(K+2), mean, na.rm = TRUE)) %>%
        dplyr::ungroup()
      Pop_anc_coeff <- Pop_anc[,c(1:(K+1))]
      ord_pop <- Pop_anc_coeff$Pop
      Pop_anc_coeff$ID <- ord_pop
      Pop_anc_coeff$ID <- factor(Pop_anc_coeff$ID, levels = ord_pop)


      ## Create the layer for individuals

      # Add coordinates to the individual ancestry data frame
      Ind_anc[,c(ncol(Ind_anc) +1,ncol(Ind_anc) +2)] <- Pops[,3:4]


      Ind_map <- scatterpie::geom_scatterpie(data = Ind_anc, ggplot2::aes(Long, Lat, r = piesize), cols = c(colnames(Ind_anc[2:(K+1)])))

      # Map populations
      Pop_map <- scatterpie::geom_scatterpie(data = Pop_anc, ggplot2::aes(Long, Lat, r = piesize), cols = c(colnames(Pop_anc[2:(K+1)])))


    } else if(plot.type == 'individual'){
    # Individual plots
    colnames(Ind_anc) <- c("Sample", paste0(rep("cluster", K), 1:K))
    Ind_anc2 <- Ind_anc
    ord <- Ind_anc2[,1]
    Ind_anc2$ID <- Ind_anc2[,1]
    Ind_anc2$ID <- factor(Ind_anc2$ID, levels = ord)


    # Add coordinates to the individual ancestry data frame
    Ind_anc[,c(ncol(Ind_anc) +1,ncol(Ind_anc) +2)] <- Pops[,3:4]

    Ind_map <- scatterpie::geom_scatterpie(data = Ind_anc, ggplot2::aes(Long, Lat, r = piesize), cols = c(colnames(Ind_anc[2:(K+1)])))

  } else if(plot.type == 'population'){

    colnames(Ind_anc) <- c("Sample", paste0(rep("cluster", K), 1:K))

    # Population plots
    Pop_anc <- Ind_anc[,-1]
    Pop_anc$Pop <- Pops$Population
    Pop_anc[,c(ncol(Pop_anc) +1,ncol(Pop_anc) +2)] <- Pops[,3:4]
    Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise(dplyr::across(1:(K+2), mean, na.rm = TRUE)) %>%
      dplyr::ungroup()
    Pop_anc_coeff <- Pop_anc[,c(1:(K+1))]
    ord_pop <- Pop_anc_coeff$Pop
    Pop_anc_coeff$ID <- ord_pop
    Pop_anc_coeff$ID <- factor(Pop_anc_coeff$ID, levels = ord_pop)

    # Map populations
    Pop_map <- scatterpie::geom_scatterpie(data = Pop_anc, ggplot2::aes(Long, Lat, r = piesize), cols = c(colnames(Pop_anc[2:(K+1)])))
  }

  else {
    stop("Please supply input for plot.type. The options are 'all', 'individual', or 'population'.")

  }


  # Generate your plots if there is no shapefile or raster
  if(is.null(shapefile)){
    if(plot.type == "all"){

      # Plot individuals
      Ind_piemap <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) + Ind_map +
        ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                   labels = c(paste('Cluster', 1:K, sep = ' ')),
                                   values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      # Plot populations
      Pop_piemap <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) + Pop_map +
        ggplot2::scale_fill_manual(breaks = c(colnames(Pop_anc[2:(K+1)])),
                                   labels = c(paste('Cluster', 1:K, sep = ' ')),
                                   values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      ### Add north arrow, scale bar if necessary
      ################### Add north arrow and scale bar
      if(isTRUE(scale_bar)){
        Ind_piemap <- Ind_piemap + ggspatial::annotation_scale()
        Pop_piemap <- Pop_piemap + ggspatial::annotation_scale()
      } else{
        Ind_piemap <- Ind_piemap
        Pop_piemap <- Pop_piemap
      }

      if(isTRUE(north_arrow)){
        if(is.null(north_arrow_position)){
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
        } else{
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
        }
      } else{
        Ind_piemap <- Ind_piemap
        Pop_piemap <- Pop_piemap
      }

      # Output map, mapping data
      out <- list(Ind_piemap, Pop_piemap, Ind_anc, Pop_anc)
      names(out) <- c("Individual_piemap", "Population_piemap", "Individal_piemap_data", "Population_piemap_data")


    } else if(plot.type == 'individual') {
      # Plot individuals
      Ind_piemap <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) + Ind_map +
        ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                   labels = c(paste('Cluster', 1:K, sep = ' ')),
                                   values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')


      ### Add north arrow, scale bar if necessary
      ################### Add north arrow and scale bar
      if(isTRUE(scale_bar)){
        Ind_piemap <- Ind_piemap + ggspatial::annotation_scale()
      } else{
        Ind_piemap <- Ind_piemap
      }

      if(isTRUE(north_arrow)){
        if(is.null(north_arrow_position)){
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
        } else{
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
        }
      } else{
        Ind_piemap <- Ind_piemap
      }

      # Output map, mapping data
      out <- list(Ind_piemap, Ind_anc)
      names(out) <- c("Individual_piemap", "Individal_piemap_data")

    } else if(plot.type == 'population'){

      # Plot populations
      Pop_piemap <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) + Pop_map +
        ggplot2::scale_fill_manual(breaks = c(colnames(Pop_anc[2:(K+1)])),
                                   labels = c(paste('Cluster', 1:K, sep = ' ')),
                                   values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      ### Add north arrow, scale bar if necessary
      ################### Add north arrow and scale bar
      if(isTRUE(scale_bar)){
        Pop_piemap <- Pop_piemap + ggspatial::annotation_scale()
      } else{
        Pop_piemap <- Pop_piemap
      }

      if(isTRUE(north_arrow)){
        if(is.null(north_arrow_position)){
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
        } else{
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
        }
      } else{
        Pop_piemap <- Pop_piemap
      }
      # Output map, mapping data
      out <- list(Pop_piemap, Pop_anc)
      names(out) <- c("Population_piemap", "Population_piemap_data")
    }

    return(out)
  }


  ### Shapefile handling, plotting if there is only a shapefile and no raster !!! Need to address plotting !!!

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
      # Otherwise, make it the Outline_col transparent
      shp_out <- rep(NA, length(shapefiles))
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

    if(!is.null(shapefile) && plot.type == "all"){

      if(shapefile_plot_position == 1 && !is.null(country_code)){
        Ind_piemap <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + all_shp_mapped +
          ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                     labels = c(paste('Cluster', 1:K, sep = ' ')),
                                     values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

        Pop_piemap <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + all_shp_mapped +
          ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      } else if (shapefile_plot_position == 1 && is.null(country_code)){
        Ind_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

        Pop_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      } else if(shapefile_plot_position == 2){
        Ind_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

        Pop_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      } else if(shapefile_plot_position == 3){
        Ind_piemap <- base_map +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + all_shp_mapped +
          ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                     labels = c(paste('Cluster', 1:K, sep = ' ')),
                                     values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

        Pop_piemap <- base_map +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Pop_map + all_shp_mapped +
          ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                     labels = c(paste('Cluster', 1:K, sep = ' ')),
                                     values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

      }

      ### Add north arrow, scale bar if necessary
      ################### Add north arrow and scale bar
      if(isTRUE(scale_bar)){
        Ind_piemap <- Ind_piemap + ggspatial::annotation_scale()
        Pop_piemap <- Pop_piemap + ggspatial::annotation_scale()
      } else{
        Ind_piemap <- Ind_piemap
        Pop_piemap <- Pop_piemap
      }

      if(isTRUE(north_arrow)){
        if(is.null(north_arrow_position)){
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)
        } else{
          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
          Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)
        }
      } else{
        Ind_piemap <- Ind_piemap
        Pop_piemap <- Pop_piemap
      }

      # Output map, mapping data
      out <- list(Ind_piemap, Pop_piemap, Ind_anc, Pop_anc)
      names(out) <- c("Individual_piemap", "Population_piemap", "Individal_piemap_data", "Population_piemap_data")

    } else if(!is.null(shapefile) && plot.type == 'individual'){
      if(shapefile_plot_position == 1 && !is.null(country_code)){

        Ind_piemap <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + all_shp_mapped +
          ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')



      } else if (shapefile_plot_position == 1 && is.null(country_code)){

        Ind_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')


      } else if(shapefile_plot_position == 2){

        Ind_piemap <- base_map + all_shp_mapped +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                               labels = c(paste('Cluster', 1:K, sep = ' ')),
                                               values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')



      } else if(shapefile_plot_position == 3){

        Ind_piemap <- base_map +
          ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
          Ind_map + all_shp_mapped +
          ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                     labels = c(paste('Cluster', 1:K, sep = ' ')),
                                     values = col) +
          ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')


      }

      ### Add north arrow, scale bar if necessary
      ################### Add north arrow and scale bar
      if(isTRUE(scale_bar)){
        Ind_piemap <- Ind_piemap + ggspatial::annotation_scale()

      } else{
        Ind_piemap <- Ind_piemap

      }

      if(isTRUE(north_arrow)){
        if(is.null(north_arrow_position)){

          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)


        } else{

          Ind_piemap <- Ind_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)

        }
      } else{

        Ind_piemap <- Ind_piemap
      }

      # Output map, mapping data
      out <- list(Ind_piemap, Ind_anc)
      names(out) <- c("Individual_piemap", "Individal_piemap_data")

  } else if(!is.null(shapefile) && plot.type == 'population'){

    if(shapefile_plot_position == 1 && !is.null(country_code)){

      Pop_piemap <- ggplot2::ggplot() + ggplot2::geom_sf(data = world_sf, fill = "#f4f4f4") + all_shp_mapped +
        ggplot2::geom_sf(data = countries_sf, fill = ggplot2::alpha("#f4f4f4", 0)) +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                             labels = c(paste('Cluster', 1:K, sep = ' ')),
                                             values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if (shapefile_plot_position == 1 && is.null(country_code)){

      Pop_piemap <- base_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                             labels = c(paste('Cluster', 1:K, sep = ' ')),
                                             values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 2){

      Pop_piemap <- base_map + all_shp_mapped +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        Pop_map + ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                             labels = c(paste('Cluster', 1:K, sep = ' ')),
                                             values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    } else if(shapefile_plot_position == 3){


      Pop_piemap <- base_map +
        ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
        Pop_map + all_shp_mapped +
        ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                   labels = c(paste('Cluster', 1:K, sep = ' ')),
                                   values = col) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = legend_pos) + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    }

    ### Add north arrow, scale bar if necessary
    ################### Add north arrow and scale bar
    if(isTRUE(scale_bar)){

      Pop_piemap <- Pop_piemap + ggspatial::annotation_scale()

    } else{

      Pop_piemap <- Pop_piemap
    }

    if(isTRUE(north_arrow)){
      if(is.null(north_arrow_position)){

        Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = "tl", style = north_arrow_style, which_north = TRUE)

      } else{

        Pop_piemap <- Pop_piemap + ggspatial::annotation_north_arrow(location = north_arrow_position, style = north_arrow_style, which_north = TRUE)

      }
    } else{

      Pop_piemap <- Pop_piemap

    }

    # Output map, mapping data
    out <- list(Pop_piemap, Pop_anc)
    names(out) <- c("Population_piemap", "Population_piemap_data")

  }

    return(out)

  }


}
