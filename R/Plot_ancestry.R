#' Plot an ancestry matrix and map of ancestry pie charts.
#'
#' @param anc.mat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the names of each sample/population, followed by the estimated contribution of each cluster to that individual/pop.
#' @param pops Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The columns should be named Sample, containing the sample IDs; Population indicating the population assignment of the individual; Long, indicating the longitude of the sample; Lat, indicating the latitude of the sample.
#' @param K Numeric.The number of genetic clusters in your data set, please contact the package authors if you need help doing this.
#' @param plot.type Character string. Options are all, individual, and population. All is default and recommended, this will plot a barchart and piechart map for both the individuals and populations.
#' @param col Character vector indicating the colors you wish to use for plotting.
#' @param countries Character vector indicating the country borders that you wish to plot on a map. Can be any country that is valid in the ne_statest function in the rnaturalearth package.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#'
#' @return A list containing your plots and the data frames used to generate the plots.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' data(Q_dat)
#' Qmat <- Q_dat[[1]]
#' Loc <- Q_dat[[2]]
#' Test_all <- Plot_ancestry(anc.mat = Qmat, pops = Loc, K = 5,
#' plot.type = 'all', col = c('red', 'orange', 'pink', 'purple', 'blue'),
#' countries = c("united states of america", "mexico"), Lat_buffer = 1, Long_buffer = 1)}
Plot_ancestry <- function(anc.mat, pops, K, plot.type = 'all', col, countries, Lat_buffer, Long_buffer){
  Pop <- aes <- Long <- Lat <- NULL
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

  # Pull coordinates
  Coords <- Pops[,3:4]

  ################### Get the data for mapping
  # Get map data
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

  ### Get coordinate ranges for our data
  Lat_Min <- min(Coords$Lat) - Lat_buffer
  Lat_Max <- max(Coords$Lat) + Lat_buffer
  Long_Min <- min(Coords$Long) - Long_buffer
  Long_Max <- max(Coords$Long) + Long_buffer

  ######################## Barcharts

  if(plot.type == 'all') {
    # Individual plots
    graphics::barplot(t(as.matrix(Ind_anc[,2:ncol(Ind_anc)])), col=col,
            space = 0, xlab="Individual", ylab = "Ancestry proportions",
            border=NA, axisnames = FALSE)
    graphics::axis(1, at = 1:nrow(Ind_anc),labels = Pops$Sample, las=2 ,cex.axis = .4)
    # Population plots
    Pop_anc <- Ind_anc[,-1]
    Pop_anc$Pop <- Pops$Population
    Pop_anc[,c(ncol(Pop_anc) +1,ncol(Pop_anc) +2)] <- Pops[,3:4]
    Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise(dplyr::across(,mean, na.rm = TRUE))

    graphics::barplot(t(as.matrix(Pop_anc[,2:(K+1)])), col=col,
            space = 0, xlab="Population", ylab = "Ancestry proportions",
            border=NA, axisnames = FALSE)
    graphics::axis(1, at = 1:nrow(Pop_anc),labels = Pop_anc$Pop, las=2 ,cex.axis = .4)

    # Map individuals

    # Add coordinates to the individual ancestry data frame
    Ind_anc[,c(ncol(Ind_anc) +1,ncol(Ind_anc) +2)] <- Pops[,3:4]


    Ind_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      scatterpie::geom_scatterpie(data = Ind_anc, ggplot2::aes(Long, Lat, r = 0.35), cols = c(colnames(Ind_anc[2:(K+1)]))) +
      ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                 labels = c(paste('Cluster', 1:K, sep = ' ')),
                                 values = c(col[1], col[2], col[3], col[4], col[5], col[6], col[7], col[8], col[9], col[10])) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "none") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    # Map populations
    Pop_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      scatterpie::geom_scatterpie(data = Pop_anc, ggplot2::aes(Long, Lat, r = 0.35), cols = c(colnames(Pop_anc[2:(K+1)]))) +
      ggplot2::scale_fill_manual(breaks = c(colnames(Pop_anc[2:(K+1)])),
                                 labels = c(paste('Cluster', 1:K, sep = ' ')),
                                 values = c(col[1], col[2], col[3], col[4], col[5], col[6], col[7], col[8], col[9], col[10])) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "none") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')
    Output_all <- list(Ind_anc, Pop_anc, Ind_map, Pop_map)
    names(Output_all) <- c("Individual Ancestry Matrix", "Population Ancestry Matrix", "Individual Map", "Population Map")
    return(Output_all)
  }
  else if(plot.type == 'individual'){
    graphics::barplot(t(as.matrix(Ind_anc[,2:ncol(Ind_anc)])), col=col,
            space = 0, xlab="Individual", ylab = "Ancestry proportions",
            border=NA, axisnames = FALSE)
    graphics::axis(1, at = 1:nrow(Ind_anc),labels = Pops$Sample, las=2 ,cex.axis = .4)

    # Add coordinates to the individual ancestry data frame
    Ind_anc[,c(ncol(Ind_anc) +1,ncol(Ind_anc) +2)] <- Pops[,3:4]

    Ind_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      scatterpie::geom_scatterpie(data = Ind_anc, ggplot2::aes(Long, Lat, r = 0.35), cols = c(colnames(Ind_anc[2:(K+1)]))) +
      ggplot2::scale_fill_manual(breaks = c(colnames(Ind_anc[2:(K+1)])),
                                 labels = c(paste('Cluster', 1:K, sep = ' ')),
                                 values = c(col[1], col[2], col[3], col[4], col[5], col[6], col[7], col[8], col[9], col[10])) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "none") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')
    Output_indanc <- list(Ind_anc, Ind_map)
    names(Output_indanc) <- c("Individual Ancestry Matrix", "Individual Piechart Map")
    return(Output_indanc)
  }

  else if(plot.type == 'population'){
    Pop_anc <- Ind_anc[,-1]
    Pop_anc$Pop <- Pops$Population
    Pop_anc[,c(ncol(Pop_anc) +1,ncol(Pop_anc) +2)] <- Pops[,3:4]
    Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise(dplyr::across(,mean, na.rm = TRUE))

    graphics::barplot(t(as.matrix(Pop_anc[,2:(K+1)])), col=col,
            space = 0, xlab="Population", ylab = "Ancestry proportions",
            border=NA, axisnames = FALSE)
    graphics::axis(1, at = 1:nrow(Pop_anc),labels = Pop_anc$Pop, las=2 ,cex.axis = .4)

    Pop_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
      scatterpie::geom_scatterpie(data = Pop_anc, ggplot2::aes(Long, Lat, r = 0.35), cols = c(colnames(Pop_anc[2:(K+1)]))) +
      ggplot2::scale_fill_manual(breaks = c(colnames(Pop_anc[2:(K+1)])),
                                 labels = c(paste('Cluster', 1:K, sep = ' ')),
                                 values = c(col[1], col[2], col[3], col[4], col[5], col[6], col[7], col[8], col[9], col[10])) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "none") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

    Output_popanc <- list(Pop_anc, Pop_map)
    names(Output_popanc) <- c("Population Ancestry Matrix", "Population Piechart Map")
    return(Output_popanc)
  }

  else {
    stop("Please supply input for plot.type. The options are 'both', 'individual', or 'population'.")

  }

}
