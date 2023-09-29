#' A function to map differentiation statistics.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. If it is a csv, the 1st row should contain the individual/population names. The columns should also be named in this fashion.
#' @param pops Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The columns should be named Sample, containing the sample IDs; Population indicating the population assignment of the individual; Long, indicating the longitude of the sample; Lat, indicating the latitude of the sample.
#' @param neighbors Numeric. The number of neighbors to plot connections with.
#' @param col Character vector indicating the colors you wish to use for plotting.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. Users should supply 3 values if custom breaks are desired.
#' @param Lat_buffer Numeric. A buffer to customize visualization.
#' @param Long_buffer Numeric. A buffer to customize visualization.
#'
#' @return A list containing the map and the matrix used to plot the map.
#' @export
#'
#' @examples
#' \donttest{
#' data(Fst_dat)
#' Fst <- Fst_dat[[1]]
#' Loc <- Fst_dat[[2]]
#' Test <- Dif_Stats_Map(dat = Fst, pops = Loc,
#' neighbors = 2,
#' col = c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'),Lat_buffer = 1, Long_buffer = 1)}
Dif_Stats_Map <- function(dat, pops, neighbors, col, breaks = NULL, Lat_buffer, Long_buffer){
  X1 <- X2 <- X3 <- X4 <- Dif <- Long <- Lat <- alpha <- NULL
  ################### Get the data for mapping
  # Get map data
  map <- spData::world["continent"]
  states <- spData::us_states["REGION"]
  ### Make a base map for the countries of interest
  base_map <- ggplot2::ggplot() + ggplot2::geom_sf(data = map, fill = "#f4f4f4") +
    ggplot2::geom_sf(data = states, fill = ggplot2::alpha("#f4f4f4", 0))


  # Read in files
  if(is.data.frame(dat) == TRUE){
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
    stop("Plesae supply a population assignment file as a dataframe or .csv file name")
  }

  # Make sure that the diagonal of the diferentiation matrix is 0
  diag(Dif_mat) <- 0
  # Pull coordinates
  Coords <- Pops[,3:4]
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
  Dif_Mapping <- Dif_Mapping[-c(which(Dif_Mapping$value == 0)),]
  Dif_Mapping$labels <- Pop_combos

  # Make sure that Dif_Mapping is in the same order as All_pts
  if(any(table(Dif_Mapping$labels == rownames(All_pts)) == FALSE)){
    stop("The matrices are not in the correct order, please contact the package authors if you see this error")
  }
  # Append the Dif statistic to the all pts matrix
  All_pts$Dif <- Dif_Mapping$value
  # Set breaks
  if(missing(col)){
    col <- c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
  }
  if(is.null(breaks) == TRUE){
    Breaks <- summary(All_pts$Dif)
    Breaks <- as.numeric(Breaks[1:5])
  }
  else{
    Breaks <- breaks
    Breaks <- as.numeric(Breaks)
  }
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
  else if(is.numeric(neighbors) == TRUE){
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
  else if(is.character(neighbors) == TRUE){
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
  # Add buffers for mapping
  Lat_Min <- min(All_pts_toplot[,3:4]) - Lat_buffer
  Lat_Max <- max(All_pts_toplot[,3:4]) + Lat_buffer
  Long_Min <- min(All_pts_toplot[,1:2]) - Long_buffer
  Long_Max <- max(All_pts_toplot[,1:2]) + Long_buffer

  # Map it
  Dif_map <- base_map + ggplot2::coord_sf(xlim = c(Long_Min, Long_Max),  ylim = c(Lat_Min, Lat_Max)) +
    ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4), color = 'black', size = 1.05) +
    ggplot2::geom_segment(data = All_pts_toplot, ggplot2::aes(x = X1, xend = X2, y = X3, yend = X4, color = Dif), size = 1) +
    ggplot2::scale_color_gradientn(colors = col, breaks = Breaks, name = "Differentiation") +
    ggplot2::geom_point(data = Mapping, ggplot2::aes(x = Long, y = Lat), size = 3, shape = 21, fill = 'gray', color = "black") +
    ggplot2::theme(panel.grid=ggplot2::element_blank(), legend.position = "right") + ggplot2::xlab('Longitude') + ggplot2::ylab('Latitude')

  Output_difmap <- list(All_pts_toplot[,c("Dif","Comp")], Dif_map)
  names(Output_difmap) <- c("Differentiation Matrix", "Differentiation Map")
  return(Output_difmap)
}
