#' A function to plot a heatmap from a symmetric matrix.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. If it is a csv, the 1st row should contain the individual/population names. The columns should also be named in this fashion.
#' @param statistic Character indicating the statistic represented in the matrix, this will be used to label the plot.
#' @param col Character vector indicating the colors to be used in plotting. The vector should contain two colors, the first will be the low value, the second will be the high value.
#' @param breaks Numeric. The breaks used to generate the color ramp when plotting. The number of breaks should match the number of colors.
#' 
#' @return A heatmap plot
#' @export
#'
#' @examples
#' \donttest{
#' #' data(Fst_dat)
#' Fst <- Fst_dat[[1]]
#' Fstat_plot <- Pairwise_heatmap(dat = Fst, statistic = 'FST')}
Pairwise_heatmap <- function(dat, statistic, col = c('#abd9e9','#2c7bb6','#ffffbf', '#fdae61','#d7191c'), breaks = NULL) {
  Var2 <- Var1 <- value <- NULL
  ### Reading in the data
  if(missing(dat)){
    stop("Please supply a data file for analysis")
  }
  else if(is.data.frame(dat) == TRUE){
    dat <- dat
  }
  else if(is.character(dat)){
    dat <- utils::read.csv(dat, header = TRUE, row.names = 1)
  }
  
  ### Colors
  if(is.null(col)) {
    
    col <- c('#abd9e9','#2c7bb6','#ffffbf', '#fdae61','#d7191c')
    
  } else{
    
    col <- col
  }
  
  ### Breaks
  if(is.null(breaks)){
    Breaks <- summary(reshape2::melt(as.matrix(dat), na.rm = TRUE)$value)
    Breaks <- as.numeric(Breaks[c(1,2,4,5,6)])
  }
  else{
    Breaks <- breaks
    Breaks <- as.numeric(Breaks)
  }
  
  ### Plotting the heatmap
  if(missing(statistic)) {
    stop("Please supply a statistic for plotting")
  }
  else {
    stat <- statistic
    tri <- as.matrix(dat)
    tri_melt <- reshape2::melt(tri, na.rm = TRUE)
    
    Heatmap <- ggplot2::ggplot(data = tri_melt, ggplot2::aes(Var2, Var1, fill = value))+
      ggplot2::geom_tile(color = "white")+
      ggplot2::scale_fill_gradientn(colors = col, breaks = Breaks, name=stat) +
      ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::coord_fixed()
  }
  
  Output <- Heatmap
  return(Output)
}
