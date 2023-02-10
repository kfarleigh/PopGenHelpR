#' A function to plot a heatmap from a symmetric matrix.
#'
#' @param dat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. If it is a csv, the 1st row should contain the individual/population names. The columns should also be named in this fashion.
#' @param statistic Character indicating the statistic represented in the matrix, this will be used to label the plot.
#' @param col Character vector indicating the colors to be used in plotting. The vector should contain two colors, the first will be the low value, the second will be the high value.
#'
#' @return A heatmap plot
#' @export
#'
#' @examples
#' \donttest{
#' #' data(Fst_dat)
#' Fst <- Fst_dat[[1]]
#' Fstat_plot <- Pairwise_heatmap(dat = Fst, statistic = 'FST')}
Pairwise_heatmap <- function(dat, statistic, col = NULL) {
  Var2 <- Var1 <- value <- NULL
  ### Reading in the data
  if(missing(dat)){
    stop("Please supply a data file for analysis")
  }
  else if(is.data.frame(dat) == TRUE){
    dat <- dat
  }
  else if(is.character(dat) == TRUE){
    dat <- utils::read.csv(dat, header = TRUE, row.names = 1)
  }

  ### Colors
  if(is.null(col) == TRUE) {
    low.col <- "yellow"
      high.col <- "red"
  }
  else{
    low.col <- col[1]
    high.col <- col[2]
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
      ggplot2::scale_fill_gradient(low = low.col, high = high.col, name=stat) +
      ggplot2::labs(x = "Locality", y = "Locality") +
      ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::coord_fixed()
  }
  return(Heatmap)
  Heatmap
}
