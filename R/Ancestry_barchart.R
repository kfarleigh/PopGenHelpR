#' Plot an ancestry matrix for individuals and(or) populations.
#'
#' @param anc.mat Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first column should be the names of each sample/population, followed by the estimated contribution of each cluster to that individual/pop.
#' @param pops Data frame or character string that supplies the input data. If it is a character string, the file should be a csv. The first two columns should indicate the sample name (first column) and the population that sample belongs to (second column). Other columns (i.e., latitude, longitude) can be present, but will not be used. 
#' @param K Numeric.The number of genetic clusters in your data set, please contact the package authors if you need help doing this.
#' @param plot.type Character string. Options are all, individual, and population. All is default and recommended, this will plot a barchart for both the individuals and populations.
#' @param col Character vector indicating the colors you wish to use for plotting.
#' @param ind.order Character vector indicating the order to plot the individuals in the individual ancestry bar chart.
#' @param pop.order Character vector indicating the order to plot the populations in the population ancestry bar chart.
#' @param legend_pos Character. The desired position of the legend. The default is "none", which removes the legend. Other options include "left", "right", "top" or "bottom". Please see the ggplot2 documentation for all of the legend placement options. 
#'
#'
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
#' Test_all <- Ancestry_barchart(anc.mat = Qmat, pops = Loc, K = 5,
#' plot.type = 'all',col = c('#d73027', '#fc8d59', '#e0f3f8', '#91bfdb', '#4575b4'))}
Ancestry_barchart <- function(anc.mat, pops, K, plot.type = 'all', col, ind.order = NULL, pop.order = NULL, legend_pos = "right"){
  Pop <- coeff <- Sample <- value <- variable <- aes  <- alpha <- ID<- NULL
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
  
  # Set the first column name of the ancestry matrix to Ind
  colnames(Ind_anc)[1] <- "Ind"
  
  # Set the column names of the population data
  colnames(Pops)[1:2] <- c("Sample", "Population")
  
  
  if(!is.null(ind.order)){
    Ind_anc_ord <- Ind_anc[order(match(Ind_anc$Ind, ind.order)),]
    Ind_anc_ord$Ind <- factor(Ind_anc_ord$Ind, levels = ind.order)
    Ind_anc <- Ind_anc_ord
  } else{
    
    # If there is no order specified then we will sort by Q, which is the standard practice if there is no specified order
    Ind_anc$grp <- colnames(Ind_anc[2:(K+1)])[apply(Ind_anc[,2:(K+1)], 1, which.max)]
    
    Ind_anc <-  Ind_anc[order(Ind_anc$grp),]
    
    # Remove the grp column
    Ind_anc <- Ind_anc[,-ncol(Ind_anc)]
    
    ind.order <- Ind_anc[,1]
    Ind_anc_ord <- Ind_anc[order(match(ind.order, Ind_anc$Ind)),]
    Ind_anc_ord$Ind <- factor(Ind_anc_ord$Ind, levels = ind.order)
    Ind_anc <- Ind_anc_ord
    
    colnames(Ind_anc)[1] <- "Sample"
  }
  
  # Convert the individual and population information to characters for plotting

  Ind_anc[,1] <- as.character(Ind_anc[,1])
  Pops[,1] <- as.character(Pops[,1])
  Pops[,2] <- as.character(Pops[,2])
  
  ######################## Barcharts

  
  if(plot.type == 'all'){
    
    # If plotting population data, make sure that the population data lists all of the individuals
    if(all(Ind_anc[,1] %in% Pops[,1])){
      print("All information needed is present, moving onto plotting.")
      
    } else{
      
      stop("The sample names in the population data (pops argument) do not match the sample names in the individual data (anc.mat argument), 
         use a command like print(paste(as.character(Ind_anc[which(Ind_anc[,1] %in% Pops[,1] == FALSE),1]))) to determine the problame samples.")
    }
    
    
    # Individual plots
    colnames(Ind_anc) <- c("Sample", paste0(rep("cluster", K), 1:K))
    Ind_anc2 <- Ind_anc
    ord <- Ind_anc2[,1]
    Ind_anc2$ID <- Ind_anc2[,1]
    Ind_anc2$ID <- factor(Ind_anc2$ID, levels = ord)
    qmatrix_melt <- reshape2::melt(Ind_anc2[2:(K+2)], id = 'ID', value = coeff)
    
    Indplot <- qmatrix_melt %>% ggplot2::ggplot(ggplot2::aes(x= ID)) +
      ggplot2::geom_bar(ggplot2::aes(y = value, fill = variable), stat = "identity", position = "fill",width = 1) +
      ggplot2::scale_fill_manual("Population", values = col[c(1:K)], labels = paste0(rep("Cluster ", K), 1:K)) +
      ggplot2::scale_color_manual(values = col[c(1:K)], guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Ancestry Proportion", x = "") +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(), legend.position = legend_pos) +
      ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      ggplot2::scale_x_discrete() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    # Reorder the population assignment file to match the order of individuals
    Pops_ord <- Pops[order(match(Pops$Sample, Ind_anc$Sample)),]
    
    # Population plots
    Pop_anc <- Ind_anc[,-1]
    Pop_anc$Pop <- Pops_ord$Population
    Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise_all(mean, na.rm = TRUE)
    Pop_anc_coeff <- Pop_anc[,c(1:(K+1))]
    
    # Order if necessary
    if(!is.null(pop.order)){
      Pop_anc_ord <-  Pop_anc_coeff[order(match(pop.order, Pop_anc_coeff$Pop)),]
      Pop_anc_ord$Pop <- factor(Pop_anc_ord$Pop, levels = pop.order)
      Pop_anc_coeff <- Pop_anc_ord
    } else{
      
      Pop_anc_coeff$grp <- colnames(Pop_anc_coeff[2:(K+1)])[apply(Pop_anc_coeff[,2:(K+1)], 1, which.max)]
      
      Pop_anc_coeff <-  Pop_anc_coeff[order(Pop_anc_coeff$grp),]
      
      #Remove the grp column
      Pop_anc_coeff <- Pop_anc_coeff[,-ncol(Pop_anc_coeff)]
      
      pop.order <- Pop_anc_coeff$Pop
      Pop_anc_ord <-  Pop_anc_coeff[order(match(pop.order, Pop_anc_coeff$Pop)),]
      Pop_anc_ord$Pop <- factor(Pop_anc_ord$Pop, levels = pop.order)
      Pop_anc_coeff <- Pop_anc_ord
    }
    
    qmatrix_melt_pop <- reshape2::melt(Pop_anc_coeff, id = 'Pop', value = coeff)
    
    Popplot <- qmatrix_melt_pop %>% ggplot2::ggplot(ggplot2::aes(x= Pop)) +
      ggplot2::geom_bar(ggplot2::aes(y = value, fill = variable), stat = "identity", position = "fill",width = 1) +
      ggplot2::scale_fill_manual("Population", values = col[c(1:K)], labels = paste0(rep("Cluster ", K), 1:K)) +
      ggplot2::scale_color_manual(values = col[c(1:K)], guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Ancestry Proportion", x = "") +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(), legend.position = legend_pos) +
      ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      ggplot2::scale_x_discrete() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Output list
    Output_all <- list(Ind_anc2, Pop_anc_coeff, Indplot, Popplot)
    names(Output_all) <- c("Individual Plotting Data Frame", "Population Plotting Data Frame", "Individual Ancestry Plot", "Population Ancestry Plot")
    
    
    print("Want to change the the text size, font, or any other formatting? See https://kfarleigh.github.io/PopGenHelpR/articles/PopGenHelpR_plotformatting.html for examples and help.")
    
    return(Output_all) 
    
    } else if(plot.type == "individual"){
      
      # Individual plots
      colnames(Ind_anc) <- c("Sample", paste0(rep("cluster", K), 1:K))
      Ind_anc2 <- Ind_anc
      ord <- Ind_anc2[,1]
      Ind_anc2$ID <- Ind_anc2[,1]
      Ind_anc2$ID <- factor(Ind_anc2$ID, levels = ord)
      qmatrix_melt <- reshape2::melt(Ind_anc2[2:(K+2)], id = 'ID', value = coeff)
      
      Indplot <- qmatrix_melt %>% ggplot2::ggplot(ggplot2::aes(x= ID)) +
        ggplot2::geom_bar(ggplot2::aes(y = value, fill = variable), stat = "identity", position = "fill",width = 1) +
        ggplot2::scale_fill_manual("Population", values = col[c(1:K)], labels = paste0(rep("Cluster ", K), 1:K)) +
        ggplot2::scale_color_manual(values = col[c(1:K)], guide = "none") +
        ggplot2::theme_minimal() +
        ggplot2::labs(y = "Ancestry Proportion", x = "") +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(), legend.position = legend_pos) +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        ggplot2::scale_x_discrete() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
      
      Output_indanc <- list(Indplot, Ind_anc2)
      names(Output_indanc) <- c("Individual Ancestry Matrix", "Individual Plotting Data Frame")
      
      
      print("Want to change the the text size, font, or any other formatting? See https://kfarleigh.github.io/PopGenHelpR/articles/PopGenHelpR_plotformatting.html for examples and help.")
      
      return(Output_indanc)
      
    } else if(plot.type == "population"){
      
      if(all(Ind_anc[,1] %in% Pops[,1])){
        print("All information needed is present, moving onto plotting.")
        
      } else{
        
        stop("The sample names in the population data (pops argument) do not match the sample names in the individual data (anc.mat argument), 
         use a command like print(paste(as.character(Ind_anc[which(Ind_anc[,1] %in% Pops[,1] == FALSE),1]))) to determine the problame samples.")
      }
      
      # Reorder the population assignment file to match the order of individuals
      Pops_ord <- Pops[order(match(Pops$Sample, Ind_anc$Sample)),]
      
      # Population plots
      Pop_anc <- Ind_anc[,-1]
      Pop_anc$Pop <- Pops_ord$Population
      Pop_anc <- Pop_anc %>% dplyr::group_by(Pop) %>% dplyr::summarise_all(mean, na.rm = TRUE)
      Pop_anc_coeff <- Pop_anc[,c(1:(K+1))]
      
      # Order if necessary
      if(!is.null(pop.order)){
        
        Pop_anc_ord <-  Pop_anc_coeff[order(match(pop.order, Pop_anc_coeff$Pop)),]
        Pop_anc_ord$Pop <- factor(Pop_anc_ord$Pop, levels = pop.order)
        Pop_anc_coeff <- Pop_anc_ord
        
      } else{
        
        Pop_anc_coeff$grp <- colnames(Pop_anc_coeff[2:(K+1)])[apply(Pop_anc_coeff[,2:(K+1)], 1, which.max)]
        
        Pop_anc_coeff <-  Pop_anc_coeff[order(Pop_anc_coeff$grp),]
        
        #Remove the grp column
        Pop_anc_coeff <- Pop_anc_coeff[,-ncol(Pop_anc_coeff)]
        
        pop.order <- Pop_anc_coeff$Pop
        Pop_anc_ord <-  Pop_anc_coeff[order(match(pop.order, Pop_anc_coeff$Pop)),]
        Pop_anc_ord$Pop <- factor(Pop_anc_ord$Pop, levels = pop.order)
        Pop_anc_coeff <- Pop_anc_ord
      }
      
      qmatrix_melt_pop <- reshape2::melt(Pop_anc_coeff, id = 'Pop', value = coeff)
      
      Popplot <- qmatrix_melt_pop %>% ggplot2::ggplot(ggplot2::aes(x= Pop)) +
        ggplot2::geom_bar(ggplot2::aes(y = value, fill = variable), stat = "identity", position = "fill",width = 1) +
        ggplot2::scale_fill_manual("Population", values = col[c(1:K)], labels = paste0(rep("Cluster ", K), 1:K)) +
        ggplot2::scale_color_manual(values = col[c(1:K)], guide = "none") +
        ggplot2::theme_minimal() +
        ggplot2::labs(y = "Ancestry Proportion", x = "") +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),legend.position = legend_pos) +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        ggplot2::scale_x_discrete() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
      
      Output_popanc <- list(Popplot, Pop_anc_coeff)
      names(Output_popanc) <- c("Population Ancestry Matrix","Population Plotting Data Frame")
      
      
      print("Want to change the the text size, font, or any other formatting? See https://kfarleigh.github.io/PopGenHelpR/articles/PopGenHelpR_plotformatting.html for examples and help.")
      
      return(Output_popanc)
      
    } else {
      
    stop("Please supply input for plot.type. The options are 'all', 'individual', or 'population'.")
    
  }
  
}