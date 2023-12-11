#' A function to estimate heterozygosity.
#'
#' @param VCF Character string indicating the name of the vcf file to be used in analysis.
#' @param pops Character string indicating the name of the population assignment file. This file should have four columns and be in the same order as your vcf file. The first column named Sample indicates the sample name. The second column named Population indicates the population assignment of each individual. The third column named Longitude indicates the longitude of the sample.  The fourth column named Latitude indicates the latitude of the sample.
#' @param ploidy Numeric. The ploidy of the data.
#' @param prefix Character string that will be appended to file output.
#' @param write Boolean. Whether or not to write the output to a file in the current working directory.
#'
#' @return A list containing the estimated diversity statistics, model output from linear regression of these statistics against latitude, and model plots.
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Div_stats(VCF = HornedLizard_VCF, pops = HornedLizard_Pop,
#' ploidy = 2, write = FALSE)}
Heterozygosity <- function(data, pops, ploidy, statistic, write = FALSE, prefix, population_col = NULL, individual_col = NULL) {
  Latitude <- Heterozygosity <- Pop <- Standard.Deviation <- NULL

  # Read in population assignment data
  if(is.data.frame(pops)){
    Pops <- pops
  } else if(is.character(pops)){
      Pops <- utils::read.csv(pops)
  } else{
      stop("Please supply a csv file or data frame for population assignment")
  }
  # Get the list of individuals from population assignment file
  if(is.null(individual_col)) {
    Inds <- Pops[,1]
  } else if(!is.null(individual_col)){
    Inds <- Pops[,individual_col]
  }

  # Read in files and convert to necessary formats
  if(missing(data)){
    stop("Please supply a file or vcfR object for analysis")
  }
  if(methods::is(data,"vcfR")){
    Dat <- data
    print("vcfR object detected, proceeding to formatting.")
    # Convert the vcf gt slot to a geno style table for calculations
    gt <- vcfR::extract.gt(Dat)
    gt[gt == "0/0"] <- 0
    gt[gt == "0/1" | gt == "1/0"] <- 1
    gt[gt == "1/1"] <- 2
    Dat <- as.data.frame(t(as.matrix(gt)))
    # Preserve individual names
    Inds <- rownames(Dat)
    Dat <- sapply(Dat, as.numeric)
  }
  else if(tools::file_ext(data) == 'vcf') {
    Dat <- vcfR::read.vcfR(data, verbose = FALSE)
    print("VCF file detected, proceeding to formatting.")
    # Convert the vcf gt slot to a geno style table for calculations
    gt <- vcfR::extract.gt(Dat)
    gt[gt == "0/0"] <- 0
    gt[gt == "0/1" | gt == "1/0"] <- 1
    gt[gt == "1/1"] <- 2
    Dat <- as.data.frame(t(as.matrix(gt)))
    # Preserve individual names
    Inds <- rownames(Dat)
    Dat <- sapply(Dat, as.numeric)
  }
  else if(tools::file_ext(data) == 'geno'){
    Dat <- utils::read.table(Data)
    print("Geno file detected, proceeding to formatting. Note, PopGenHelpR assumes that your individuals in the geno file and
          popmap are in the same order, please check to avoid erroneous inferences.")
  }
  else {
      stop("Please supply a geno file, vcf file, or vcfR object for analysis")
  }

  P <- Pops
  ### Check to see if the individuals are in the same order in the vcf data and the population assignment file
  # Make sure that the individuals in the vcf/genlight are in the same order as your popmap
  if(methods::is(data,"vcfR") || tools::file_ext(data) == 'vcf') {
    if(any(Dat[,1] != P[,1])){
      warning("Sample names in the VCF and Population file may not be in the same order or samples are missing,
       The sample(s) in question are: ",
              print(paste(Dat[,1][(Dat[,1] != P[,1])], collapse = ' ')))  }
    if(is.null(population_col)){
      Pops <- as.factor(Pops[,2])
    }
    else if(!is.null(population_col)){
      Pops <- as.factor(Pops[,population_col])
    }
  } else if(tools::file_ext(data) == 'geno'){
      if(is.null(population_col)){
      Pops <- as.factor(Pops[,2])
    }   else if(!is.null(population_col)){
        Pops <- as.factor(Pops[,population_col])
    }
  }

  P <- Pops
  Dat <- cbind.data.frame(Inds, P, Dat)

  # Break into list with populations for each element
  Dat_perpop <- list()
  for(i in unique(P)){
    Dat_perpop[[i]] <- Dat[which(Dat[,2] == i),]
  }

  message('Formatting has finished, moving onto calculations')
  #######################################
  ##### Heterozygosity calculations #####
  #######################################
  # Find observed heterozygosity (1-homozygosity)
 ObsHet <- function(Dat){
    ObsHet_perloc <- 1-(colSums(Dat[3:ncol(Dat)] != 1)/nrow(Dat))
    return(ObsHet_perloc)
  }

  ObsHet_res_perloc  <- lapply(Dat_perpop, ObsHet)
  Obs_Het_res_avg <- lapply(ObsHet_res_perloc, mean)

 ExpHe <- function(Dat){
   p <- ((colSums(Dat[3:ncol(Dat)]== 0)*2) + colSums(Dat[3:ncol(Dat)]== 1))/(nrow(Dat)*2)
   q <- ((colSums(Dat[3:ncol(Dat)]== 2)*2) + colSums(Dat[3:ncol(Dat)]== 1))/(nrow(Dat)*2)
   He_perloc <- 1-(p^2)-(q^2)
   return(He_perloc)
  }

 ExpHet_res_perloc  <- lapply(Dat_perpop, ExpHe)
 ExpHet_res_avg  <- lapply(ExpHet_res_perloc, mean)

 PropHt <-
 StHe <-
 StHo <-
 IR <-
 HL <-








  ##########################
  ##### Heterozygosity #####
  ##########################
  # Calculate heterozygosity and standard deviation for each population
  Het <- hierfstat::basic.stats(Hstat)

  # Get per populations heterozygosity
  H_all <- data.frame(colMeans(Het$Ho, na.rm = TRUE))
  H_all$Pop <- rownames(H_all)
  colnames(H_all) <- c('Heterozygosity', 'Pop')

  # Standard deviation
  H_all$StandardDeviation <- stats::sd(H_all$Heterozygosity)

  # Attach coordinates
  # First we check to make sure that the populations are in the right order
  if(any(H_all$Pop != unique(P[,2]))){
    stop("Populations are not in the correct order, if you see this please email the package authors")
  }

  # Pull populations coordinates
  Pop_coords <- P[!duplicated(P$Population),]

  # Append to heterozygosity estimates
  H_all[,4:5] <- Pop_coords[,3:4]
  colnames(H_all) <- c('Heterozygosity', 'Pop', 'Standard.Deviation','Longitude', 'Latitude')

  # Is there a statistical relationship between heterozygosity and latitude
  Het_model <- stats::lm(H_all$Heterozygosity ~ H_all$Latitude)

  message('Heterozygosity calculated, moving to private alleles')


  ##########################
  ##### Visualizations #####
  ##########################
  ##### Heterozygosity #####
  # We will make a plot of heterozygostiy estimates and an interpolated map of heterozygosity values for each population

  # Plot of heterozygosity values (y-axis) and latitude (x-axis)
  Het_plot <- ggplot2::ggplot(data = H_all, ggplot2::aes(x = Latitude, y = Heterozygosity)) + ggplot2::geom_point(ggplot2::aes(color = Pop),size = 3) +
    ggplot2::geom_errorbar(data = H_all, ggplot2::aes(x = Latitude, ymin = Heterozygosity-Standard.Deviation, ymax = Heterozygosity+Standard.Deviation,
                                                      color = Pop)) +
    ggplot2::stat_smooth(method = 'lm', formula = y ~ x, size =1, colour = 'black') +
    ggplot2::labs(x = 'Latitude', y = 'Heterozygosity') +
    ggplot2::ggtitle(expression(atop("Obersved Heterozygosity of Localities"))) +
    ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                              legend.title = ggplot2::element_blank())


  Output <- list(H_all, Het_model, Het_plot)
  names(Output) <- c('Heterozygosity_calculations', "Heterozygosity_vs_Latitude_Model",
                     "Heterozygosity_vs_Latitude_plot")
  if(write == TRUE){
    # Write out heterozygosity results
    utils::write.csv(H_all, file = paste(as.character(prefix), '_Heterozygosity.csv', sep = ''), row.names = FALSE)
    summary(Het_model)
    # Write out linear regression results heterozygosity ~ latitude
    sink(paste(as.character(prefix), '_Heterozygosity_lm.txt', sep = ''))
    summary(Het_model)
    sink()
  }

  message("Calculations have finished, the packages used to perform file formatting and calculations were
  vcfR, adegenet, and dartR for formatting, hierfstat to calculate heterozygosity, and poppr to calculate private alleles")

  return(Output)
}
