#' A function to estimate heterozygosity and the number of private alleles from a vcf file.
#'
#' @param VCF Character string indicating the name of the vcf file to be used in analysis.
#' @param pops Character string indicating the name of the population assignment file. This file should have four columns and be in the same order as your vcf file. The first column named Sample indicates the sample name. The second column named Population indicates the population assignment of each individual. The third column named Long indicates the longitude of the sample.  The fourth column named Lat indicates the latitude of the sample.
#' @param ploidy Numeric. The ploidy of the data.
#' @param prefix Character string that will be appended to file output.
#'
#' @return A list containing the estimated diversity statistics, model output from linear regression of these statistics against latitude, and model plots.
#' @export
#'
#' @examples
#' \dontrun{
#' Test <- Div_stats('my_vcf.vcf', 'my_popmap.csv',
#' ploidy = 2, prefix = 'Test')}
Div_stats <- function(VCF, pops, ploidy, prefix) {
  Latitude <- Heterozygosity <- Pop <- Standard.Deviation <- Private.Alleles <- NULL
  # Read in files and convert to necessary formats
  if(missing(VCF)){
    stop("Please supply a vcf file for analysis")
  }
  Dat <- vcfR::read.vcfR(VCF, verbose = FALSE)
  Glight <- vcfR::vcfR2genlight(Dat)
  # Set the ploidy of the genlight
  adegenet::ploidy(Glight) <- ploidy
  # Make sure that the individuals in the vcf/genlight are in the same order as your popmap

  Pops <- utils::read.csv(pops)
  P <- Pops
  if(any(Glight@ind.names %in% P[,1] == FALSE)){
    stop("Sample names in the VCF and Population file are not in the same order or samples are missing,
       The sample(s) in question are: ",
         paste(Glight@ind.names[!(Glight@ind.names %in% P[,1])], collapse = ' '))  }
  # Set the populations in the genlight
  Pops <- as.factor(Pops[,2])
  Glight@pop <- Pops

  # Make an genind object and convert to format for heterozygosity calculations
  Genind <- dartR::gl2gi(Glight, v = 0)
  Hstat <- hierfstat::genind2hierfstat(Genind)

  message('Formatting has finished, moving onto calculations')

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

  # Write out to a file
  utils::write.csv(H_all, file = paste(as.character(prefix), '_Heterozygosity.csv', sep = ''), row.names = FALSE)

  # Is there a statistical relationship between heterozygosity and latitude
  Het_model <- stats::lm(H_all$Heterozygosity ~ H_all$Latitude)

  summary(Het_model)
  sink(paste(as.character(prefix), '_Heterozygosity_lm.txt', sep = ''))
  summary(Het_model)
  sink()

  message('Heterozygosity calculated, moving to private alleles')

  ###########################
  ##### Private Alleles #####
  ###########################
  # Get private alleles per population
  Private_alleles <- data.frame(rowSums(poppr::private_alleles(Genind)))

  Private_alleles$Pop <- rownames(Private_alleles)
  colnames(Private_alleles) <- c('PrivateAlleles', 'Pop')
  # Attach coordinates
  # First we check to make sure that the populations are in the right order
  if(any(Private_alleles$Pop != unique(P[,2]))){
    stop("Populations are not in the correct order, if you see this please email the package authors")
  }
  Private_alleles[,3:4] <- Pop_coords[,3:4]

  colnames(Private_alleles) <- c('Private.Alleles', 'Pop','Longitude', 'Latitude')

  # Write out to a file
  utils::write.csv(Private_alleles, file = paste(as.character(prefix), '_PrivateAlleles.csv', sep = ''), row.names = FALSE)

  # Is there a statistical relationship between the number of private alleles and latitude
  PA_model <- stats::lm(Private_alleles$Private.Alleles ~ Private_alleles$Latitude)

  summary(PA_model)
  sink(paste(as.character(prefix), '_Private.Alleles_lm.txt', sep = ''))
  summary(PA_model)
  sink()

  message('Private Alleles have been calculated, moving onto plotting')

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

  ##### Private Alleles #####
  # We will make a map (interpolated and not interpolated) of the number of private alleles in each population

  # Plot of private allele values (y-axis) and latitude (x-axis)
  Paf_plot <- ggplot2::ggplot(data = Private_alleles, ggplot2::aes(x = Latitude, y = Private.Alleles)) + ggplot2::geom_point(ggplot2::aes(color = Pop), size = 3) +
    ggplot2::stat_smooth(method = 'lm', formula = y ~ x, size =1, colour = 'black') +
    ggplot2::labs(x = 'Latitude', y = 'Private Alleles') +
    ggplot2::ggtitle(expression(atop("Private Alleles of Localities"))) +
    ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                              legend.title = ggplot2::element_blank())

  Output <- list(H_all, Private_alleles, Het_model, PA_model, Het_plot, Paf_plot)
  names(Output) <- c('Heterozygosity_calculations', 'Private_Allele_Calculations', "Heterozygosity_vs_Latitude_Model",
                     "Private_Allleles_vs_Latitude_Model", "Heterozygosity_vs_Latitude_plot",
                     "PA_vs_Latitude_plot")
  return(Output)

  message("Calculations have finished, the packages used to perform file formatting and calculations were
  vcfR, adegenet, and dartR for formatting, hierfstat to calculate heterozygosity, and poppr to calculate private alleles")


}
