#' A function to calculate differentiation statistics and perform significance testing with a vcf file.
#'
#' @param VCF Character string indicating the name of the vcf file to be used in analysis.
#' @param pops Character string indicating the name of the population assignment file. This file should have four columns and be in the same order as your vcf file. The first column named Sample indicates the sample name. The second column named Population indicates the population assignment of each individual. The third column named Long indicates the longitude of the sample.  The fourth column named Lat indicates the latitude of the sample.
#' @param ploidy Numeric. The ploidy of the data.
#' @param statistic Character string. Options are both, FST, and NeisD.
#' @param boots Numeric. The number of boostraps to use to evaluate statistical significance. Only relevant for FST estimation.
#' @param prefix Character string that will be appended to file output.
#' @param write Boolean. Whether or not to write the output to a file in the current working directory.
#'
#' @return A list contianing data frames for the requested statistic.
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Dif_stats(VCF = HornedLizard_VCF, pops = HornedLizard_Pop,
#' ploidy = 2, statistic = "both", boots = 10, write = FALSE)}
Dif_stats <- function(VCF, pops, ploidy, statistic = "both", boots, write = FALSE, prefix = NULL) {
  # Read in files and convert to necessary formats
  if(missing(VCF)){
    stop("Please supply a vcf file for analysis")
  }
  if(methods::is(VCF, "vcfR")){
    Dat <- VCF
    Glight <- vcfR::vcfR2genlight(Dat)
  }
  else if(is.character(VCF)) {
    Dat <- vcfR::read.vcfR(VCF, verbose = FALSE)
    Glight <- vcfR::vcfR2genlight(Dat)
  }
  else{
    stop("Please supply a vcf file or object of class vcfR for analysis")
  }
  # Set the ploidy of the genlight
  adegenet::ploidy(Glight) <- ploidy
  # Make sure that the individuals in the vcf/genlight are in the same order as your popmap
  if(is.data.frame(pops)){
    Pops <- pops
  }
  else if(is.character(pops)){
    Pops <- utils::read.csv(pops)
  }
  else{
    stop("Please supply a csv file or data frame for population assignment")
  }
  P <- Pops
  if(any(Glight@ind.names %in% P[,1] == FALSE)){
    stop("Sample names in the VCF and Population file are not in the same order or samples are missing,
       The sample(s) in question are: ",
         paste(Glight@ind.names[!(Glight@ind.names %in% P[,1])], collapse = ' '))  }
  # Set the populations in the genlight
  Pops <- as.factor(Pops[,2])
  Glight@pop <- Pops

  message('Formatting has finished, moving onto calculations')
  if(statistic == "FST" && write == FALSE) {
    Fst <- StAMPP::stamppFst(Glight, nboots = boots)
    Fst_raw <- Fst$Fsts
    Fst_pval <- Fst$Pvalues

    # Set anything to 0 to < 1/# of bootstraps
    Fst_pval[Fst_pval == 0] <- paste('P-value < ', (round((1/boots),3)), collapse = '')
    Fst_est <- list(Fst_raw, Fst_pval)
    names(Fst_est) <- c("Fst estimates", "Fst p-values")
    return(Fst_est)}
  if(statistic == "FST" && write == TRUE) {
    Fst <- StAMPP::stamppFst(Glight, nboots = boots)
    Fst_raw <- Fst$Fsts
    Fst_pval <- Fst$Pvalues

    # Set anything to 0 to < 1/# of bootstraps
    Fst_pval[Fst_pval == 0] <- paste('P-value < ', (round((1/boots),3)), collapse = '')
    # Write out to files
    utils::write.csv(Fst_raw, file = paste(prefix, '_Fst.csv', sep = ''), row.names = TRUE)
    utils::write.csv(Fst_pval, file = paste(prefix, '_Fstpvalues.csv', sep = ''), row.names = TRUE)
    Fst_est <- list(Fst_raw, Fst_pval)
    names(Fst_est) <- c("Fst estimates", "Fst p-values")
    return(Fst_est)}
  else if(statistic == "NeisD" && write == FALSE) {
    # Neis D
    NeisD <- StAMPP::stamppNeisD(Glight)
    colnames(NeisD) <- rownames(NeisD)
    names(Nei_est) <- c("Neis D")
    return(Nei_est)}
  else if(statistic == "NeisD" && write == TRUE) {
    # Neis D
    NeisD <- StAMPP::stamppNeisD(Glight)
    colnames(NeisD) <- rownames(NeisD)

    # Write out to files
    utils::write.csv(NeisD, file = paste(prefix, '_NeisD.csv', sep = ''), row.names = TRUE)
    Nei_est <- list(NeisD)
    names(Nei_est) <- c("Neis D")
    return(Nei_est)}
  else if(statistic == "both" && write == FALSE) {
    Fst <- StAMPP::stamppFst(Glight, nboots = boots)
    Fst_raw <- Fst$Fsts
    Fst_pval <- Fst$Pvalues

    # Set anything to 0 to < 1/# of bootstraps
    Fst_pval[Fst_pval == 0] <- paste('P-value < ', (round((1/boots),3)), collapse = '')

    NeisD <- StAMPP::stamppNeisD(Glight)
    colnames(NeisD) <- rownames(NeisD)

    Stat_est <- list(Fst_raw, Fst_pval, NeisD)
    names(Stat_est) <- c("Fst estimates", "Fst p-values", "Neis D")
    return(Stat_est)
  }
  else if(statistic == "both" && write == TRUE) {
    Fst <- StAMPP::stamppFst(Glight, nboots = boots)
    Fst_raw <- Fst$Fsts
    Fst_pval <- Fst$Pvalues

    # Set anything to 0 to < 1/# of bootstraps
    Fst_pval[Fst_pval == 0] <- paste('P-value < ', (round((1/boots),3)), collapse = '')

    NeisD <- StAMPP::stamppNeisD(Glight)
    colnames(NeisD) <- rownames(NeisD)

    utils::write.csv(Fst_raw, file = paste(prefix, '_Fst.csv', sep = ''), row.names = TRUE)
    utils::write.csv(Fst_pval, file = paste(prefix, '_Fstpvalues.csv', sep = ''), row.names = TRUE)
    utils::write.csv(NeisD, file = paste(prefix, '_NeisD.csv', sep = ''), row.names = TRUE)
    Stat_est <- list(Fst_raw, Fst_pval, NeisD)
    names(Stat_est) <- c("Fst estimates", "Fst p-values", "Neis D")
    return(Stat_est)
  }

  message("Calculations have finished, the packages used to perform file formatting and calculations were
  vcfR for formatting and StaMPP to calculate differentiation statistics")

}
