#' A function to estimate seven measures of heterozygosity using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param statistic Character. String or vector indicating the statistic to calculate. Options are any of: all; all of the statistics; Fst, Weir and Cockerham (1984) Fst; NeisD, Nei's D statistic; JostD, Jost's D; StHe, heterozygosity standardized by the average expected heterozygosity; StHo, heterozygosity standardized by the average observed heterozygosity; IR, internal relatedness; HL, homozygosity by locus.
#' @param missing_value Character. String indicating missing data in the input data. It is assumed to be NA, but that may not be true (is likely not) in the case of geno files.
#' @param write Boolean. Whether or not to write the output to files in the current working directory. There will be one or two files for each statistic. Files will be named based on their statistic such as Ho_perpop.csv or Ho_perloc.csv.
#' @param prefix Character. Optional argument. String that will be appended to file output. Please provide a prefix if write is set to TRUE.
#' @param population_col Numeric. Optional argument (a number) indicating the column that contains the population assignment information.
#' @param individual_col Numeric. Optional argument (a number) indicating the column that contains the individuals (i.e., sample name) in the data.

#' @return A list containing the estimated heterozygosity statistics. The per pop values are calculated by taking the average of the per locus estimates.
#'
#' @references
#' \bold{For Fst:}
#'
#' \href{https://www.jstor.org/stable/2408641}{Weir, B. S., & Cockerham, C. C. (1984)}. Estimating F-statistics for the analysis of population structure. evolution, 1358-1370.
#'
#' \bold{For Nei's D:}
#'
#' \href{https://www.journals.uchicago.edu/doi/abs/10.1086/282771}{Nei, M. (1972)}. Genetic distance between populations. The American Naturalist, 106(949), 283-292.
#'
#'
#' \bold{For Jost's D:}
#'
#' \href{https://doi.org/10.1111/j.1365-294X.2008.03887.x}{Jost L (2008)}. GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015–4026..
#'
#' @author \email{Keaka Farleigh (farleik@@miamioh.edu)}
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)}
Differentiation <- function(data, pops, statistic = 'all', missing_value = NA, write = FALSE, prefix = NULL, population_col = NULL, individual_col = NULL) {
  Statistic <- NULL

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
    Inds_popmap <- Pops[,1]
  } else if(!is.null(individual_col)){
    Inds_popmap <- Pops[,individual_col]
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

    # Transpose the numeric gt matrix
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

    # Transpose the numeric gt matrix
    Dat <- as.data.frame(t(as.matrix(gt)))

    # Preserve individual names
    Inds <- rownames(Dat)
    Dat <- sapply(Dat, as.numeric)
  }
  else if(tools::file_ext(data) == 'geno'){
    Dat <- utils::read.table(data)
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
    if(any(Inds != Inds_popmap)){
      warning("Sample names in the VCF and Population file may not be in the same order or samples are missing,
       The sample(s) in question are: ",
              print(paste(Inds[(Inds != Inds_popmap)], collapse = ' ')))  }
    if(is.null(population_col)){
      Pops <- as.factor(Pops[,2])
    }
    else if(!is.null(population_col)){
      Pops <- as.factor(Pops[,population_col])
    }
  } else if(tools::file_ext(data) == 'geno'){
      if(is.null(population_col)){
      Pops <- as.factor(Pops[,2])
    } else if(!is.null(population_col)){
        Pops <- as.factor(Pops[,population_col])
    }
  }

  # Replace missing data value with NA
  if(is.na(missing_value) == FALSE){
    Dat[Dat == missing_value] <- NA
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
  ##### Differentiation calculations #####
  #######################################

 if(write == TRUE && !is.null(prefix)){
    res_write <- which(Stat %in% statistic)
    Output2write <- Output[which(Stat_idx %in% res_write)]
    for (i in 1:length(Output2write)){
      utils::write.csv(Output2write, file = paste(names(Output2write[i]), ".csv", sep = "_"))
    }
  } else if(write == TRUE && is.null(prefix)){
      utils::write.csv(Output2write, file = paste(prefix, '_', names(Output2write[i]), ".csv", sep = "_"))
    }
}
