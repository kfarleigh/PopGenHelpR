#' A function to estimate heterozygosity.
#'
#' @param data Character. String indicating the name of the vcf file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param statistic Character. String indicating the statistic to calculate. Options are ...
#' @param write Boolean. Whether or not to write the output to a file in the current working directory.
#' @param prefix Character. String that will be appended to file output.
#' @param population_col Numeric. Optional argument (a number) indicating the column that contains the population assignment information.
#' @param individual_col Numeric. Optional argument (a number) indicating the column that contains the individuals (i.e., sample name) in the data.

#' @return A list containing the estimated diversity statistics, model output from linear regression of these statistics against latitude, and model plots.
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Div_stats(VCF = HornedLizard_VCF, pops = HornedLizard_Pop,
#' ploidy = 2, write = FALSE)}
Heterozygosity <- function(data, pops, statistic, write = FALSE, prefix, population_col = NULL, individual_col = NULL) {
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

 PropHt <- NULL
 StHe <- NULL
 StHo <- NULL
 IR <- NULL
 HL <- NULL

 Output <- list(Obs_Het_res_avg, ObsHet_res_perloc)

return(Output)
}
