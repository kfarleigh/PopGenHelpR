#' A function to perform principal component analysis on genetic data.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param missing_value Character. String indicating missing data in the input data. It is assumed to be NA, but that may not be true (is likely not) in the case of geno files.
#' @param write Boolean. Whether or not to write the output to files in the current working directory. There will be one or two files for each statistic. Files will be named based on their statistic such as Ho_perpop.csv or Ho_perloc.csv.
#' @param prefix Character. Optional argument. String that will be appended to file output. Please provide a prefix if write is set to TRUE.
#' @return A list containing the estimated heterozygosity statistics. The per pop values are calculated by taking the average of the per locus estimates.
#'
#' @author Keaka Farleigh
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)}
PCA <- function(data, missing_value = NA, write = FALSE, prefix = NULL) {
  Statistic <- NULL
  
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
  
  # Replace missing data value with NA
  if(is.na(missing_value) == FALSE){
    Dat[Dat == missing_value] <- NA
  }
  
  