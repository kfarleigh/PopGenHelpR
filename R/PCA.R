#' A function to perform principal component analysis (PCA) on genetic data. Loci with missing data will be removed prior to PCA.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param center Boolean. Whether or not to center the data before principal component analysis.
#' @param scale Boolean. Whether or not to scale the data before principal component analysis.
#' @param missing_value Character. String indicating missing data in the input data. It is assumed to be NA, but that may not be true (is likely not) in the case of geno files.
#' @param write Boolean. Whether or not to write the output to files in the current working directory. There will be two files, one for the individual loadings and the other for the percent variance explained by each axis.
#' @param prefix Character. Optional argument. String that will be appended to file output. Please provide a prefix if write is set to TRUE.
#' @return A list containing two elements: the loadings of individuals on each principal component and the variance explained by each principal component.
#'
#' @author Keaka Farleigh
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_VCF")
#' Test <- PCA(data = HornedLizard_VCF)}
PCA <- function(data, center = TRUE, scale = FALSE, missing_value = NA, write = FALSE, prefix = NULL) {

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

  if(any(is.na(Dat))){
    n_NA <- colSums(is.na(Dat))
    Loc2rem <- which(n_NA > 0)

    print(paste(colnames(Dat)[Loc2rem], " contain NAs which cannot occur in PCA, they have been removed from subsequent analysis. Please impute these sites if you wish to include them. Alternatively, see the R packages LEA and adegenet.", sep = ''))

    Dat <- Dat[,-Loc2rem]
  }

  # Scale the data
  Dat_scaled <- scale(Dat, center = center, scale = scale)

  # Run the pca
  Dat_pca <- stats::prcomp(Dat_scaled)

  # Calculate variance explained by each principal component
  Dat_pc_var <- summary(Dat_pca)$importance[2:3,]*100
  rownames(Dat_pc_var) <- c("Percent variance explained", "Cumulative percent variance explained")

  Dat_loadings <- Dat_pca$x
  rownames(Dat_loadings) <- Inds

  Output <- list(Dat_loadings, Dat_pc_var)

  names(Output) <- c("Individual Loadings", "Variance Explained")

  if(write & !is.null(prefix)){

    utils::write.csv(t(Dat_pc_var), paste(prefix, "_pc_variance_explained.csv", sep = '_'))
    utils::write.csv(Dat_loadings, paste(prefix, "_pc_loadings.csv", sep = '_'))

    } else if(write & is.null(prefix)){

      utils::write.csv(t(Dat_pc_var), "Pc_variance_explained.csv")
      utils::write.csv(Dat_loadings, "Pc_loadings.csv")

    }

  return(Output)

}

