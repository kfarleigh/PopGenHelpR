#' A function to estimate seven measures of heterozygosity using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param statistic Character. String or vector indicating the statistic to calculate. Options are any of: all; all of the statistics; Fst, Weir and Cockerham (1984) Fst; NeisD, Nei's D statistic; JostsD, Jost's D.
#' @param missing_value Character. String indicating missing data in the input data. It is assumed to be NA, but that may not be true (is likely not) in the case of geno files.
#' @param write Boolean. Whether or not to write the output to files in the current working directory. There will be one or two files for each statistic. Files will be named based on their statistic such as Fst_perpop.csv.
#' @param prefix Character. Optional argument. String that will be appended to file output. Please provide a prefix if write is set to TRUE.
#' @param population_col Numeric. Optional argument (a number) indicating the column that contains the population assignment information.
#' @param individual_col Numeric. Optional argument (a number) indicating the column that contains the individuals (i.e., sample name) in the data.

#' @return A list containing the estimated heterozygosity statistics. The per pop values are calculated by taking the average of the per locus estimates.
#'
#' @references
#' \bold{Fst:}
#'
#' \href{https://doi.org/10.1111/1755-0998.12129}{Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013)}. StAMPP: An R package for calculation of genetic differentiation and structure of mixed‐ploidy level populations. Molecular ecology resources, 13(5), 946-952.
#'
#' \href{https://www.jstor.org/stable/2408641}{Weir, B. S., & Cockerham, C. C. (1984)}. Estimating F-statistics for the analysis of population structure. evolution, 1358-1370.
#'
#' \bold{Nei's D:}
#'
#' \href{https://www.journals.uchicago.edu/doi/abs/10.1086/282771}{Nei, M. (1972)}. Genetic distance between populations. The American Naturalist, 106(949), 283-292.
#'
#' \href{https://doi.org/10.1111/1755-0998.12129}{Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013)}. StAMPP: An R package for calculation of genetic differentiation and structure of mixed‐ploidy level populations. Molecular ecology resources, 13(5), 946-952.
#'
#' \bold{Jost's D:}
#'
#' \href{https://doi.org/10.1111/j.1365-294X.2008.03887.x}{Jost L (2008)}. GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015–4026.
#'
#' @author Keaka Farleigh
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
  } else if(tools::file_ext(data) == 'vcf') {
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
  } else if(tools::file_ext(data) == 'geno'){
    Dat <- utils::read.table(data)
    print("Geno file detected, proceeding to formatting. Note, PopGenHelpR assumes that your individuals in the geno file and
          popmap are in the same order, please check to avoid erroneous inferences.")
  } else {
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
  ########################################
  ##### Differentiation calculations #####
  ########################################

  ### Fst
  # We will use the Dat_perpop object


  Fst <- function(Dat){
    # Dat is the Dat_perpop object
    # Create a couple of matrices to store data
    ind.nomissing <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))
    q.freq <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))
    het.freq <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))

    # Calculate the alternate allele frequency and the frequency of heterozygotes per population per locus
    for(i in 1:length(Dat)){
    tmp <- Dat[[i]]
    # Get frequency of alternate alleles
    q.freq[i,] <- (((colSums(tmp[3:ncol(tmp)] == 2, na.rm = T))*2) + colSums(tmp[3:ncol(tmp)] == 1, na.rm = T))/(2*colSums(tmp[3:ncol(tmp)] != "NA"))
    # Get number of heterozygotes in each population
    het.freq[i,] <- colSums(tmp[3:ncol(tmp)] == 1)/colSums(tmp[3:ncol(tmp)] != "NA")
    # Get the number of individuals with data
    ind.nomissing[i,] <- colSums(tmp[3:ncol(tmp)] != "NA")
    }
    # Set the names of the matrices
    row.names(q.freq) <- row.names(het.freq) <- row.names(ind.nomissing) <- names(Dat)


    # Use StamPP's indexing
    index2_a <- index1_a <- NULL
    step <- 2
    step2 <- 1
    for(i in step:length(unique(Pops))){
      index1_a <- c(index1_a, i:length(unique(Pops)))
      index2_a <- c(index2_a, rep(step2, length(step:length(unique(Pops)))))
      step=step+1
      step2=step2+1
    }

    # Set r (number of pops)
    r <- 2


    nPop1 <- ind.nomissing[index1_a,]
    nPop2 <- ind.nomissing[index2_a,]

    q1 <- q.freq[index1_a,]
    q2 <- q.freq[index2_a,]

    h1 <- het.freq[index1_a,]
    h2 <- het.freq[index2_a,]

    ### WC Fst, following Weir and Cockerham 1984

    # Calculate nbar, average sample size
    n.bar <- (nPop1+nPop2)/r

    # Squared coefficient of variation of sample sizes
    nc <- (r*n.bar)-(((nPop1^2)+(nPop2^2)) / (r*n.bar))

    # Average sample frequency of allele A (q in this case)
    p.bar <- ((nPop1*q1)/(r*n.bar)) + ((nPop2*q2)/(r*n.bar))

    # Sample variance of allele A frequency over populations
    s.square <- ((nPop1*((q1-p.bar)^2))/n.bar) + ((nPop2*((q2-p.bar)^2))/n.bar)

    # Average heterozygote frequency of Allele A
    h.bar <- ((nPop1*h1)/(r*n.bar)) + ((nPop2*h2)/(r*n.bar))


    a <- (n.bar/nc) * (s.square - (1/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - ((1/4)*h.bar) ))
    b <- (n.bar/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - (((2*n.bar-1)/(4*n.bar))*h.bar))
    c <- (1/2)*h.bar

    # a, b, and c can be infinite, we will set these values to NA
    if(any(is.finite(a) == FALSE) | any(is.finite(b) == FALSE) | any(is.finite(c) == FALSE)){
      a[which(is.infinite(a))] <- NA
      b[which(is.infinite(b))] <- NA
      c[which(is.infinite(c))] <- NA
    }

    ### StAMPP's code to consolidate Fst estimates

    if(ncol(q.freq)>1){ #if there is more than 1 locus in the genotype dataset

      if(length(unique(Pops))>2){ #if there are greater than 2 populations, ie. greater than 1 pairwise comparision

        fst <- rowSums(a, na.rm=TRUE) / (rowSums(a, na.rm=TRUE) + rowSums(b, na.rm=TRUE) + rowSums(c, na.rm=TRUE)) #Fst results

      } else{ #if there is only 2 populations, ie. 1 pairwise comparision

        fst <- sum(a, na.rm=TRUE) / (sum(a, na.rm=TRUE) + sum(b, na.rm=TRUE) + sum(c, na.rm=TRUE)) #Fst results

      }

    } else{ # if there is only 1 locus in the genotype dataset

      fst <- a/(a+b+c) #Fst results

    }

    fstmat <- matrix(NA, nrow=length(unique(Pops)), ncol=length(unique(Pops)), dimnames=list(names(Dat), names(Dat)))
    fstmat[lower.tri(fstmat)]=fst
    return(fstmat)
  }

  Fst_wc <- Fst(Dat_perpop)


  ### Nei's D
  ## Calculates individual and population Nei's D

  NeiD <- function(Dat){
    nind <- nrow(Dat)
    npop <- length(Dat)

    # Calculate population Nei's D

    # Create a matrix to store the alternate allele frequency for each population
    q.freq <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))

    # First we calculate allele frequencies of the alternate allele in each population
    for(i in 1:length(Dat)){
      tmp <- Dat[[i]]
      # Get frequency of alternate alleles
      q.freq[i,] <- (((colSums(tmp[3:ncol(tmp)] == 2, na.rm = T))*2) + colSums(tmp[3:ncol(tmp)] == 1, na.rm = T))/(2*colSums(tmp[3:ncol(tmp)] != "NA"))
    }
    # Set the names of the matrices
    row.names(q.freq)  <- names(Dat)

    # Make the allele frequencies a numeric matrix
    q.freq <- as.matrix(q.freq)

    p.freq <- as.matrix(1-q.freq)

    freq_comb <- cbind(q.freq, p.freq)

    freq_comb2 <- freq_comb^2/r

    Jx <- sum(freq_comb2[1,])
    Jy <- sum(freq_comb2[2,])

    Jxy <- freq_comb[1,]*freq_comb[2,]/r
    Jxy <- sum(Jxy)

    ND <- Jxy/(sqrt(Jx*Jy))
    -log(D)

    # POPPR method, ours works above
    IDMAT2 <- freq_comb %*% t(freq_comb)
    vec2 <- sqrt(diag(IDMAT2))
    IDMAT2 <- IDMAT2/vec2[col(IDMAT2)]
    IDMAT2 <- IDMAT2/vec2[row(IDMAT2)]

    ND <- -log(IDMAT2)

    return(ND)
  }

  ### Jost's D
  JostD <- function(Dat){

    return(JD)
  }


  Output <- list(Fst_wc, ND, JD)

  names(Output) <- c("Fst", "NeisD", "JostsD")
  ### Write output
  Stat <- c("Fst", "NeisD", "JostsD")
  Stat_idx <- c(1,2,3)

  if(length(statistic) == 1 && statistic ==  "all"){
    return(Output)
  } else {
    res <- which(Stat %in% statistic)
    Output_final<- Output[which(Stat_idx %in% res)]
    return(Output_final)
  }

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