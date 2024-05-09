#' A function to estimate three measures of genetic differentiation using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.
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
#' Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013). StAMPP: An R package for calculation of genetic differentiation and structure of mixed‐ploidy level populations. Molecular ecology resources, 13(5), 946-952.\doi{10.1111/1755-0998.12129}
#'
#' \href{https://www.jstor.org/stable/2408641}{Weir, B. S., & Cockerham, C. C. (1984)}. Estimating F-statistics for the analysis of population structure. evolution, 1358-1370.
#'
#' \bold{Nei's D:}
#'
#' Nei, M. (1972). Genetic distance between populations. The American Naturalist, 106(949), 283-292.\doi{10.1086/282771} 
#'
#' \doi{10.1111/1755-0998.12129} Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013). StAMPP: An R package for calculation of genetic differentiation and structure of mixed‐ploidy level populations. Molecular ecology resources, 13(5), 946-952.
#'
#' \bold{Jost's D:}
#'
#' Jost L (2008). GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015–4026.\doi{10.1111/j.1365-294X.2008.03887.x}
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
  nind <- nrow(Dat)
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

  if("Fst" %in% statistic | statistic == "all"){
    Fst_wc <- Fst(Dat_perpop)
  } else{
    Fst_wc <- NULL
  }


  ### Nei's D
  ## Calculates individual and population Nei's D
  # Assumes that Dat is a list, such as Dat_perpop

  NeisD <- function(Dat){
    nind <- nind
    npop <- length(Dat)

    # Calculate population Nei's D

    # Create a matrix to store the alternate allele frequency for each population
    q.freq <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))

    # First we calculate allele frequencies of the alternate allele in each population
    for(i in 1:length(Dat)){
      tmp <- Dat[[i]]
      # Get frequency of alternate alleles
      q.freq[i,] <- (((colSums(tmp[3:ncol(tmp)] == 2, na.rm = T))*2) + colSums(tmp[3:ncol(tmp)] == 1, na.rm = T))/(colSums(!is.na(tmp[3:ncol(tmp)]))*2)
    }
    # Set the names of the matrices
    row.names(q.freq)  <- names(Dat)

    # Get the comparisons
    Comps <- utils::combn(rownames(q.freq), m = 2)

    # Set the results matrix
    ND_res <- matrix(ncol = nrow(q.freq), nrow = nrow(q.freq))
    rownames(ND_res) <- colnames(ND_res) <- row.names(q.freq)

    # Make the allele frequencies a numeric matrix
    q.freq <- as.matrix(q.freq)

    p.freq <- as.matrix(1-q.freq)

    freq_comb <- cbind(q.freq, p.freq)


    for(i in 1:ncol(Comps)){
      Comp <- Comps[,i]
      Pops2comp <- which(rownames(q.freq) %in% Comp)

      # Remove any NaNs, calcualte r
      freq_comb_comp <- freq_comb[Pops2comp,]
      freq_comb_comp <- t(stats::na.omit(t(freq_comb_comp)))

      r <- ncol(freq_comb_comp)/2
      freq_comb_comp2 <- freq_comb_comp^2/r


      Jx <- sum(freq_comb_comp2[1,])
      Jy <- sum(freq_comb_comp2[2,])
      Jxy <- freq_comb_comp[1,]*freq_comb_comp[2,]/r
      Jxy <- sum(Jxy)
      ND <- -log(Jxy/(sqrt(Jx*Jy)))

      row.idx <- which(rownames(ND_res) == Comp[2])
      col.idx <- which(rownames(ND_res) == Comp[1])

      ND_res[row.idx, col.idx] <- ND
      diag(ND_res) <- 0
    }

    return(ND_res)
  }

  if("NeisD" %in% statistic | statistic == "all"){
    ND_pop <- NeisD(Dat_perpop)

    ## Run Nei's D per individual, we will treat each individual as a population
    Dat_perind <- Dat
    Dat_perind$P <- Dat_perind$Inds

    Dat_perind2 <- list()
    for(i in unique(Dat_perind$P)){
      Dat_perind2[[i]] <- Dat_perind[which(Dat_perind[,2] == i),]
    }

    ND_ind <- NeisD(Dat_perind2)
  } else{
    ND_pop <- ND_ind <- NULL
  }

  ### Jost's D
  # Equation 11 from Jost, 2008
  # Dat is Dat_perpop

  JostD <- function(Dat){

    # Get the number of individuals per population and the number of populations
    n.perpop <- as.data.frame(lapply(Dat, nrow))
    n.pop <- as.numeric(length(Dat))

    ### Get allele frequency at each locus

    # Create a matrix to store the alternate allele frequency for each population
    q.freq <- matrix(nrow = length(Dat), ncol = length(Dat[3:ncol(Dat[[1]])]))

    # First we calculate allele frequencies of the alternate allele in each population
    for(i in 1:length(Dat)){
      tmp <- Dat[[i]]
      # Get frequency of alternate alleles
      q.freq[i,] <- (((colSums(tmp[3:ncol(tmp)] == 2, na.rm = T))*2) + colSums(tmp[3:ncol(tmp)] == 1, na.rm = T))/(colSums(!is.na(tmp[3:ncol(tmp)]))*2)
    }
    # Set the names of the matrices
    row.names(q.freq)  <- names(Dat)

    # Make the allele frequencies a numeric matrix
    q.freq <- as.matrix(q.freq)

    p.freq <- as.matrix(1-q.freq)

    # Combine allele frequency matrices
    freq_comb <- cbind(q.freq, p.freq)

    # Squared allele frequencies
    freq_comb2 <- freq_comb^2


    ### Calculate Hs and Ht per locus for D calculations
    # freqs is the freq_comb object
    # freqs2 is the freq_comb2 object
    # pop.sizes is the n.perpop object
    D_calc <- function(freqs, freqs2, pop.sizes){
      Loc_Hets <- list('Hs' = c(), 'Ht' = c(), 'D' = c())
      nloc <- ncol(freqs)/2
      ntot <- nrow(freqs)
      Nharm <- pop.sizes[,which(rownames(freqs) %in% colnames(pop.sizes))]
      Nharm <- 1/mean(as.numeric(1/Nharm))
      for(i in 1:nloc){

        ## Within population measures
        # Remove any NAs from freqs2, recalculate the harmonic mean, recount the number of populations
        freqs2_noNA <- stats::na.omit(freqs2[,c(i,(i+nloc))])

        Nharm <- pop.sizes[,which(rownames(freqs2_noNA) %in% colnames(pop.sizes))]
        Nharm <- 1/mean(as.numeric(1/Nharm))

        n <- nrow(freqs2_noNA)

        Hs_exp <- mean(1 - rowSums(freqs2_noNA))

        Hs <- (2*Nharm/(2*Nharm-1))*Hs_exp

        Loc_Hets[["Hs"]] <- c(Loc_Hets[["Hs"]], Hs)


        # Remove NA from freqs
        freqs_noNA <- stats::na.omit(freqs[,c(i,(i+nloc))])

        ## Overall measures
        Jt <- 1 - sum(colMeans(freqs_noNA)^2)
        Ht <- Jt + Hs/(2*Nharm*n)

        Loc_Hets[["Ht"]] <- c(Loc_Hets[["Ht"]], Ht)


        # D per locus
        D <- (Ht-Hs)/(1-Hs) * (n/(n-1))

        Loc_Hets[["D"]] <- c(Loc_Hets[["D"]], D)
      }

      ## Calculate global D, do not consider and NA's
      Hs.bar <-  mean(Loc_Hets$Hs, na.rm = T)
      Ht.bar <- mean(Loc_Hets$Ht, na.rm = T)
      D.bar <- (Ht.bar-Hs.bar)/(1-Hs.bar)*(ntot/(ntot-1))

      Final <- list(D.bar, Loc_Hets)
      return(Final)
    }

    # Get the comparisons
    Comps <- utils::combn(rownames(q.freq), m = 2)

    # Set the results matrix
    JD_res <- matrix(ncol = nrow(q.freq), nrow = nrow(q.freq))
    rownames(JD_res) <- colnames(JD_res) <- row.names(q.freq)

    ## Apply the function we wrote above to calculate pairwise D

    for(i in 1:ncol(Comps)){
      # Get the comparison
      Comp <- Comps[,i]

      # Isolate population names
      Pops2comp <- which(rownames(freq_comb) %in% Comp)

      # Get allele frequencies and squared allele frequencies for each population
      comp_freq <- freq_comb[Pops2comp,]
      comp_freq2 <- freq_comb2[Pops2comp,]

      # Calculate Jost's D
      D <- D_calc(freqs = comp_freq, freqs2 = comp_freq2, pop.sizes = n.perpop)
      D <- D[[1]]

      # Get indices to store the results
      row.idx <- which(rownames(JD_res) == Comp[2])
      col.idx <- which(rownames(JD_res) == Comp[1])

      # Store the results
      JD_res[row.idx, col.idx] <- D

    }
    diag(JD_res) <- 0

    return(JD_res)
  }

  if("JostsD" %in% statistic | statistic == "all"){
    JD_pop <- JostD(Dat_perpop)
  } else{
    JD_pop <- NULL
  }

  Output <- list(Fst_wc, ND_pop, ND_ind, JD_pop)

  names(Output) <- c("Fst", "NeisD_pop", "NeisD_ind", "JostsD")
  ### Write output
  Stat <- c("Fst", "NeisD", "JostsD")
  Stat_idx <- c(1,2,2,3)

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
