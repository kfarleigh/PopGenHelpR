#' A function to estimate seven measures of heterozygosity using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param statistic Character. String or vector indicating the statistic to calculate. Options are any of: all; all of the statistics; Ho, observed heterozygosity; He, expected heterozygosity; PHt, proportion of heterozygous loci; Hs_exp, heterozygosity standardized by the average expected heterozygosity; Hs_obs, heterozygosity standardized by the average observed heterozygosity; IR, internal relatedness; HL, homozygosity by locus.
#' @param missing_value Character. String indicating missing data in the input data. It is assumed to be NA, but that may not be true (is likely not) in the case of geno files.
#' @param write Boolean. Whether or not to write the output to files in the current working directory. There will be one or two files for each statistic. Files will be named based on their statistic such as Ho_perpop.csv or Ho_perloc.csv.
#' @param prefix Character. Optional argument. String that will be appended to file output. Please provide a prefix if write is set to TRUE.
#' @param population_col Numeric. Optional argument (a number) indicating the column that contains the population assignment information.
#' @param individual_col Numeric. Optional argument (a number) indicating the column that contains the individuals (i.e., sample name) in the data.
#' @return A list containing the estimated heterozygosity statistics. The per pop values are calculated by taking the average of the per locus estimates.
#'
#' @references
#' \bold{Expected (He) and observed heterozygosity (Ho):}
#'
#' Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#' \bold{Homozygosity by locus (HL) and internal relatedness (IR):}
#'
#' \href{https://onlinelibrary.wiley.com/doi/full/10.1111/j.1755-0998.2010.02830.x?casa_token=QiNcMSJyunkAAAAA%3Agv-CK7GrUn1bHSgz4qZSOcB2nyHDeR8B1Wtm9bM7q7vZCAcJhNkhTWnpM0EfkSCb2EvkRrr2ArMzC7v7}{Alho, J. S., Välimäki, K., & Merilä, J. (2010)}. Rhh: an R extension for estimating multilocus heterozygosity and heterozygosity–heterozygosity correlation. Molecular ecology resources, 10(4), 720-722.
#'
#' Amos, W., Worthington Wilmer, J., Fullard, K., Burg, T. M., Croxall, J. P., Bloch, D., & Coulson, T. (2001). The influence of parental relatedness on reproductive success. Proceedings of the Royal Society of London. Series B: Biological Sciences, 268(1480), 2021-2027.\doi{10.1098/rspb.2001.1751}
#'
#' \href{https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2006.03111.x}{Aparicio, J. M., Ortego, J., & Cordero, P. J. (2006)}. What should we weigh to estimate heterozygosity, alleles or loci?. Molecular Ecology, 15(14), 4659-4665.
#'
#' \bold{Heterozygosity standardized by expected (Hs_exp) and observed heterozygosity (Hs_obs):}
#'
#' Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J. M. (1999). Parasite‐mediated selection against Inbred Soay sheep in a free‐living island population. Evolution, 53(4), 1259-1267.\doi{10.1111/j.1558-5646.1999.tb04538.x}
#'
#' @author Keaka Farleigh
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)}
Heterozygosity <- function(data, pops, statistic = 'all', missing_value = NA, write = FALSE, prefix = NULL, population_col = NULL, individual_col = NULL) {
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
    # Detect the seperator, "/" is generally used for unphased genotypes, "|" for phased, and SLiM.
    if(grepl("/", Dat@gt[1,2])){

      # Convert the vcf gt slot to a geno style table for calculations
      gt <- vcfR::extract.gt(Dat)
      gt[gt == "0/0"] <- 0
      gt[gt == "0/1" | gt == "1/0"] <- 1
      gt[gt == "1/1"] <- 2
    } else{
      gt <- vcfR::extract.gt(Dat)
      gt[gt == "0|0"] <- 0
      gt[gt == "0|1" | gt == "1|0"] <- 1
      gt[gt == "1|1"] <- 2
    }

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
    # Detect the seperator, "/" is generally used for unphased genotypes, "|" for phased, and SLiM.
    if(grepl("/", Dat@gt[1,2])){

      # Convert the vcf gt slot to a geno style table for calculations
      gt <- vcfR::extract.gt(Dat)
      gt[gt == "0/0"] <- 0
      gt[gt == "0/1" | gt == "1/0"] <- 1
      gt[gt == "1/1"] <- 2
    } else{
      gt <- vcfR::extract.gt(Dat)
      gt[gt == "0|0"] <- 0
      gt[gt == "0|1" | gt == "1|0"] <- 1
      gt[gt == "1|1"] <- 2
    }
    # Transpose the numeric gt matrix
    Dat <- as.data.frame(t(as.matrix(gt)))

    # Preserve individual names
    Inds <- rownames(Dat)
    Dat <- sapply(Dat, as.numeric)
  }
  else if(tools::file_ext(data) == 'geno'){
    Dat <- utils::read.table(data)
    Dat <- sapply(Dat, as.numeric)
    Inds <- Inds_popmap
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
  ##### Heterozygosity calculations #####
  #######################################
  # Estimate observed heterozygosity (1-homozygosity), accounts for NA
  ObsHet <- function(Dat){
    ObsHet_perloc <- 1-(colSums(Dat[3:ncol(Dat)] != 1, na.rm = T)/(nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)]))))
    return(ObsHet_perloc)
  }

  if("Ho" %in% statistic | statistic == "all"){
    # Calculate observed heterozygosity
    ObsHet_res  <- lapply(Dat_perpop, ObsHet)

    # Calculate the average per population
    Obs_Het_res_avg <- lapply(ObsHet_res, stats::na.omit)
    Obs_Het_res_avg <- lapply(Obs_Het_res_avg, mean)

    Obs_Het_res_avg <- as.data.frame(do.call("rbind", Obs_Het_res_avg))
    Obs_Het_res_avg$Pop <- rownames(Obs_Het_res_avg)
    colnames(Obs_Het_res_avg)[1] <- "Observed.Heterozygosity"
    Obs_Het_res_avg[,1] <- as.numeric(Obs_Het_res_avg[,1])

    # Make the results per locus easier to view
    ObsHet_res_perloc <- as.data.frame(do.call("rbind", ObsHet_res))
    ObsHet_res_perloc[ObsHet_res_perloc == "NaN"] <- NA
    rownames(ObsHet_res_perloc) <- names(ObsHet_res)

  } else{
    ObsHet_res_perloc <- Obs_Het_res_avg <- NULL
  }

  # Estimate expected heterozygosity
  ExpHe <- function(Dat){
    p <- as.data.frame(((colSums(Dat[3:ncol(Dat)]== 0, na.rm = T)*2) + colSums(Dat[3:ncol(Dat)]== 1, na.rm = T))/((nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)])))*2))
    p <- t(p)
    rownames(p) <- NULL

    colnames(p) <- paste(colnames(p), ".a", sep = "")

    q <- as.data.frame(((colSums(Dat[3:ncol(Dat)]== 2, na.rm = T)*2) + colSums(Dat[3:ncol(Dat)]== 1, na.rm = T))/((nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)])))*2))
    q <- t(q)
    rownames(q) <- NULL

    colnames(q) <- paste(colnames(q), ".b", sep = "")

    # Remove loci with only missing data
    loc2rem <- which(colSums(is.na(Dat[3:ncol(Dat)])) == nrow(Dat))

    # Set those values to NA
    p[loc2rem] <- NA
    q[loc2rem] <- NA

    # Combine the two
    pq <- cbind(p, q)

    # Square and add the allele frequencies, removing NA
    pq_comb <- apply(pq^2, 1, sum, na.rm = TRUE)

    # Calculate He per loc
    He_perloc <- 1-(p^2)-(q^2)

    colnames(He_perloc) <- gsub(".a", "", colnames(He_perloc))

    # How many loci were genotyped in the population?
    loc_gtyped <- length(p) - length(loc2rem)

    # Get the average He, divided by the number of loci genotyped in that population
    He_avg <- 1-(pq_comb/loc_gtyped)

    res_list <- list()

    res_list[[1]] <- He_perloc
    res_list[[2]] <- He_avg

    return(res_list)
  }

  if("He" %in% statistic | statistic == "all"){
    ExpHet_res  <- lapply(Dat_perpop, ExpHe)

    ExpHet_res_avg <- as.data.frame(do.call(rbind, lapply(ExpHet_res, function(x) x[[2]])))
    colnames(ExpHet_res_avg) <- "Expected.Heterozygosity"

    ExpHet_res_avg$Pop <- rownames(ExpHet_res_avg)
    ExpHet_res_avg[,1] <- as.numeric(ExpHet_res_avg[,1])

    # Make it into a dataframe so that it is easier for users to view
    ExpHet_res_perloc <- lapply(ExpHet_res, function(x) x[[1]])

    ExpHet_res_perloc <- as.data.frame(do.call("rbind", ExpHet_res_perloc))

    rownames(ExpHet_res_perloc) <- names(ExpHet_res)

  } else{
    ExpHet_res_perloc <- ExpHet_res_avg <- NULL
  }

  # Estimate the proportion of heterozygous loci per individual
  PropHt <- function(Dat){
    PHt <- t(as.data.frame(rowSums(Dat[3:ncol(Dat)] == 1, na.rm = T)/((ncol(Dat)-2) - rowSums(is.na(Dat[3:ncol(Dat)])))))
    rownames(PHt) <- 'PHt'
    colnames(PHt) <- Dat[,1]

    return(PHt)
  }

  if("PHt" %in% statistic | statistic == "all"){
    PropHt_res_perind <- lapply(Dat_perpop, PropHt)
    PropHt_res_perind <- mapply(rbind, PropHt_res_perind, "Pop"=names(PropHt_res_perind), SIMPLIFY=F)

    # Make it into a dataframe so that it is easier for users to view
    PropHt_res_perind <- t(as.data.frame(do.call("cbind", PropHt_res_perind)))
    PropHt_res_perind <- as.data.frame(PropHt_res_perind)
    PropHt_res_perind[,1] <- as.numeric(PropHt_res_perind[,1])
  } else{
    PropHt_res_perind <- NULL
  }

  # Estimate standardized heterozygosity based on the expected heterozygosity
  StHe <- function(Dat){

    ExpHet_res <- ExpHe(Dat)
    ExpHet_res <- t(as.data.frame(ExpHet_res[[1]]))

    ExpHet_res_avg <- mean(ExpHet_res[,1], na.rm = TRUE)

    St_He <- PropHt(Dat)/ExpHet_res_avg
    rownames(St_He) <- 'Hs_exp'
    return(St_He)
  }
  if("Hs_exp" %in% statistic | statistic == "all"){
    StHe_res_perind <- lapply(Dat_perpop, StHe)
    StHe_res_perind <- mapply(rbind, StHe_res_perind, "Pop"=names(StHe_res_perind), SIMPLIFY=F)

    # Make it into a dataframe so that it is easier for users to view
    StHe_res_perind <- t(as.data.frame(do.call("cbind", StHe_res_perind)))
    StHe_res_perind <- as.data.frame(StHe_res_perind)
    StHe_res_perind[,1] <- as.numeric(StHe_res_perind[,1])
  } else{
    StHe_res_perind <- NULL
  }
  # Estimate standardized heterozygosity based on the observed heterozyosity
  StHo <- function(Dat){

    ObsHet_res_perloc  <- ObsHet(Dat)
    ObsHet_res_avg  <- stats::na.omit(ObsHet_res_perloc)
    ObsHet_res_avg  <- mean(ObsHet_res_avg, na.rm = TRUE)
    St_Ho <- PropHt(Dat)/ObsHet_res_avg
    rownames(St_Ho) <- 'Hs_obs'
    return(St_Ho)
  }
  if("Hs_obs" %in% statistic | statistic == "all"){
    StHo_res_perind <- lapply(Dat_perpop, StHo)
    StHo_res_perind <- mapply(rbind, StHo_res_perind, "Pop"=names(StHo_res_perind), SIMPLIFY=F)

    # Make it into a dataframe so that it is easier for users to view
    StHo_res_perind <- t(as.data.frame(do.call("cbind", StHo_res_perind)))
    StHo_res_perind <- as.data.frame(StHo_res_perind)
    StHo_res_perind[,1] <- as.numeric(StHo_res_perind[,1])
  } else{
    StHo_res_perind <- NULL
  }

  IR <- function(Dat){
    # Extract genetic data and convert to count the frequency of allele s
    tmp <- Dat
    tmp[tmp == 0] <- "0/0"
    tmp[tmp == 1] <- "0/2"
    tmp[tmp == 2] <- "2/2"
    tmp[is.na(tmp)] <- NA/NA

    if(any(colSums(is.na(tmp)) == nrow(tmp))){
      print(paste(names(which(colSums(is.na(tmp)) == nrow(tmp))), " has no data (everything is NA), please check your data.", sep = ''))
    }

    rownames(tmp) <- Inds

    gtypes <- tmp[,3:ncol(tmp)]

    Loc_IR_formatted <- data.frame(matrix(NA, nrow = nrow(tmp), ncol = (ncol(gtypes)*2)))
    rownames(Loc_IR_formatted) <- tmp[,1]

    for(i in (1:ncol(gtypes))){
      # Get the locus name
      Loc_nam <- colnames(gtypes)[i]

      # Split the genotype
      Loc_split <- strsplit(gtypes[,i], "/", fixed = T)
      # Set the name of individuals
      names(Loc_split) <- Dat[,1]

      # Bind into a data frame
      Loc_split_df <- do.call(rbind, Loc_split)

      # Set an index to position the genotypes correctly
      idx <- i*2-1

      Loc_IR_formatted[,idx:(idx+1)] <- Loc_split_df
      colnames(Loc_IR_formatted)[idx:(idx+1)] <- c(paste(Loc_nam, "a", sep = ""), paste(Loc_nam, "b", sep = ""))
    }
    # Convert to a matrix for calculations
    Loc_IR_mat <- as.matrix(Loc_IR_formatted)

    # Bind together with individual and population information
    Loc_IR_mat_comb <- cbind(tmp[,1:2], Loc_IR_mat)
    Loc_IR_mat_comb <- as.matrix(Loc_IR_mat_comb)

    # Get the number of individuals
    Individuals <- nrow(Loc_IR_mat)

    # Get the number of loci
    Nloc <- ncol(Loc_IR_mat)/2

    Nloc <- Nloc + 1

    # Set up a results table
    res_tab <- data.frame(IR = rep(NA, Individuals))
    rownames(res_tab) <- Inds

    # Get the counts of alleles for each locus
    Counts <- list()

    # Count the occurrences of each allele (0 and 2) at each locus
    for(i in 1:Nloc) {
      # Set the same index as above
      idx <- 2*i-1
      idx2 <- 2*i
      Counts[[i]] <- table(Loc_IR_mat_comb[,idx:idx2])
    }

    ### Calculate IR for each individual
    # Formula from the archived Rhh package https://cran.r-project.org/web/packages/Rhh/index.html
    for(i in 1:Individuals){

      H <- 0
      N <- 0
      f <- 0

      for(j in 1:Nloc){
        # Set our index again
        idx1 <- 2*j-1
        idx2 <- 2*j

        if((!is.na(Loc_IR_mat_comb[i,idx1])) && (!is.na(Loc_IR_mat_comb[i,idx2]))){
          N <- N + 1

          if(Loc_IR_mat_comb[i,idx1] == Loc_IR_mat_comb[i,idx2]){
            H <- H + 1
            # Which allele is the individual homozygous for
            Hom_Allele <- as.character(Loc_IR_mat_comb[i,idx1])
            f <- f + (2 * Counts[[j]][[Hom_Allele]] - 2)/(sum(Counts[[j]]) - 2)

          } else{
            # If they are heterozygous that means that they are contributing two alleles, using just a sum of 1 leads to NaN
            Het_Allele1 <- as.character(Loc_IR_mat_comb[i,idx1])
            f <- f + (Counts[[j]][[Het_Allele1]] - 1)/(sum(Counts[[j]]) - 2)
            Het_Allele2 <- as.character(Loc_IR_mat_comb[i,idx2])
            f <- f + (Counts[[j]][[Het_Allele2]] - 1)/(sum(Counts[[j]]) - 2)
          }
        }
      }
      # Calculate internal relatedness
      res_tab[i,1] <- (2 * H - f) / (2 * N - f)
    }
    return(res_tab)
  }

  if("IR" %in% statistic | statistic == "all"){
    IR_perind <- IR(Dat)
  } else{
    IR_perind <- NULL
  }


  # Estimate homozygosity by locus
  HL <- function(Dat){
    # Extract genetic data and convert to count the frequency of allele s
    tmp <- Dat
    tmp[tmp == 0] <- "0/0"
    tmp[tmp == 1] <- "0/2"
    tmp[tmp == 2] <- "2/2"
    tmp[is.na(tmp)] <- NA/NA

    gtypes <- tmp[,3:ncol(tmp)]


    if(any(colSums(is.na(tmp)) == nrow(tmp))){
      print(paste(names(which(colSums(is.na(tmp)) == nrow(tmp))), " has no data (everything is NA), please check your data.", sep = ''))
    }

    Loc_HL_formatted <- data.frame(matrix(NA, nrow = nrow(tmp), ncol = (ncol(gtypes)*2)))
    rownames(Loc_HL_formatted) <- tmp[,1]


    for(i in (1:ncol(gtypes))){
      # Get the locus name
      Loc_nam <- colnames(gtypes)[i]

      # Split the genotype
      Loc_split <- strsplit(gtypes[,i], "/", fixed = T)
      # Set the name of individuals
      names(Loc_split) <- Dat[,1]

      # Bind into a data frame
      Loc_split_df <- do.call(rbind, Loc_split)

      # Set an index to position the genotypes correctly
      idx <- i*2-1

      Loc_HL_formatted[,idx:(idx+1)] <- Loc_split_df
      colnames(Loc_HL_formatted)[idx:(idx+1)] <- c(paste(Loc_nam, "a", sep = ""), paste(Loc_nam, "b", sep = ""))
    }


    Loc_HL_mat <- as.matrix(Loc_HL_formatted)

    # Bind together with individual and population information
    Loc_HL_mat_comb <- cbind(tmp[,1:2], Loc_HL_mat)
    Loc_HL_mat_comb <- as.matrix(Loc_HL_mat_comb)

    # Get the number of individuals
    Individuals <- nrow(Loc_HL_mat_comb)

    # Get the number of loci
    Nloc <- ncol(Loc_HL_mat)/2
    Nloc <- Nloc + 1

    # Set up a results table
    res_tab <- data.frame(HL = matrix(NA, nrow = Individuals, ncol = 1), row.names = Inds)

    E <- array(Nloc)
    Counts <- array(Nloc)

    # Count the occurrences of each allele (0 and 2) at each locus
    for(i in 1:Nloc) {
      E[i] <- 1
      # Set the same index as above
      idx <- i*2-1
      Counts[i] <- list(table(Loc_HL_mat_comb[,idx:(idx+1)]))
      E[i] <- 1 - sum((Counts[[i]] / sum(Counts[[i]]))^2)
    }
    # Formula from the archived Rhh package https://cran.r-project.org/web/packages/Rhh/index.html
    for (i in 1:Individuals) {

      sum.Eh <- 0
      sum.Ej <- 0

      for (j in 1:Nloc) {

        idx1 <- 2*j-1
        idx2 <- 2*j

        if ((!is.na(Loc_HL_mat_comb[i, idx1])) && (!is.na(Loc_HL_mat_comb[i, idx2]))) {
          if (Loc_HL_mat_comb[i, idx1] == Loc_HL_mat_comb[i, idx2]) {
            sum.Eh <- sum.Eh + E[j]
          }
          else {
            sum.Ej <- sum.Ej + E[j]
          }
        }
      }

      res_tab[i,1] <- sum.Eh / (sum.Eh + sum.Ej)

    }

    return(res_tab)

  }

  if("HL" %in% statistic | statistic == "all"){
    HL_perind <- HL(Dat)
  } else{
    HL_perind <- NULL
  }


  Output <- list(Obs_Het_res_avg, ObsHet_res_perloc, ExpHet_res_avg, ExpHet_res_perloc,
                 PropHt_res_perind, StHe_res_perind, StHo_res_perind, IR_perind, HL_perind)

  names(Output) <- c("Ho_perpop", "Ho_perloc", "He_perpop", "He_perloc", "PHt", "Hs_exp", "Hs_obs", "IR", "HL")

  # Set list of possible statistics
  Stat <- c("Ho", "He", "PHt", "Hs_exp", "Hs_obs", "IR", "HL")
  Stat_idx <- c(1,1,2,2,3,4,5,6,7)

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
