#' A function to estimate seven measures of heterozygosity using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.
#'
#' @param data Character. String indicating the name of the vcf file, geno file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param statistic Character. String or vector indicating the statistic to calculate. Options are any of: all; all of the statistics; Ho, observed heterozygosity; He, expected heterozygosity; PHt, proportion of heterozygous loci; StHe, heterozygosity standardized by the average expected heterozygosity; StHo, heterozygosity standardized by the average observed heterozygosity; IR, internal relatedness; HL, homozygosity by locus.
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
#' \href{https://royalsocietypublishing.org/doi/abs/10.1098/rspb.2001.1751?casa_token=aXK15yYNjpEAAAAA%3ATR6IXyuIl03yygS_CWET4ymKFthc0dxwkpgCRGGZL250Am_Ssbxi9bnTDIQmpNnshM7H6vTZ1v83PD0}{Amos, W., Worthington Wilmer, J., Fullard, K., Burg, T. M., Croxall, J. P., Bloch, D., & Coulson, T. (2001)}. The influence of parental relatedness on reproductive success. Proceedings of the Royal Society of London. Series B: Biological Sciences, 268(1480), 2021-2027.
#'
#' \href{https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2006.03111.x}{Aparicio, J. M., Ortego, J., & Cordero, P. J. (2006)}. What should we weigh to estimate heterozygosity, alleles or loci?. Molecular Ecology, 15(14), 4659-4665.
#'
#' \bold{Heterozygosity standardized by expected (StHe) and observed heterozygosity (StHo):}
#'
#' \href{https://academic.oup.com/evolut/article/53/4/1259/6757148}{Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J. M. (1999)}. Parasite‐mediated selection against Inbred Soay sheep in a free‐living island populaton. Evolution, 53(4), 1259-1267.
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
  ##### Heterozygosity calculations #####
  #######################################
 # Estimate observed heterozygosity (1-homozygosity), accounts for NA
 ObsHet <- function(Dat){
    ObsHet_perloc <- 1-(colSums(Dat[3:ncol(Dat)] != 1, na.rm = T)/(nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)]))))
    return(ObsHet_perloc)
  }

  ObsHet_res_perloc  <- lapply(Dat_perpop, ObsHet)
  Obs_Het_res_avg <- lapply(ObsHet_res_perloc, mean)

  # Make it into a dataframe so that it is easier for users to view
  ObsHet_res_perloc <- mapply(cbind, ObsHet_res_perloc, "Pop"=names(ObsHet_res_perloc), SIMPLIFY=F)
  ObsHet_res_perloc <- as.data.frame(do.call("rbind", ObsHet_res_perloc))
  colnames(ObsHet_res_perloc)[1] <- "Observed.Heterozygosity"
  ObsHet_res_perloc[,1] <- as.numeric(ObsHet_res_perloc[,1])

  Obs_Het_res_avg <- as.data.frame(do.call("rbind", Obs_Het_res_avg))
  Obs_Het_res_avg$Pop <- rownames(Obs_Het_res_avg)
  colnames(Obs_Het_res_avg)[1] <- "Observed.Heterozygosity"
  Obs_Het_res_avg[,1] <- as.numeric(Obs_Het_res_avg[,1])

  # Estimate expected heterozygosity
 ExpHe <- function(Dat){
  p <- ((colSums(Dat[3:ncol(Dat)]== 0, na.rm = T)*2) + colSums(Dat[3:ncol(Dat)]== 1, na.rm = T))/((nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)])))*2)
  q <- ((colSums(Dat[3:ncol(Dat)]== 2, na.rm = T)*2) + colSums(Dat[3:ncol(Dat)]== 1, na.rm = T))/((nrow(Dat)- colSums(is.na(Dat[3:ncol(Dat)])))*2)
   He_perloc <- 1-(p^2)-(q^2)
   return(He_perloc)
  }

  ExpHet_res_perloc  <- lapply(Dat_perpop, ExpHe)
  ExpHet_res_avg  <- lapply(ExpHet_res_perloc, mean)

  # Make it into a dataframe so that it is easier for users to view
  ExpHet_res_perloc <- mapply(cbind, ExpHet_res_perloc, "Pop"=names(ExpHet_res_perloc), SIMPLIFY=F)
  ExpHet_res_perloc <- as.data.frame(do.call("rbind", ExpHet_res_perloc))
  colnames(ExpHet_res_perloc)[1] <- "Expected.Heterozygosity"
  ExpHet_res_perloc[,1] <- as.numeric(ExpHet_res_perloc[,1])

  ExpHet_res_avg <- as.data.frame(do.call("rbind", ExpHet_res_avg))
  ExpHet_res_avg$Pop <- rownames(ExpHet_res_avg)
  colnames(ExpHet_res_avg)[1] <- "Expected.Heterozygosity"
  ExpHet_res_avg[,1] <- as.numeric(ExpHet_res_avg[,1])

  # Estimate the proportion of heterozygous loci per individual
 PropHt <- function(Dat){
   PHt <- t(as.data.frame(rowSums(Dat[3:ncol(Dat)] == 1, na.rm = T)/((ncol(Dat)-2) - rowSums(is.na(Dat[3:ncol(Dat)])))))
   rownames(PHt) <- 'PHt'
   colnames(PHt) <- Dat[,1]

   return(PHt)
 }

 PropHt_res_perind <- lapply(Dat_perpop, PropHt)
 PropHt_res_perind <- mapply(rbind, PropHt_res_perind, "Pop"=names(PropHt_res_perind), SIMPLIFY=F)

 # Make it into a dataframe so that it is easier for users to view
 PropHt_res_perind <- t(as.data.frame(do.call("cbind", PropHt_res_perind)))
 PropHt_res_perind <- as.data.frame(PropHt_res_perind)
 PropHt_res_perind[,1] <- as.numeric(PropHt_res_perind[,1])

 # Estimate standardized heterozygosity based on the expected heterozygosity
 StHe <- function(Dat){
   St_He <- PropHt(Dat)/mean(ExpHe(Dat))
   rownames(St_He) <- 'StHe'
   return(St_He)
 }

 StHe_res_perind <- lapply(Dat_perpop, StHe)
 StHe_res_perind <- mapply(rbind, StHe_res_perind, "Pop"=names(StHe_res_perind), SIMPLIFY=F)

 # Make it into a dataframe so that it is easier for users to view
 StHe_res_perind <- t(as.data.frame(do.call("cbind", StHe_res_perind)))
 StHe_res_perind <- as.data.frame(StHe_res_perind)
 StHe_res_perind[,1] <- as.numeric(StHe_res_perind[,1])

 # Estimate standardized heterozygosity based on the observed heterozyosity
 StHo <- function(Dat){
   St_Ho <- PropHt(Dat)/mean(ObsHet(Dat))
   rownames(St_Ho) <- 'StHo'
   return(St_Ho)
 }

 StHo_res_perind <- lapply(Dat_perpop, StHo)
 StHo_res_perind <- mapply(rbind, StHo_res_perind, "Pop"=names(StHo_res_perind), SIMPLIFY=F)

 # Make it into a dataframe so that it is easier for users to view
 StHo_res_perind <- t(as.data.frame(do.call("cbind", StHo_res_perind)))
 StHo_res_perind <- as.data.frame(StHo_res_perind)
 StHo_res_perind[,1] <- as.numeric(StHo_res_perind[,1])


 IR <- function(Dat){
     # Extract genetic data and convert to count the frequency of allele s
     tmp <- Dat[,3:ncol(Dat)]
     tmp[tmp == "0"] <- "0/0"
     tmp[tmp == "1"] <- "0/2"
     tmp[tmp == "2"] <- "2/2"

     rownames(tmp) <- Inds

     Loc_IR_formatted <- data.frame(matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp)*2), row.names = rownames(tmp))

     for(i in 1:ncol(tmp)){
       # Get the locus name
       Loc_nam <- colnames(tmp)[i]

       # Split the genotype
       Loc_split <- strsplit(tmp[,i], "/", fixed = T)
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

     # Get the number of individuals
     Individuals <- nrow(Loc_IR_mat)

     # Get the number of loci
     Nloc <- ncol(Loc_IR_mat)/2

     # Set up a results table
     res_tab <- data.frame(IR = matrix(NA, nrow = Individuals, ncol = 1), row.names = Inds)

     # Get the counts of alleles for each locus
     Counts <- list()

     # Count the occurrences of each allele (0 and 2) at each locus
     for(i in 1:Nloc) {
       # Set the same index as above
       idx <- i*2-1
       Counts[[i]] <- table(Loc_IR_mat[,idx:(idx+1)])
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

         if((!is.na(Loc_IR_mat[i,idx1])) && (!is.na(Loc_IR_mat[i,idx2]))){
           N <- N + 1
         }

         if(Loc_IR_mat[i,idx1] == Loc_IR_mat[i,idx2]){
           H <- H + 1
           # Which allele is the individual homozygous for
           Hom_Allele <- as.character(Loc_IR_mat[i,idx1])
           f <- f + (2 * Counts[[j]][[Hom_Allele]] - 2)/(sum(Counts[[j]]) - 2)

         } else if(Loc_IR_mat[i,idx1] != Loc_IR_mat[i,idx2]){
           # If they are heterozygous that means that they are contributing two alleles, using just a sum of 1 leads to NaN
           Het_Allele1 <- as.character(Loc_IR_mat[i,idx1])
           f <- f + (Counts[[j]][[Het_Allele1]] - 1)/(sum(Counts[[j]]) - 2)
           Het_Allele2 <- as.character(Loc_IR_mat[i,idx2])
           f <- f + (Counts[[j]][[Het_Allele2]] - 1)/(sum(Counts[[j]]) - 2)
         }
       }
       # Calculate internal relatedness
       res_tab[i,1] <- (2 * H - f) / (2 * N - f)
     }
     return(res_tab)
 }

 IR_perind <- IR(Dat)

 # Estimate homozygosity by locus
 HL <- function(Dat){
   # Extract genetic data and convert to count the frequency of allele s
   tmp <- Dat[,3:ncol(Dat)]
   tmp[tmp == "0"] <- "0/0"
   tmp[tmp == "1"] <- "0/2"
   tmp[tmp == "2"] <- "2/2"

   rownames(tmp) <- Inds

   Loc_HL_formatted <- data.frame(matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp)*2), row.names = rownames(tmp))

   for(i in 1:ncol(tmp)){
     # Get the locus name
     Loc_nam <- colnames(tmp)[i]

     # Split the genotype
     Loc_split <- strsplit(tmp[,i], "/", fixed = T)
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

   # Get the number of individuals
   Individuals <- nrow(Loc_HL_mat)

   # Get the number of loci
   Nloc <- ncol(Loc_HL_mat)/2

   # Set up a results table
   res_tab <- data.frame(HL = matrix(NA, nrow = Individuals, ncol = 1), row.names = Inds)

   E <- array(Nloc)
   Counts <- array(Nloc)

   # Count the occurrences of each allele (0 and 2) at each locus
   for(i in 1:Nloc) {
      E[i] <- 1
      # Set the same index as above
      idx <- i*2-1
      Counts[i] <- list(table(Loc_HL_mat[,idx:(idx+1)]))
      E[i] <- 1 - sum((Counts[[i]] / sum(Counts[[i]]))^2)
   }
  # Formula from the archived Rhh package https://cran.r-project.org/web/packages/Rhh/index.html
   for (i in 1:Individuals) {

     sum.Eh <- 0
     sum.Ej <- 0

     for (j in 1:Nloc) {

       idx1 <- 2*j-1
       idx2 <- 2*j

       if ((!is.na(Loc_HL_mat[i, idx1])) && (!is.na(Loc_HL_mat[i, idx2]))) {
         if (Loc_HL_mat[i, idx1] == Loc_HL_mat[i, idx2]) {
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

 HL_perind <- HL(Dat)


 Output <- list(Obs_Het_res_avg, ObsHet_res_perloc, ExpHet_res_avg, ExpHet_res_perloc,
                PropHt_res_perind, StHe_res_perind, StHo_res_perind, IR_perind, HL_perind)

 names(Output) <- c("Ho_perpop", "Ho_perloc", "He_perpop", "He_perloc", "PHt", "StHe", "StHo", "IR", "HL")

 # Set list of possible statistics
 Stat <- c("Ho", "He", "PHt", "StHe", "StHo", "IR", "HL")
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
