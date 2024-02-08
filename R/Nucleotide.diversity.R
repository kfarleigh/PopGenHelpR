#' A function to estimate the number of private alleles in each population.
#'
#' @param data Character. String indicating the name of the vcf file or vcfR object to be used in the analysis.
#' @param pops Character. String indicating the name of the population assignment file or dataframe containing the population assignment information for each individual in the data. This file must be in the same order as the vcf file and include columns specifying the individual and the population that individual belongs to. The first column should contain individual names and the second column should indicate the population assignment of each individual. Alternatively, you can indicate the column containing the individual and population information using the individual_col and population_col arguments.
#' @param write Boolean. Optional argument indicating Whether or not to write the output to a file in the current working directory. This will output to files; 1) the table of private allele counts per population (named prefix_PrivateAlleles_countperpop) and 2) metadata associated with the private alleles (named prefix_PrivateAlleles_metadata). Please supply a prefix it you write files to your working directory as a best practice.
#' @param prefix Character. Optional argument indicating a string that will be appended to file output. Please set a prefix if write is TRUE.
#' @param population_col Numeric. Optional argument (a number) indicating the column that contains the population assignment information.
#' @param individual_col Numeric. Optional argument (a number) indicating the column that contains the individuals (i.e., sample name) in the data.
#' @return A list containing the count of private alleles in each population and the metadata for those alleles. The metadata is a list that contains the private allele and locus name for each population.
#' @export
#'
#' @author Keaka Farleigh
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Private.alleles(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)}
Nucleotide.diversity <- function(data, pops, missing_value = NA, write = FALSE, prefix = NULL, population_col = NULL, individual_col = NULL) {
  Pop <- Standard.Deviation <- Data <- V1 <- V2 <- loc <- NULL

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
    gt <- vcfR::extract.gt(Dat, return.alleles = TRUE)
    gt_freq <- vcfR::extract.gt(Dat)
    Dat <- as.data.frame(t(as.matrix(gt)))


    gt_freq[gt_freq == "0/0"] <- 0
    gt_freq[gt_freq == "0/1" | gt_freq == "1/0"] <- 1
    gt_freq[gt_freq == "1/1"] <- 2
    Dat_frq <- as.data.frame(t(as.matrix(gt_freq)))
    # Preserve individual names
    Inds <- rownames(Dat)
  } else if(tools::file_ext(data) == 'vcf') {
    Dat <- vcfR::read.vcfR(data, verbose = FALSE)
    print("VCF file detected, proceeding to formatting.")
    # Convert the vcf gt slot to a geno style table for calculations
    gt <- vcfR::extract.gt(Dat, return.alleles = TRUE)
    gt_freq <- vcfR::extract.gt(Dat)
    Dat <- as.data.frame(t(as.matrix(gt)))


    gt_freq[gt_freq == "0/0"] <- 0
    gt_freq[gt_freq == "0/1" | gt_freq == "1/0"] <- 1
    gt_freq[gt_freq == "1/1"] <- 2
    Dat_frq <- as.data.frame(t(as.matrix(gt_freq)))
    # Preserve individual names
    Inds <- rownames(Dat)
  } else {
    stop("Please supply a vcf file or vcfR object for analysis")
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
  }

  P <- Pops
  Dat <- cbind.data.frame(Inds, P, Dat)
  Dat_frq <- cbind.data.frame(Inds, P, Dat_frq)

  # Break into list with populations for each element
  Dat_perpop <- list()
  for(i in unique(P)){
    Dat_perpop[[i]] <- Dat[which(Dat[,2] == i),]
  }

  Dat_perpop_frq <- list()
  for(i in unique(P)){
    Dat_perpop_frq[[i]] <- Dat_frq[which(Dat_frq[,2] == i),]
  }

  message('Formatting has finished, moving onto calculations')
#############################################
##### Nucleotide diversity calculations #####
#############################################
#http://www.columbia.edu/cu/biology/courses/c3020/solutions-2.html

#Freq_a1*Freq_a2*prop_nonidentical sites

  # Dat is Dat_perpop, Dat_frq is Dat_perpop_frq
  nuc.div <- function(Dat, Dat_frq){
    Dat.pop <- Dat[,2]
    Dat.gen <- Dat[,3:ncol(Dat)]
    rownames(Dat.gen) <- Dat[,1]

  # Seperate alleles and create sequences for each individual
  Alleles <- list()
  for(i in 1:ncol(Dat.gen)){
    tmp <- Dat.gen[,i]
    tmp2 <- strsplit(tmp, split = '/')
    Alleles[[i]] <- do.call(rbind,tmp2)
    remove(tmp, tmp2)
  }
  Seqs <- do.call(cbind, Alleles)
  rownames(Seqs) <- rownames(Dat.gen)

  # Get pairwise comparisons
  Comps <- combn(rownames(Seqs), m = 2)

  # Determine how many differences there are per site per comparison
  # Use the same structure as the Comps object
  Difs <- as.data.frame(t(Comps))
  Difs$Differences <- NA

  # Perform each comparison
  for(i in 1:ncol(Comps)){
    Inds2comp <- Comps[,i]
    Comp_tmp <- Seqs[which(rownames(Seqs) %in% Inds2comp),]
    tmp_res <- as.data.frame(table((Comp_tmp[1,] == Comp_tmp[2,])))

    # Make sure that we select the right table entry
    Dif_row <- which(tmp_res[,1] == "FALSE")
    PW_difs <- as.numeric(tmp_res[Dif_row,2])

    # Store results
    if(!rlang::is_empty(PW_difs)){
    Difs[i,3] <- PW_difs
    } else{
      Difs[i,3] <- 0
    }

    remove(Inds2comp, Comp_tmp, tmp_res, Dif_row, PW_difs)

  }

  ### Calculate nucleotide diversity (pi)
  ## Get allele frequencies                     #### MIGHT NEED TO BE REWORKED #######
  nloc <- length(3:ncol(Dat_frq))
  nind <- nrow(Dat_frq)
  q.freq <- matrix(nrow = 1, ncol = nloc)
  q.freq[1,] <- (((colSums(Dat_frq[3:ncol(Dat_frq)] == 2, na.rm = T))*2) + colSums(Dat_frq[3:ncol(Dat_frq)] == 1, na.rm = T))/(2*colSums(Dat_frq[3:ncol(Dat_frq)] != "NA"))
  p.freq <- 1-q.freq

  p.freq.bar <- mean(p.freq)
  q.freq.bar <- mean(q.freq)

  Difs$Prop <- Difs$Differences/nloc
  Difs$Prop.Pq <- Difs$Prop*(p.freq.bar*q.freq.bar)

  Nucdiv <- sum(Difs$Prop.Pq)
  }

}




