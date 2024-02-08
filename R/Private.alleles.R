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
Private.alleles <- function(data, pops, write = FALSE, prefix = NULL, population_col = NULL, individual_col = NULL) {
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
    Dat <- as.data.frame(t(as.matrix(gt)))
    # Preserve individual names
    Inds <- rownames(Dat)
  } else if(tools::file_ext(data) == 'vcf') {
    Dat <- vcfR::read.vcfR(data, verbose = FALSE)
    print("VCF file detected, proceeding to formatting.")
    # Convert the vcf gt slot to a geno style table for calculations
    gt <- vcfR::extract.gt(Dat, return.alleles = TRUE)
    Dat <- as.data.frame(t(as.matrix(gt)))
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

# Break into list with populations for each element
Dat_perpop <- list()
for(i in unique(P)){
  Dat_perpop[[i]] <- Dat[which(Dat[,2] == i),]
}

message('Formatting has finished, moving onto calculations')
#######################################
##### Private allele calculations #####
#######################################


Uniq_alleles <- function(x) {
  tmp_res <- list()
  for(i in 3:ncol(x)){
    tmp <- unique(unlist(strsplit(x[,i], split = '/')))
    tmp_res[[i]] <- tmp
    remove(tmp)
  }
  tmp_res <- as.data.frame(do.call("rbind", tmp_res))
  rownames(tmp_res) <- colnames(x[3:ncol(x)])
  return(tmp_res)
}

PA_perpop <- lapply(Dat_perpop, Uniq_alleles)
PA_perpop <- mapply(cbind, PA_perpop, "Pop"=names(PA_perpop), SIMPLIFY=F)


PA_test_df <- do.call("rbind", PA_perpop)
PA_test_df$loc <- gsub('^.*\\.', "", rownames(PA_test_df))

# Get locus names
locnames <- unique(PA_test_df$loc)

P_uniq <- unique(P)

PA_test <- list()
for(i in 1:length(P_uniq)){
  PA <- list()
  for(j in 1:length(locnames)){
# Isolate test population
P_test <- PA_test_df[which(PA_test_df$Pop == P_uniq[i]),]
# Isolate remaining populations
Rem_pop <- PA_test_df[-c(which(PA_test_df$Pop == P_uniq[i])),]
# For each allele, find any alleles that are in the population P_test but not the remaining populations (Rem_pop)
Al_11 <- suppressMessages(dplyr::anti_join(P_test[which(P_test$loc == locnames[j]),c(1,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(1,4)]))
Al_22 <- suppressMessages(dplyr::anti_join(P_test[which(P_test$loc == locnames[j]),c(2,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(2,4)]))
Al_21 <- suppressMessages(dplyr::anti_join(P_test[which(P_test$loc == locnames[j]),c(2,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(1,4)], by = dplyr::join_by(V2 == V1, loc)))
Al_12 <- suppressMessages(dplyr::anti_join(P_test[which(P_test$loc == locnames[j]),c(1,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(2,4)], by = dplyr::join_by(V1 == V2, loc)))

# Check for cross comparisons, set to be empty if either allele was erroneously identified as private
# For example if we have A/G in P_test and G/A in Rem_pop; Al_11 and Al_22 would identify them as private because A != G and G !=A
# But they are indeed not private
if(nrow(Al_12) == 0 | nrow(Al_11) == 0){
  Al_11 <- Al_12 <- data.frame(Allele = character(), Locus = character())
}
if(nrow(Al_21) == 0 | nrow(Al_22) == 0){
  Al_22 <- Al_21 <- data.frame(Allele = character(), Locus = character())
}

colnames(Al_11) <- colnames(Al_22) <- colnames(Al_21) <- colnames(Al_12) <-c('Allele', "Locus")
PA_df <- rbind(Al_11, Al_22, Al_21, Al_12)

# Only keep unique private alleles
PA_df <- PA_df %>% dplyr::distinct()

if(nrow(PA_df) != 0){
  PA[[j]] <- PA_df
    }
  }
  if(rlang::is_empty(PA)){
    PA_test[[i]] <- list()
  } else{
  PA_test[[i]] <- PA
  }
  print(paste("Finished private allele calculations for ", P_uniq[i], sep = ""))
}

names(PA_test) <- P_uniq

for(i in 1:length(PA_test)) {
  if(length(PA_test[[i]]) > 0){
  PA_test[[i]] <- do.call('rbind', PA_test[[i]])
  } else{
    PA_test[[i]] <- data.frame(Allele = character(), Locus = character())
  }
}

# Count the number of private alleles in each population
PA_count <- list()
for(i in 1:length(PA_test)) {
  PA_count[[i]] <- nrow(PA_test[[i]])
}
names(PA_count) <- names(PA_test)

PA_count <- do.call("rbind", PA_count)

PA_count <- as.data.frame(PA_count)
PA_count$Standard.deviation <- stats::sd(PA_count[,1])
colnames(PA_count) <- c("Number.Private.Alleles", "Standard.deviation")

if(write == TRUE && !is.null(prefix)){
 PA_test2 <- mapply(cbind, PA_test, "Pop"=names(PA_perpop), SIMPLIFY=F)
 PA_mdata <- do.call("rbind", PA_test2)
 utils::write.csv(PA_mdata, paste(prefix, "PrivateAlleles_metadata.csv", sep = '_'), row.names = F)
 utils::write.csv(PA_count, paste(prefix, "PrivateAlleles_countperpop.csv", sep = '_'), row.names = F)
} else if(write == TRUE && is.null(prefix)) {
  PA_test2 <- mapply(cbind, PA_test, "Pop"=names(PA_perpop), SIMPLIFY=F)
  PA_mdata <- do.call("rbind", PA_test2)
  utils::write.csv(PA_mdata, "PrivateAlleles_metadata.csv", row.names = F)
  utils::write.csv(PA_count, "PrivateAlleles_countperpop.csv", row.names = F)
}

Final_res <- list(PA_count, PA_test)
names(Final_res) <- c("Private.Allele.Count", "Private.Allele.Metadata")

return(Final_res)

}

