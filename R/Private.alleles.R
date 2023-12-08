#' A function to estimate heterozygosity.
#'
#' @param VCF Character string indicating the name of the vcf file to be used in analysis.
#' @param pops Character string indicating the name of the population assignment file. This file should have four columns and be in the same order as your vcf file. The first column named Sample indicates the sample name. The second column named Population indicates the population assignment of each individual. The third column named Longitude indicates the longitude of the sample.  The fourth column named Latitude indicates the latitude of the sample.
#' @param ploidy Numeric. The ploidy of the data.
#' @param prefix Character string that will be appended to file output.
#' @param write Boolean. Whether or not to write the output to a file in the current working directory.
#'
#' @return A list containing the estimated diversity statistics, model output from linear regression of these statistics against latitude, and model plots.
#' @export
#'
#' @examples
#' \donttest{
#' data("HornedLizard_Pop")
#' data("HornedLizard_VCF")
#' Test <- Div_stats(VCF = HornedLizard_VCF, pops = HornedLizard_Pop,
#' ploidy = 2, write = FALSE)}
Private.alleles <- function(data, pops, ploidy, statistic, write = FALSE, prefix, population_col = NULL, individual_col = NULL) {
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
    gt <- vcfR::extract.gt(Dat, return.alleles = TRUE)
    Dat <- as.data.frame(t(as.matrix(gt)))
    # Preserve individual names
    Inds <- rownames(Dat)
  }
  else if(tools::file_ext(data) == 'vcf') {
    Dat <- vcfR::read.vcfR(data, verbose = FALSE)
    print("VCF file detected, proceeding to formatting.")
    # Convert the vcf gt slot to a geno style table for calculations
    gt <- vcfR::extract.gt(Dat, return.alleles = TRUE)
    Dat <- as.data.frame(t(as.matrix(gt)))
    # Preserve individual names
    Inds <- rownames(Dat)
  }
 } else {
  stop("Please supply a vcf file or vcfR object for analysis")
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

# Get locus names
locnames <- rownames(PA_perpop[[1]])

PA_test_df <- do.call("rbind", PA_perpop)
PA_test_df$loc <- gsub('^.*\\.', "", rownames(PA_test_df))

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
Al_1 <- suppressMessages(anti_join(P_test[which(P_test$loc == locnames[j]),c(1,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(1,4)]))
Al_2 <- suppressMessages(anti_join(P_test[which(P_test$loc == locnames[j]),c(2,4)], Rem_pop[which(Rem_pop$loc == locnames[j]),c(2,4)]))

colnames(Al_1) <- colnames(Al_2) <- c('Allele', "Locus")
PA_df <- rbind(Al_1, Al_2)

if(nrow(PA_df) != 0){
  PA[[j]] <- PA_df
    }
  }
  PA_test[[i]] <- PA
}

names(PA_test) <- P_uniq

for(i in 1:length(PA_test)) {
  PA_test[[i]] <- do.call('rbind', PA_test[[i]])
}

# Count the number of private alleles in each population
PA_count <- list()
for(i in 1:length(PA_test)) {
  PA_count[[i]] <- nrow(PA_test[[i]])
}
names(PA_count) <- names(PA_test)

PA_count <- do.call("rbind", PA_count)

colnames(PA_count) <- "Number.Private.Alleles"

##########################
##### Heterozygosity #####
##########################
# Calculate heterozygosity and standard deviation for each population
Het <- hierfstat::basic.stats(Hstat)

# Get per populations heterozygosity
H_all <- data.frame(colMeans(Het$Ho, na.rm = TRUE))
H_all$Pop <- rownames(H_all)
colnames(H_all) <- c('Heterozygosity', 'Pop')

# Standard deviation
H_all$StandardDeviation <- stats::sd(H_all$Heterozygosity)

# Attach coordinates
# First we check to make sure that the populations are in the right order
if(any(H_all$Pop != unique(P[,2]))){
  stop("Populations are not in the correct order, if you see this please email the package authors")
}

# Pull populations coordinates
Pop_coords <- P[!duplicated(P$Population),]

# Append to heterozygosity estimates
H_all[,4:5] <- Pop_coords[,3:4]
colnames(H_all) <- c('Heterozygosity', 'Pop', 'Standard.Deviation','Longitude', 'Latitude')

# Is there a statistical relationship between heterozygosity and latitude
Het_model <- stats::lm(H_all$Heterozygosity ~ H_all$Latitude)

message('Heterozygosity calculated, moving to private alleles')


##########################
##### Visualizations #####
##########################
##### Heterozygosity #####
# We will make a plot of heterozygostiy estimates and an interpolated map of heterozygosity values for each population

# Plot of heterozygosity values (y-axis) and latitude (x-axis)
Het_plot <- ggplot2::ggplot(data = H_all, ggplot2::aes(x = Latitude, y = Heterozygosity)) + ggplot2::geom_point(ggplot2::aes(color = Pop),size = 3) +
  ggplot2::geom_errorbar(data = H_all, ggplot2::aes(x = Latitude, ymin = Heterozygosity-Standard.Deviation, ymax = Heterozygosity+Standard.Deviation,
                                                    color = Pop)) +
  ggplot2::stat_smooth(method = 'lm', formula = y ~ x, size =1, colour = 'black') +
  ggplot2::labs(x = 'Latitude', y = 'Heterozygosity') +
  ggplot2::ggtitle(expression(atop("Obersved Heterozygosity of Localities"))) +
  ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                            legend.title = ggplot2::element_blank())


Output <- list(H_all, Het_model, Het_plot)
names(Output) <- c('Heterozygosity_calculations', "Heterozygosity_vs_Latitude_Model",
                   "Heterozygosity_vs_Latitude_plot")
if(write == TRUE){
  # Write out heterozygosity results
  utils::write.csv(H_all, file = paste(as.character(prefix), '_Heterozygosity.csv', sep = ''), row.names = FALSE)
  summary(Het_model)
  # Write out linear regression results heterozygosity ~ latitude
  sink(paste(as.character(prefix), '_Heterozygosity_lm.txt', sep = ''))
  summary(Het_model)
  sink()
}

message("Calculations have finished, the packages used to perform file formatting and calculations were
  vcfR, adegenet, and dartR for formatting, hierfstat to calculate heterozygosity, and poppr to calculate private alleles")

return(Output)
}
