# A function to estimate the number of private alleles in each population.

A function to estimate the number of private alleles in each population.

## Usage

``` r
Private.alleles(
  data,
  pops,
  write = FALSE,
  prefix = NULL,
  population_col = NULL,
  individual_col = NULL
)
```

## Arguments

- data:

  Character. String indicating the name of the vcf file or vcfR object
  to be used in the analysis.

- pops:

  Character. String indicating the name of the population assignment
  file or dataframe containing the population assignment information for
  each individual in the data. This file must be in the same order as
  the vcf file and include columns specifying the individual and the
  population that individual belongs to. The first column should contain
  individual names and the second column should indicate the population
  assignment of each individual. Alternatively, you can indicate the
  column containing the individual and population information using the
  individual_col and population_col arguments.

- write:

  Boolean. Optional argument indicating Whether or not to write the
  output to a file in the current working directory. This will output to
  files; 1) the table of private allele counts per population (named
  prefix_PrivateAlleles_countperpop) and 2) metadata associated with the
  private alleles (named prefix_PrivateAlleles_metadata). Please supply
  a prefix it you write files to your working directory as a best
  practice.

- prefix:

  Character. Optional argument indicating a string that will be appended
  to file output. Please set a prefix if write is TRUE.

- population_col:

  Numeric. Optional argument (a number) indicating the column that
  contains the population assignment information.

- individual_col:

  Numeric. Optional argument (a number) indicating the column that
  contains the individuals (i.e., sample name) in the data.

## Value

A list containing the count of private alleles in each population and
the metadata for those alleles. The metadata is a list that contains the
private allele and locus name for each population.

## Author

Keaka Farleigh

## Examples

``` r
# \donttest{
data("HornedLizard_Pop")
data("HornedLizard_VCF")
Test <- Private.alleles(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)# }
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations
#> [1] "Finished private allele calculations for East"
#> [1] "Finished private allele calculations for South"
#> [1] "Finished private allele calculations for West"
```
