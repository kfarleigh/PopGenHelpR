# A function to estimate three measures of genetic differentiation using geno files, vcf files, or vcfR objects. Data is assumed to be bi-allelic.

A function to estimate three measures of genetic differentiation using
geno files, vcf files, or vcfR objects. Data is assumed to be
bi-allelic.

## Usage

``` r
Differentiation(
  data,
  pops,
  statistic = "all",
  missing_value = NA,
  write = FALSE,
  prefix = NULL,
  population_col = NULL,
  individual_col = NULL
)
```

## Arguments

- data:

  Character. String indicating the name of the vcf file, geno file or
  vcfR object to be used in the analysis.The genotypes within the vcf
  should be seperated by a "/" or "\|". This normally indicates unphased
  and phased genotypes, respectively. Please reach out to PopGenHelpR
  authors if you have questions.

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

- statistic:

  Character. String or vector indicating the statistic to calculate.
  Options are any of: all; all of the statistics; Fst, Weir and
  Cockerham (1984) Fst; NeisD, Nei's D statistic; JostsD, Jost's D.

- missing_value:

  Character. String indicating missing data in the input data. It is
  assumed to be NA, but that may not be true (is likely not) in the case
  of geno files.

- write:

  Boolean. Whether or not to write the output to files in the current
  working directory. There will be one or two files for each statistic.
  Files will be named based on their statistic such as Fst_perpop.csv.

- prefix:

  Character. Optional argument. String that will be appended to file
  output. Please provide a prefix if write is set to TRUE.

- population_col:

  Numeric. Optional argument (a number) indicating the column that
  contains the population assignment information.

- individual_col:

  Numeric. Optional argument (a number) indicating the column that
  contains the individuals (i.e., sample name) in the data.

## Value

A list containing the estimated heterozygosity statistics. The per pop
values are calculated by taking the average of the per locus estimates.

## References

**Fst:**

Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013). StAMPP: An R
package for calculation of genetic differentiation and structure of
mixed‐ploidy level populations. Molecular ecology resources, 13(5),
946-952.[doi:10.1111/1755-0998.12129](https://doi.org/10.1111/1755-0998.12129)

[Weir, B. S., & Cockerham, C. C.
(1984)](https://www.jstor.org/stable/2408641). Estimating F-statistics
for the analysis of population structure. evolution, 1358-1370.

**Nei's D:**

Nei, M. (1972). Genetic distance between populations. The American
Naturalist, 106(949),
283-292.[doi:10.1086/282771](https://doi.org/10.1086/282771)

[doi:10.1111/1755-0998.12129](https://doi.org/10.1111/1755-0998.12129)
Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013). StAMPP: An R
package for calculation of genetic differentiation and structure of
mixed‐ploidy level populations. Molecular ecology resources, 13(5),
946-952.

**Jost's D:**

Jost L (2008). GST and its relatives do not measure differentiation.
Molecular Ecology, 17,
4015–4026.[doi:10.1111/j.1365-294X.2008.03887.x](https://doi.org/10.1111/j.1365-294X.2008.03887.x)

## Author

Keaka Farleigh

## Examples

``` r
# \donttest{
data("HornedLizard_Pop")
data("HornedLizard_VCF")
Test <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)# }
#> Loading required namespace: vcfR
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations
```
