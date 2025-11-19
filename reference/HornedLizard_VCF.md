# A vcfR object to be used in `Heterozygosity` and `Differentiation`.

Data frame containing 4 columns and 72 rows

## Usage

``` r
data(HornedLizard_Pop)
```

## Format

A vcfR object

- vcfR object:

  A vcfR object containing genotype and sample informaiton for 72
  individuals.

## Source

Farleigh, K., Vladimirova, S. A., Blair, C., Bracken, J. T., Koochekian,
N., Schield, D. R., ... & Jezkova, T. (2021). The effects of climate and
demographic history in shaping genomic variation across populations of
the Desert Horned Lizard (Phrynosoma platyrhinos). Molecular Ecology,
30(18), 4481-4496.

## Examples

``` r
 # \donttest{
data("HornedLizard_Pop")
data("HornedLizard_VCF")
Test <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)# }
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations

```
