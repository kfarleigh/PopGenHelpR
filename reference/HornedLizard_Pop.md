# A population assignment data frame to be used in `Heterozygosity` and `Differentiation`.

Data frame containing 4 columns and 72 rows

## Usage

``` r
data(HornedLizard_Pop)
```

## Format

A data frame with 4 columns and 72 rows:

- Sample:

  Sample Name

- Population:

  Population assignment according to sNMF results (see citation)

- Longitude:

  Longitude

- Latitude:

  Latitude

## Source

Coordinates and population names taken from Farleigh, K., Vladimirova,
S. A., Blair, C., Bracken, J. T., Koochekian, N., Schield, D. R., ... &
Jezkova, T. (2021). The effects of climate and demographic history in
shaping genomic variation across populations of the Desert Horned Lizard
(Phrynosoma platyrhinos). Molecular Ecology, 30(18), 4481-4496.

## Examples

``` r
 # \donttest{
data("HornedLizard_Pop")
data("HornedLizard_VCF")
Test <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, write = FALSE)# }
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations

```
