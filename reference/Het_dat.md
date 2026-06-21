# A data frame of hypothetical heterozygosity data produced by Heterozygosity.

Data frame containing 5 columns and 3 rows

## Usage

``` r
data(Het_dat)
```

## Format

A data frame with 5 columns and 3 rows:

- Heterozygosity:

  Estimated heterozygosity

- Pop:

  Population assignment

- Standard.Deviation:

  standard deviation

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
data(Het_dat)
Test <- Point_map(Het_dat, statistic = "Heterozygosity")# }
#> Cached as: /tmp/Rtmpb9RhkE/gadm/gadm36_adm0_r5_pk.rds
```
