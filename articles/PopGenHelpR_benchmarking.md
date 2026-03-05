# Evaluating PopGenHelpR with adegenet, hierfstat, mmod, and StAMPP

## Purpose

We will compare the performance of `PopGenHelpR` with other R packages
available on CRAN. Below, we list the packages we will compare
`PopGenHelpR` with and the statistics in each comparison.

**F_(st) and Nei’s D**

- [StAMPP](https://cran.r-project.org/web/packages/StAMPP/index.html)
  ([Pembleton et al., 2013](https://doi.org/10.1111/1755-0998.12129))

**Jost’s D**

- [mmod](https://cran.r-project.org/web/packages/mmod/index.html)
  (Winter et al., 2017)

**Expected and Observed Heterozygosity**

- [hierfstat](https://cran.r-project.org/web/packages/hierfstat/index.html)
  ([Goudet, 2005](https://doi.org/10.1111/j.1471-8286.2004.00828.x))

- [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html)
  ([Jombart, 2008](https://doi.org/10.1093/bioinformatics/btn129))

## Notes

We also performed additional testing using external data to ensure that
PopGenHelpR works with real world data. You can find those tests at
<https://github.com/kfarleigh/PopGenHelpR_testing>.

### Let’s Begin

First we will load the packages and the data in `PopGenHelpR`. The data
comes from [Farleigh et al. (2021)](https://doi.org/10.1111/mec.16070).
We also load `vcfR` to convert between data formats ([Knaus & Grunwald,
2017](https://doi.org/10.1111/1755-0998.12549)).

``` r
# Load the packages
library(PopGenHelpR)
library(adegenet)
#> Loading required package: ade4
#> 
#>    /// adegenet 2.1.11 is loaded ////////////
#> 
#>    > overview: '?adegenet'
#>    > tutorials/doc/questions: 'adegenetWeb()' 
#>    > bug reports/feature requests: adegenetIssues()
library(hierfstat)
#> 
#> Attaching package: 'hierfstat'
#> The following objects are masked from 'package:adegenet':
#> 
#>     Hs, read.fstat
library(StAMPP)
#> Loading required package: pegas
#> Loading required package: ape
#> 
#> Attaching package: 'ape'
#> The following objects are masked from 'package:hierfstat':
#> 
#>     pcoa, varcomp
#> Registered S3 method overwritten by 'pegas':
#>   method      from
#>   print.amova ade4
#> 
#> Attaching package: 'pegas'
#> The following object is masked from 'package:ape':
#> 
#>     mst
#> The following object is masked from 'package:ade4':
#> 
#>     amova
library(mmod)
library(vcfR)
#> 
#>    *****       ***   vcfR   ***       *****
#>    This is vcfR 1.16.0 
#>      browseVignettes('vcfR') # Documentation
#>      citation('vcfR') # Citation
#>    *****       *****      *****       *****
#> 
#> Attaching package: 'vcfR'
#> The following objects are masked from 'package:pegas':
#> 
#>     getINFO, write.vcf

# Load the data 
data("HornedLizard_VCF")
data("HornedLizard_Pop")
```

## *F_(ST)* and *Nei’s D* Comparison

We will compare `PopGenHelpR` and `StAMPP`. Both packages use the
formulas from Weir and Cockerham (1984) ane Nei (1972) to calculate
*F_(ST)* and *Nei’s D*, respectively but we want to make sure that our
estimates are consistent across packages.

First, we need to format the data for `StAMPP`

Now we can calculate our statistics. Let’s start with *F_(ST)*.

``` r
PGH_fst <- Differentiation(dat = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Fst")
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations

Stmp_fst <- stamppFst(Glight, nboots = 0)
```

Let’s inspect the results.

``` r
PGH_fst$Fst
#>            East     South West
#> East         NA        NA   NA
#> South 0.2511135        NA   NA
#> West  0.3905512 0.3029886   NA
Stmp_fst
#>            East     South West
#> East         NA        NA   NA
#> South 0.2511135        NA   NA
#> West  0.3905512 0.3029886   NA

# Is there a difference between the two?
Fst_comparison <- PGH_fst$Fst-Stmp_fst

summary(Fst_comparison)
#>       East       South        West    
#>  Min.   :0   Min.   :0   Min.   : NA  
#>  1st Qu.:0   1st Qu.:0   1st Qu.: NA  
#>  Median :0   Median :0   Median : NA  
#>  Mean   :0   Mean   :0   Mean   :NaN  
#>  3rd Qu.:0   3rd Qu.:0   3rd Qu.: NA  
#>  Max.   :0   Max.   :0   Max.   : NA  
#>  NA's   :1   NA's   :2   NA's   :3
```

Now we move onto *Nei’s D*. We can use the same genlight that we created
for the *F_(ST)* calculations. We will calculate *Nei’s D* for both
population’s and individual’s.

``` r
PGH_ND <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "NeisD")
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations

# StAMPP population Nei's D
Stmp_popND <- stamppNeisD(Glight)

# StAMPP individual Nei's D
Stmp_indND <- stamppNeisD(Glight, pop = FALSE)
```

Compare the results like we did with *F_(ST)*. Note that `PopGenHelpR`
only reports result on the lower triangular element so we will set the
upper triangular element of the `Stmp_popND` and `Stmp_indND` objects to
NA.

``` r
# Population comparison
PGH_ND$NeisD_pop
#>             East     South West
#> East  0.00000000        NA   NA
#> South 0.09005846 0.0000000   NA
#> West  0.19806009 0.1148848    0
Stmp_popND
#>           [,1]     [,2]     [,3]
#> East  0.000000 0.090058 0.198060
#> South 0.090058 0.000000 0.114885
#> West  0.198060 0.114885 0.000000

# Set StAMPP upper diagnoals to NA
Stmp_popND[upper.tri(Stmp_popND)] <- NA
Stmp_indND[upper.tri(Stmp_indND)] <- NA

popND_comparison <- PGH_ND$NeisD_pop-Stmp_popND
summary(popND_comparison)
#>       East               South               West  
#>  Min.   :0.000e+00   Min.   :-1.9e-07   Min.   :0  
#>  1st Qu.:4.518e-08   1st Qu.:-1.4e-07   1st Qu.:0  
#>  Median :9.036e-08   Median :-9.0e-08   Median :0  
#>  Mean   :1.840e-07   Mean   :-9.0e-08   Mean   :0  
#>  3rd Qu.:2.760e-07   3rd Qu.:-5.0e-08   3rd Qu.:0  
#>  Max.   :4.615e-07   Max.   : 0.0e+00   Max.   :0  
#>                      NA's   :1          NA's   :2

# Get the mean difference
mean(popND_comparison, na.rm = T)
#> [1] 6.038476e-08

# Individual comparison, uncomment if you want to see it
#PGH_ND$NeisD_ind
#Stmp_indND

indND_comparison <- PGH_ND$NeisD_ind - Stmp_indND
mean(indND_comparison, na.rm = T)
#> [1] 5.807753e-09
```

We see that the difference is very small and is because of rounding.
Let’s move onto *Jost’s D* comparison with `mmod`.

## *Jost’s D* Comparison

We will compare `PopGenHelpR` and `mmod`. Both packages use the formulas
from Jost (2008). `mmod` uses genind objects so we have to do some
format conversion first.

``` r
Genind <- vcfR2genind(HornedLizard_VCF)
Genind@pop <- as.factor(HornedLizard_Pop$Population)
ploidy(Genind) <- 2

# Calculate Jost's D
PGH_JD <- Differentiation(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "JostsD")
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations
mmod_JD <- pairwise_D(Genind)

PGH_JD$JostsD
#>             East     South West
#> East  0.00000000        NA   NA
#> South 0.08135043 0.0000000   NA
#> West  0.17370519 0.1044025    0
mmod_JD
#>             East      South
#> South 0.08135043           
#> West  0.17370519 0.10440251

# Compare differences mathematically
PGH_JD$JostsD[2,1] - mmod_JD[1]
#> [1] -6.938894e-17
PGH_JD$JostsD[3,1] - mmod_JD[2]
#> [1] 0
PGH_JD$JostsD[3,2] - mmod_JD[3]
#> [1] 6.938894e-17
```

Estimate’s are very similar between `PopGenHelpR` and `mmod`, we will
move onto heterozygosity.

## *Expected* and *Observed Heterozygosity* Comparison

We will compare `PopGenHelpR`, `hierfstat`, and `adegenet`. Again, the
packages use the same formula’s, so we expect similar if not identical
results. `hierfstat` uses it’s own format, so we will convert the data
before calculations. Luckily we can convert the genind object from
*Jost’s D* comparisons.

``` r
Hstat <- genind2hierfstat(Genind)

### Calculate heterozygosities
# Expected
PGH_He <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "He")
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations
# Observed
PGH_Ho <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "Ho")
#> [1] "vcfR object detected, proceeding to formatting."
#> Formatting has finished, moving onto calculations

Hstat_hets <- basic.stats(Hstat)
Hstat_Ho <- colMeans(Hstat_hets$Ho)

He_adnet <- Hs(Genind)

PGH_He$He_perpop$Expected.Heterozygosity-He_adnet
#>         East        South         West 
#> -0.006854206 -0.003143262 -0.006625364
PGH_Ho$Ho_perpop$Observed.Heterozygosity-Hstat_Ho
#>          East         South          West 
#>  7.511576e-06 -1.790920e-06  0.000000e+00
```

We again see very small differences between the estimates.

Please reach out to Keaka Farleigh (<farleik@miamioh.edu>) if you have
any questions, and please see the references and acknowledgments below.

## References

Farleigh, K., Vladimirova, S. A., Blair, C., Bracken, J. T., Koochekian,
N., Schield, D. R., … & Jezkova, T. (2021). The effects of climate and
demographic history in shaping genomic variation across populations of
the Desert Horned Lizard (Phrynosoma platyrhinos). Molecular Ecology,
30(18), 4481-4496.

Goudet, J. (2005). hierfstat, a package for R to compute and test
hierarchical F‐statistics. Molecular ecology notes, 5(1), 184-186.

Jost, L. (2008). GST and its relatives do not measure differentiation.
Molecular ecology, 17(18), 4015-4026.

Knaus, B. J., & Grünwald, N. J. (2017). vcfr: a package to manipulate
and visualize variant call format data in R. Molecular ecology
resources, 17(1), 44-53.

Nei, M. (1972). Genetic distance between populations. The American
Naturalist, 106(949), 283-292.

Pembleton, L. W., Cogan, N. O., & Forster, J. W. (2013). St AMPP: An R
package for calculation of genetic differentiation and structure of
mixed‐ploidy level populations. Molecular ecology resources, 13(5),
946-952.

Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
analysis of population structure. evolution, 1358-1370.

Winter, D., Green, P., Kamvar, Z., & Gosselin, T. (2017). mmod: modern
measures of population differentiation (Version 1.3.3).

## Acknowledgements

We thank the authors of `hierfstat`, `mmod`, `StAMPP`, and all of the
package dependencies. They provided inspiration for `PopGenHelpR` and
their commitment to open science made it possible to develop and
benchmark our package.
