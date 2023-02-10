## This is a resubmission, I have included the comments from the CRAN maintainer and my changes below. 
```
Please remove the redundant "A Package to" from your title.
```
This has been removed and the title has been changed to "Streamline Population Genomic and Genetic Analyses""
```
Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'PopGenHelpR'
Please note that package names are case sensitive.
```
The package name has been placed in single quotes in the description.
```
If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")
```
There is no paper for the package at this time.
```
\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or create
additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer
than 5 sec to be executed, but it would be preferable to have automatic
checks for functions. Otherwise, you can also write some tests.)
```
All examples have been changed from \dontrun to \donttest.
```
Please ensure that your functions do not write by default or in your
examples/vignettes/tests in the user's home filespace (including the
package directory and getwd()). This is not allowed by CRAN policies.
Please omit any default path in writing functions. In your
examples/vignettes/tests you can write to tempdir().
```
All vignette examples have been set to eval = FALSE so that they do not write anything out. Examples have been set to donttest so that they do not write into teh user's home filespace. 

## R CMD check results

0 errors | 0 warnings | 1 note
```
* This is a new release. 
```

## check_mac_release results
There were no ERRORs, WARNINGs, or NOTEs.

## check_win_release results
There was 1 NOTE:
```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Keaka Farleigh <keakafarleigh@gmail.com>'

New submission
```

## check_rhub results
```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Keaka Farleigh <keakafarleigh@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    PopGenHelpR (9:132)
    csv (9:157)
    vcf (9:149)
```
The words are spelled correctly

```
* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
  
    'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

## check_win_devel results
There was 1 NOTE:
```
* checking CRAN incoming feasibility ... [8s] NOTE
Maintainer: 'Keaka Farleigh <keakafarleigh@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Genomic (2:43)
  PopGenHelpR (9:132)
  csv (9:157)
  genomic (9:48)
  vcf (9:149)
```
The words are spelled correctly
