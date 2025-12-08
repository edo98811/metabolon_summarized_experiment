<!-- badges: start -->
[![R-CMD-check](https://github.com/edo98811/metabolon_summarized_experiment/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edo98811/helper-functions-rbioinfo/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Introduction

This is a package developed as interface from standard bioinformatics pipelines to the metabolon software.

It can load metabolon data in a `SummarizedExperiment` object, which is a common data structure used in bioinformatics for storing high-throughput assay data along with associated metadata.

Metabolon offers a platform for multiomics analysis. This package also offers an interface to it by allowing export of data in the required format. There are functions to export data from a summarized experiment as "counts" file that can be uploaded to the metabolon platform.

Furthermore, results from different standard bioinfomratics pipelines, for example differential analysis using `DESeq2` or `limma`, can be exported in the required format to be uploaded to the metabolon platform. The results can also be just a data frame with the required columns. The obrained tables can be directly uploaded to the myMetabolon platform.

# Installation 

[https://github.com/edo98811/metabolon_summarized_experiment.git](https://github.com/edo98811/metabolon_summarized_experiment.git)

To install this package, you can use the `remotes` package in R:

```r
# Install the remotes package if you haven't already
install.packages("remotes")

# Use remotes to install the package from GitHub
remotes::install_github("edo98811/metabolon_summarized_experiment")
```
