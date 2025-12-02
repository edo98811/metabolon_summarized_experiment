# metabolonIO

## setup

``` r
library("metabolonIO")
```

``` r
suppressPackageStartupMessages({
  library("macrophage")
  library("SummarizedExperiment")
  library("limma")
  library("edgeR")
  library("DESeq2")
})
```

## Path to example client data table.

This is data saved in the package now, it will be needed later. The cdt
path is th e path to the Client data table, in this case we took the
model of a client data table and filled it up with example data

``` r
cdt_path <- system.file("extdata", "example_data_tables.XLSX", package = "metabolonIO")

data(gse)
rownames(gse) <- substr(rownames(gse), 1, 15) # per togliere tutto quello dopo il punto
```

## Setting up example data.

First step is setting up a differential expression analysis using DESeq2
and limma on the macrophage data included in the macrophage package as
example, this data will be used to show how to export it to metabolon.

Setting up the data

``` r
dds <- DESeqDataSet(gse, design = ~condition)
#> using counts and average transcript lengths from tximeta

keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

dds <- DESeq(dds)
#> estimating size factors
#> using 'avgTxLength' from assays(dds), correcting for library size
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
```

``` r
condition <- factor(colData(gse)[, "condition_name"])

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

contrast.matrix <- makeContrasts(
  IFNg_vs_naive = IFNg - naive,
  levels = design
)

d0 <- DGEList(assay(gse))
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]

matrix_to_fit_model_to <- voom(d, design, plot = F)
# matrix_to_fit_model_to <- assay(se) # Use the assay data from the SummarizedExperiment object

fit <- lmFit(matrix_to_fit_model_to, design) #
fit2 <- contrasts.fit(fit, contrast.matrix)

# View(as.data.frame(colData(se)))
fit2 <- eBayes(fit2) 
```

``` r
res_macrophage_IFNg_vs_naive_dds <- results(dds,
  contrast = c("condition", "IFNg", "naive"),
  lfcThreshold = 1, alpha = 0.05
)
summary(res_macrophage_IFNg_vs_naive_dds)
#> 
#> out of 17806 with nonzero total read count
#> adjusted p-value < 0.05
#> LFC > 1.00 (up)    : 416, 2.3%
#> LFC < -1.00 (down) : 141, 0.79%
#> outliers [1]       : 95, 0.53%
#> low counts [2]     : 0, 0%
#> (mean count < 3)
#> [1] see 'cooksCutoff' argument of ?results
#> [2] see 'independentFiltering' argument of ?results
```

``` r
res_macrophage_IFNg_vs_naive_limma <- topTable(fit2, coef = "IFNg_vs_naive", adjust = "fdr", number = Inf, confint = TRUE)
```

## Package usage

### Converting SummarizedExperiment or DE results to Metabolon format

The requirement for this package to work is that the `CLIENT_SAMPLE_ID`
is the same as the colnames of the summarized experiment.

You can call these two functions to convert your data to a format that
can be uploaded to the Metabolon platform in the multiomics section. It
automatically formats the different expression results to the format
requested by Metabolon.

It works with DESeq2 and limma results, you can also provide a custom
data frame with the required columns. The `cdt` argument is the path to
the client data table file that you want to use to get the sample
metadata.

In addition to this you can also convert the assay data of a summarized
experiment object to the format requested by Metabolon using the
function `se_to_metabolon`.

Important note: for metabolon to be able to run the pathway analysis the
feature annotation must be in the format requested by them. For
transcriptomics data the recommended format is ENSEMBL gene ids, for
proteomics data the recommended format is UNIPROT ids.

``` r
table1 <- se_to_metabolon(gse,
  cdt = cdt_path,
  save_file = F
)

print(head(table1[,c(1:10)]))
#>                PARENT_SAMPLE_NAME ENSG00000000003 ENSG00000000005
#> SAMEA103885102               <NA>          83.818               0
#> SAMEA103885347               <NA>         144.156               0
#> SAMEA103885043               <NA>          71.239               0
#> SAMEA103885392               <NA>         130.372               0
#> SAMEA103885182               <NA>         241.527               0
#> SAMEA103885136               <NA>         140.124               0
#>                ENSG00000000419 ENSG00000000457 ENSG00000000460 ENSG00000000938
#> SAMEA103885102        1462.001         616.732         913.913         210.999
#> SAMEA103885347        1380.000        1058.364         330.360         273.001
#> SAMEA103885043        1117.000         436.315         292.047         635.000
#> SAMEA103885392         961.000         703.578         119.337        1210.000
#> SAMEA103885182        1399.000         580.188         492.361         528.000
#> SAMEA103885136        1276.000        1036.716         173.410         669.000
#>                ENSG00000000971 ENSG00000001036 ENSG00000001084
#> SAMEA103885102         508.000        7715.294        6363.152
#> SAMEA103885347        8355.997        6376.500       15601.710
#> SAMEA103885043         657.999        6231.000        2356.371
#> SAMEA103885392        6407.001        5307.074        8890.247
#> SAMEA103885182         140.999        5099.999        5275.775
#> SAMEA103885136        6670.997        4314.000       14894.069

res_1 <- results_to_metabolon(res_macrophage_IFNg_vs_naive_dds,
  format = "DESeqResults",
  omics = "transcriptomics",
  custom_colnames = NULL,
  save_file = F
)

print(head(table1[,c(1:10)]))
#>                PARENT_SAMPLE_NAME ENSG00000000003 ENSG00000000005
#> SAMEA103885102               <NA>          83.818               0
#> SAMEA103885347               <NA>         144.156               0
#> SAMEA103885043               <NA>          71.239               0
#> SAMEA103885392               <NA>         130.372               0
#> SAMEA103885182               <NA>         241.527               0
#> SAMEA103885136               <NA>         140.124               0
#>                ENSG00000000419 ENSG00000000457 ENSG00000000460 ENSG00000000938
#> SAMEA103885102        1462.001         616.732         913.913         210.999
#> SAMEA103885347        1380.000        1058.364         330.360         273.001
#> SAMEA103885043        1117.000         436.315         292.047         635.000
#> SAMEA103885392         961.000         703.578         119.337        1210.000
#> SAMEA103885182        1399.000         580.188         492.361         528.000
#> SAMEA103885136        1276.000        1036.716         173.410         669.000
#>                ENSG00000000971 ENSG00000001036 ENSG00000001084
#> SAMEA103885102         508.000        7715.294        6363.152
#> SAMEA103885347        8355.997        6376.500       15601.710
#> SAMEA103885043         657.999        6231.000        2356.371
#> SAMEA103885392        6407.001        5307.074        8890.247
#> SAMEA103885182         140.999        5099.999        5275.775
#> SAMEA103885136        6670.997        4314.000       14894.069
```

### Converting DE results to Metabolon format

The package works also with limma results using topTable format:

``` r
res_2 <- results_to_metabolon(res_macrophage_IFNg_vs_naive_limma,
  format = "topTable",
  omics = "transcriptomics",
  custom_colnames = NULL,
  save_file = F
)

print(head(table1[,c(1:10)]))
#>                PARENT_SAMPLE_NAME ENSG00000000003 ENSG00000000005
#> SAMEA103885102               <NA>          83.818               0
#> SAMEA103885347               <NA>         144.156               0
#> SAMEA103885043               <NA>          71.239               0
#> SAMEA103885392               <NA>         130.372               0
#> SAMEA103885182               <NA>         241.527               0
#> SAMEA103885136               <NA>         140.124               0
#>                ENSG00000000419 ENSG00000000457 ENSG00000000460 ENSG00000000938
#> SAMEA103885102        1462.001         616.732         913.913         210.999
#> SAMEA103885347        1380.000        1058.364         330.360         273.001
#> SAMEA103885043        1117.000         436.315         292.047         635.000
#> SAMEA103885392         961.000         703.578         119.337        1210.000
#> SAMEA103885182        1399.000         580.188         492.361         528.000
#> SAMEA103885136        1276.000        1036.716         173.410         669.000
#>                ENSG00000000971 ENSG00000001036 ENSG00000001084
#> SAMEA103885102         508.000        7715.294        6363.152
#> SAMEA103885347        8355.997        6376.500       15601.710
#> SAMEA103885043         657.999        6231.000        2356.371
#> SAMEA103885392        6407.001        5307.074        8890.247
#> SAMEA103885182         140.999        5099.999        5275.775
#> SAMEA103885136        6670.997        4314.000       14894.069
```

### Converting Metabolon client data table to SummarizedExperiment

A central functionality is the ability to convert Metabolon Client data
table files to SummarizedExperiment objects. You can do this using the
function `cdt_to_se`, the output will be a SummarizedExperiment object
that you can use in your analysis pipelines. It automatically merges the
compound information and the subject metdata and adds them respectively
to the rowData and colData of The SummarizedExperiment object.

You can specify which annotation type you want to use for the rownames
of the SummarizedExperiment object using the argument `rowdata_key`.
Note: the compounds that do not have the specified annotation type will
be removed from the final object, as the rownames cannot be duplicated
or empty.

The `data_type` argument allows you to select which type of data you
want to load from the client data table.

You can also save the output SummarizedExperiment object as an RDS file
using the `output_file` argument. If it is set to `FALSE`, the object is
just returned without saving, the output_file argument in this case is
of course ignored.

``` r
se <- cdt_to_se(
  cdt = cdt_path,
  output_file = "test/cdt_to_se_bni.rds",
  save_file = F,
  data_type = "batch_norm_imputed",
  rowdata_key = "INCHIKEY"
)
se
#> class: SummarizedExperiment 
#> dim: 6 24 
#> metadata(0):
#> assays(1): counts
#> rownames(6): FAKEINCHIKEY00001 FAKEINCHIKEY00002 ... FAKEINCHIKEY00005
#>   FAKEINCHIKEY00006
#> rowData names(21): CHEM_ID LIB_ID ... PUBCHEM PLATFORM
#> colnames(24): TEST-1 TEST-2 ... TEST-23 TEST-24
#> colData names(34): PARENT_SAMPLE_NAME CLIENT_IDENTIFIER ... TREATMENT
#>   WT_KO
```

A small trick: you can convert the SummarizedExperiment object to a
SingleCellExperiment object to be able to store additional information
such as dimensionality reductions in it.

``` r
library(SingleCellExperiment)
#> 
#> Attaching package: 'SingleCellExperiment'
#> The following objects are masked from 'package:edgeR':
#> 
#>     cpm, tpm

sce <- as(se, "SingleCellExperiment")

# As you can see here we can add dimensionality reductions (as many as we need) to the object using the capabilities of a single cell summarized experiment
pca_result <- prcomp(t(assay(se)))

reducedDim(sce, "PCA") <- pca_result$x
```

## Session Information

This document has been generated with…

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.32.0 DESeq2_1.50.2              
#>  [3] edgeR_4.8.0                 limma_3.66.0               
#>  [5] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [7] GenomicRanges_1.62.0        Seqinfo_1.0.0              
#>  [9] IRanges_2.44.0              S4Vectors_0.48.0           
#> [11] BiocGenerics_0.56.0         generics_0.1.4             
#> [13] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [15] macrophage_1.26.0           metabolonIO_0.2.0          
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.54           bslib_0.9.0        
#>  [4] ggplot2_4.0.1       lattice_0.22-7      vctrs_0.6.5        
#>  [7] tools_4.5.2         parallel_4.5.2      tibble_3.3.0       
#> [10] pkgconfig_2.0.3     Matrix_1.7-4        RColorBrewer_1.1-3 
#> [13] S7_0.2.1            desc_1.4.3          readxl_1.4.5       
#> [16] lifecycle_1.0.4     compiler_4.5.2      farver_2.1.2       
#> [19] textshaping_1.0.4   statmod_1.5.1       codetools_0.2-20   
#> [22] htmltools_0.5.8.1   sass_0.4.10         yaml_2.3.11        
#> [25] pkgdown_2.2.0       pillar_1.11.1       jquerylib_0.1.4    
#> [28] BiocParallel_1.44.0 DelayedArray_0.36.0 cachem_1.1.0       
#> [31] abind_1.4-8         locfit_1.5-9.12     digest_0.6.39      
#> [34] fastmap_1.2.0       grid_4.5.2          cli_3.6.5          
#> [37] SparseArray_1.10.3  magrittr_2.0.4      S4Arrays_1.10.0    
#> [40] scales_1.4.0        rmarkdown_2.30      XVector_0.50.0     
#> [43] cellranger_1.1.0    ragg_1.5.0          evaluate_1.0.5     
#> [46] knitr_1.50          rlang_1.1.6         Rcpp_1.1.0         
#> [49] glue_1.8.0          jsonlite_2.0.0      R6_2.6.1           
#> [52] systemfonts_1.3.1   fs_1.6.6
```
