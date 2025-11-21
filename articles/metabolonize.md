# metabolonize

## setup

``` r
library("metabolonize")
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

``` r
cdt_path <- system.file("extdata", "example_data_tables.XLSX", package = "metabolonize")

data(gse)
rownames(gse) <- substr(rownames(gse), 1, 15) # per togliere tutto quello dopo il punto
```

The requirement for this package to work is that the CLIENT_SAMPLE_ID is
the same as the colnames of the summarized experiment.

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

You can call these two functions to convert your data to a format that
can be uploaded to the Metabolon platform in the multiomics section. It
automatically formats the results from DESeq2 or limma to the required
format.

In addition to this you can also convert the assay data of a summarized
experiment object to the format requested by Metabolon using the
function `se_to_metabolon`.

``` r
table1 <- se_to_metabolon(gse,
  cdt = cdt_path,
  input_features = "ensembl_id",
  organism = "Hs",
  save_file = F,
  sample_id_column = "CLIENT_SAMPLE_ID"
)
#> 

print(head(table1[,c(1:10)]))
#>                ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457
#> SAMEA103885102          83.818               0        1462.001         616.732
#> SAMEA103885347         144.156               0        1380.000        1058.364
#> SAMEA103885043          71.239               0        1117.000         436.315
#> SAMEA103885392         130.372               0         961.000         703.578
#> SAMEA103885182         241.527               0        1399.000         580.188
#> SAMEA103885136         140.124               0        1276.000        1036.716
#>                ENSG00000000460 ENSG00000000938 ENSG00000000971 ENSG00000001036
#> SAMEA103885102         913.913         210.999         508.000        7715.294
#> SAMEA103885347         330.360         273.001        8355.997        6376.500
#> SAMEA103885043         292.047         635.000         657.999        6231.000
#> SAMEA103885392         119.337        1210.000        6407.001        5307.074
#> SAMEA103885182         492.361         528.000         140.999        5099.999
#> SAMEA103885136         173.410         669.000        6670.997        4314.000
#>                ENSG00000001084 ENSG00000001167
#> SAMEA103885102        6363.152            1605
#> SAMEA103885347       15601.710            1614
#> SAMEA103885043        2356.371            1603
#> SAMEA103885392        8890.247            1628
#> SAMEA103885182        5275.775            1765
#> SAMEA103885136       14894.069            1400

res_1 <- results_to_metabolon(res_macrophage_IFNg_vs_naive_dds,
  format = "DESeqResults",
  omics = "transcriptomics",
  custom_colnames = NULL,
  input_features = "ensembl_id",
  organism = "Hs",
  save_file = F
)

print(head(table1[,c(1:10)]))
#>                ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457
#> SAMEA103885102          83.818               0        1462.001         616.732
#> SAMEA103885347         144.156               0        1380.000        1058.364
#> SAMEA103885043          71.239               0        1117.000         436.315
#> SAMEA103885392         130.372               0         961.000         703.578
#> SAMEA103885182         241.527               0        1399.000         580.188
#> SAMEA103885136         140.124               0        1276.000        1036.716
#>                ENSG00000000460 ENSG00000000938 ENSG00000000971 ENSG00000001036
#> SAMEA103885102         913.913         210.999         508.000        7715.294
#> SAMEA103885347         330.360         273.001        8355.997        6376.500
#> SAMEA103885043         292.047         635.000         657.999        6231.000
#> SAMEA103885392         119.337        1210.000        6407.001        5307.074
#> SAMEA103885182         492.361         528.000         140.999        5099.999
#> SAMEA103885136         173.410         669.000        6670.997        4314.000
#>                ENSG00000001084 ENSG00000001167
#> SAMEA103885102        6363.152            1605
#> SAMEA103885347       15601.710            1614
#> SAMEA103885043        2356.371            1603
#> SAMEA103885392        8890.247            1628
#> SAMEA103885182        5275.775            1765
#> SAMEA103885136       14894.069            1400
```

The package works also with limma results using topTable format:

``` r
res_2 <- results_to_metabolon(res_macrophage_IFNg_vs_naive_limma,
  format = "topTable",
  omics = "transcriptomics",
  custom_colnames = NULL,
  input_features = "ensembl_id",
  organism = "Hs",
  save_file = F
)

print(head(table1[,c(1:10)]))
#>                ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457
#> SAMEA103885102          83.818               0        1462.001         616.732
#> SAMEA103885347         144.156               0        1380.000        1058.364
#> SAMEA103885043          71.239               0        1117.000         436.315
#> SAMEA103885392         130.372               0         961.000         703.578
#> SAMEA103885182         241.527               0        1399.000         580.188
#> SAMEA103885136         140.124               0        1276.000        1036.716
#>                ENSG00000000460 ENSG00000000938 ENSG00000000971 ENSG00000001036
#> SAMEA103885102         913.913         210.999         508.000        7715.294
#> SAMEA103885347         330.360         273.001        8355.997        6376.500
#> SAMEA103885043         292.047         635.000         657.999        6231.000
#> SAMEA103885392         119.337        1210.000        6407.001        5307.074
#> SAMEA103885182         492.361         528.000         140.999        5099.999
#> SAMEA103885136         173.410         669.000        6670.997        4314.000
#>                ENSG00000001084 ENSG00000001167
#> SAMEA103885102        6363.152            1605
#> SAMEA103885347       15601.710            1614
#> SAMEA103885043        2356.371            1603
#> SAMEA103885392        8890.247            1628
#> SAMEA103885182        5275.775            1765
#> SAMEA103885136       14894.069            1400
```

## Converting Metabolon data to SummarizedExperiment

You can also convert a Metabolon Client data table file to a
SummarizedExperiment object using the function `cdt_to_se`, to be able
to easily import the data in R and continue with your analysis.

``` r
se <- cdt_to_se(
  cdt = cdt_path,
  output_file = "test/cdt_to_se_bni.rds",
  save_file = F,
  data_type = "batch_norm_imputed"
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
