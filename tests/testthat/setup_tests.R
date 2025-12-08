# ---- Setup ----
# Load required libraries
library(metabolonIO)
suppressPackageStartupMessages({
  library(macrophage)
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(DESeq2)
})

# ---- Load Dataset ----
cdt_path <- system.file("extdata", "example_data_tables.XLSX", package = "metabolonIO")

data(gse)
# Trim rownames to first 15 characters
rownames(gse) <- substr(rownames(gse), 1, 15)

# Note: CLIENT_SAMPLE_ID must match colnames of the SummarizedExperiment for metabolonIO functions

# ---- DESeq2 Analysis ----
dds <- DESeqDataSet(gse, design = ~condition)

# Filter genes with low counts
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Differential expression for IFNg vs naive
res_macrophage_IFNg_vs_naive_dds <- results(
  dds,
  contrast = c("condition", "IFNg", "naive"),
  lfcThreshold = 1,
  alpha = 0.05
)
summary(res_macrophage_IFNg_vs_naive_dds)

# ---- limma Analysis ----
condition <- factor(colData(gse)[, "condition_name"])
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

contrast.matrix <- makeContrasts(
  IFNg_vs_naive = IFNg - naive,
  levels = design
)

# Prepare data for limma
d0 <- DGEList(assay(gse))
d0 <- calcNormFactors(d0)

# Filter lowly expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]

# Apply voom transformation
matrix_to_fit_model_to <- voom(d, design, plot = FALSE)

# Fit linear model and apply contrasts
fit <- lmFit(matrix_to_fit_model_to, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top DE genes
res_macrophage_IFNg_vs_naive_limma <- topTable(
  fit2,
  coef = "IFNg_vs_naive",
  adjust = "fdr",
  number = Inf,
  confint = TRUE
)