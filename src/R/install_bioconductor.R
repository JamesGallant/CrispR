if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures")
if (!require("BSgenome")) BiocManager::install("BSgenome")
if (!require("CRISPRseek")) BiocManager::install("CRISPRseek")
