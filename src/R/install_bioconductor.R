if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
BiocManager::install("BSgenome")
BiocManager::install("CRISPRseek")

