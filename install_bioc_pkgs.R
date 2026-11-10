if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "BiocParallel", "S4Vectors", "IRanges", "GenomicRanges",
  "GenomicFeatures", "GenomicAlignments", "Biostrings",
  "txdbmaker", "DESeq2", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "GO.db"
), update = FALSE)
