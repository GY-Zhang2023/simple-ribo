.libPaths(Sys.getenv("R_LIBS_USER"))

env_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (!dir.exists(env_lib)) dir.create(env_lib, recursive = TRUE)
.libPaths(env_lib)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "BiocParallel", "S4Vectors", "IRanges", "GenomicRanges",
  "GenomicFeatures", "GenomicAlignments", "Biostrings",
  "txdbmaker", "DESeq2", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "GO.db"
), update = FALSE)
install.packages("~/simple-ribo/riboWaltz-master",
                 repos = NULL,
                 type = "source",
                 lib = env_lib,
                 dependencies = FALSE)
