# Environment setup and package installation

# Install BiocManager
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Packages
packages <- c(
  "jsonlite",
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "EnhancedVolcano",
  "pheatmap",
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "sva"
)

BiocManager::install(packages, update = FALSE)
print("Setup successfully completed! All required packages are installed.")