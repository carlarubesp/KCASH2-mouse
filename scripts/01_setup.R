# ==============================================================================
# SCRIPT 01: Environment setup and package installation
# ==============================================================================

# Install BiocManager if it's not already installed
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# List of required packages based on the analysis scripts
packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "EnhancedVolcano",
  "pheatmap",
  "ggplot2",
  "RColorBrewer"
)

# Install packages using BiocManager
BiocManager::install(packages, update = FALSE)

print("All required packages are successfully installed.")