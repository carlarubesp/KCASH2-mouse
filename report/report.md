# KCASH2 Mouse RNA-seq Analysis — Report

## **Status:** work in progress (added to the readme &)
- [ ] Project overview
- [ ] Study design
- [ ] Aims, which questions are we making about this project.
- [x] Environment setup
- [ ] Data Preprocessing and Quality Assessment
- [ ] Differential expression analysis
- [ ] GO and pathway analysis

## 1. Environment Setup

### 1.1 Package manager

All packages are installed using **BiocManager**, the official package manager for Bioconductor, one of the two main R repositories for bioinformatics, specifically designed for genomic data analysis.

Using BiocManager instead of base `install.packages()` ensures that all packages are mutually compatible with the current versions of R and Bioconductor, avoiding version conflicts between CRAN and Bioconductor dependencies. For packages hosted on CRAN, BiocManager calls `install.packages()` internally, so there is no difference in the result — only in the version management.

This analysis uses **R version 4.6.0** and **Bioconductor version 3.23**.

### 1.2 Packages

The packages installed for the analysis are:

| Package | Repository | Purpose |
|---|---|---|
| **jsonlite** | CRAN | Parsing and handling JSON files; used for reading configuration and metadata files |
| **readxl** | CRAN | Reading Excel files (.xlsx) into R |
| **DESeq2** | Bioconductor | Main tool for differential gene expression analysis |
| **edgeR** | Bioconductor | Used here for gene filtering (`filterByExpr`); also used to validate DESeq2 results |
| **clusterProfiler** | Bioconductor | Gene Ontology (GO) and KEGG pathway enrichment analysis |
| **org.Mm.eg.db** | Bioconductor | Mouse gene annotation database; maps gene symbols to Entrez IDs |
| **AnnotationDbi** | Bioconductor | Interface for querying annotation databases such as org.Mm.eg.db |
| **EnhancedVolcano** | Bioconductor | Volcano plots for differential expression results |
| **pheatmap** | CRAN | Clustered heatmaps for visualising expression patterns |
| **ggplot2** | CRAN | General-purpose data visualisation |
| **dplyr** | CRAN | Data manipulation and filtering |
| **tidyr** | CRAN | Reshaping data frames for analysis and plotting |
| **RColorBrewer** | CRAN | Colour palettes for figures |
| **sva** | Bioconductor | Batch effect detection and correction |

&nbsp;

## 2. Data Preprocessing and Quality Assessment
