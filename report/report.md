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

### Methods & DEGs Summary

Differential expression analysis was performed using the cleaned RNA-seq count matrix and the corresponding sample metadata. The count matrix contained gene-level read counts, while the metadata described the sample type, tissue condition, genotype, and animal ID for each sample. The cleaned dataset was generated during preprocessing and used as input for all differential expression comparisons.

The analysis was performed with DESeq2. For each comparison, the count matrix was subset to the relevant sample groups based on the metadata column `type`. The comparisons T_KO vs T_WT and N_KO vs N_WT were analyzed using an unpaired design with `type` as the explanatory variable. These comparisons evaluate the effect of KCASH2 knockout within tumor tissue and normal tissue, respectively.

The comparisons T_KO vs N_KO and T_WT vs N_WT were analyzed using a paired design. For these comparisons, only animals with both tumor and normal samples were retained. Animal ID was included as a blocking factor in the DESeq2 design to account for baseline expression differences between individual mice. This resulted in 6 complete KO pairs and 5 complete WT pairs.

Differentially expressed genes were extracted using the thresholds padj < 0.05 and |log2FC| ≥ 1. Genes with a positive log2 fold change were counted as upregulated in the first group of the comparison, while genes with a negative log2 fold change were counted as downregulated in the first group.

| Comparison   | Design                       | Upregulated DEGs | Downregulated DEGs | Total DEGs |
| ------------ | ---------------------------- | ---------------: | -----------------: | ---------: |
| T_KO vs T_WT | unpaired                     |                0 |                  0 |          0 |
| N_KO vs N_WT | unpaired                     |                1 |                  1 |          2 |
| T_KO vs N_KO | paired, blocked by animal ID |             1306 |               1146 |       2452 |
| T_WT vs N_WT | paired, blocked by animal ID |             1726 |               1497 |       3223 |

The KO vs WT comparisons showed few or no differentially expressed genes, whereas the paired tumor vs normal comparisons showed substantially larger DEG sets. Overall, the strongest transcriptomic differences were observed between tumor and normal tissue within the same genotype.
