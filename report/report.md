# KCASH2 Mouse RNA-seq Analysis — Report

## **Status:** work in progress (added to the readme &)
- [x] Project overview
- [x] Study design
- [ ] Aims, which questions are we making about this project.
- [x] Environment setup
- [ ] Data Preprocessing and Quality Assessment
- [ ] Differential expression analysis
- [ ] GO and pathway analysis

## Project Overview

### Background

The *KCASH2* gene acts as a tumor suppressor and plays a critical role in colorectal cancer, primarily by negatively regulating the Hedgehog (Hh) signaling pathway. While its biochemical function is documented, understanding how the loss of this gene globally affects the cellular transcriptome is essential to uncover its underlying mechanisms in tumor development and progression.

### Objective
This project aims to investigate the transcriptomic impact of *KCASH2* using an in vivo mouse model. Specifically, we analyze RNA-sequencing (RNA-seq) data to identify Differentially Expressed Genes (DEGs) and altered biological pathways. The study evaluates expression profiles across two main dimensions: tissue condition (Tumor vs. Normal) and genetic background (Wild Type [WT] vs. *KCASH2* Knockout [KO]).

&nbsp;

## Study Design

The experimental design of this study is structured to evaluate the transcriptomic changes driven by both tumor development and the genetic deletion of *KCASH2*.

### Experimental model and sample collection

The study utilizes a genetically modified mouse model, divided into two distinct genotype groups:
* **WT (Wild-Type)**: Mice with normal KCASH2 expression.
* **KO (Knockout)**: Mice genetically engineered to lack KCASH2 expression.

From these animals, two types of tissue samples were harvested for bulk RNA-sequencing: **Tumor tissue (T)** and adjacent **Normal tissue (N)**.

### Paired strategy and sample size

A major strength of this study design is its paired nature. Whenever possible, both tumor and normal tissues were extracted from the exact same animal. This allows us to control for inter-individual genetic background noise (baseline variability between different mice), significantly increasing the statistical power of the analysis.

After initial quality control and the removal of technical outliers, the refined dataset used for the core paired analyses consists of **11 complete animal pairs** (6 complete KO pairs and 5 complete WT pairs), yielding a robust dataset for differential expression modeling.

### Statistical contrasts

Based on this sample structure, the study is designed to explore two main biological axes:
1. The Genotype Effect (Unpaired Analysis): Evaluating how the absence of KCASH2 alters gene expression within a specific tissue environment.
    * Tumor KO vs. Tumor WT (T_KO vs. T_WT)
    * Normal KO vs. Normal WT (N_KO vs. N_WT)

2. The Tumorigenesis Effect (Paired Analysis): Evaluating the transcriptomic shift from healthy to tumor tissue, utilizing the animal ID as a blocking factor to isolate the disease effect.
    * Tumor KO vs. Normal KO (T_KO vs. N_KO)
    * Tumor WT vs. Normal WT (T_WT vs. N_WT)

&nbsp;

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
| **DESeq2** | Bioconductor | Main tool for differential gene expressionproj analysis |
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
