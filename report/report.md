# KCASH2 Mouse RNA-seq Analysis — Report

## Project Overview

### Background

The *KCASH2* gene acts as a tumor suppressor and plays a critical role in colorectal cancer, primarily by negatively regulating the Hedgehog (Hh) signaling pathway. While its biochemical function is documented, understanding how the loss of this gene globally affects the cellular transcriptome is essential to uncover its underlying mechanisms in tumor development and progression.

### Objective
This project aims to investigate the transcriptomic impact of *KCASH2* using an in vivo mouse model. Specifically, we analyze RNA-sequencing (RNA-seq) data to identify Differentially Expressed Genes (DEGs) and altered biological pathways. The study evaluates expression profiles across two main dimensions: tissue condition (Tumor vs. Normal) and genetic background (Wild Type [WT] vs. *KCASH2* Knockout [KO]).

&nbsp;

## Research Aims

This project was designed to answer the following research questions:

1. Does the loss of `KCASH2` produce detectable transcriptomic differences between KO and WT mice within the same tissue type?
2. Which genes are differentially expressed between tumor and adjacent normal tissue within each genotype?
3. Does the paired tumor-versus-normal design reveal stronger transcriptional changes than the direct KO-versus-WT comparisons?
4. Which biological processes and molecular pathways are enriched among the differentially expressed genes?
5. Are Hedgehog-related genes or Hedgehog-associated pathways altered in tumor tissue?
6. Does KCASH2 loss produce a detectable pathway-level activation of Hedgehog signaling at the transcriptomic level?

These questions separate two main biological effects: the genotype effect caused by KCASH2 deletion and the tumorigenesis effect observed when tumor tissue is compared with adjacent normal tissue. The final aim is to determine whether the loss of KCASH2 is associated with broad transcriptomic changes and whether these changes support alterations in Hedgehog signaling or other cancer-related pathways.

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

## Methods

### Dataset

The input data for this analysis consists of two Excel files:
* `raw_counts_mouse.xlsx`: a gene-level read count matrix where rows represent genes and columns represent individual samples. Each cell contains the number of sequencing reads mapped to a given gene in a given sample.
* `metadata_mouse.xlsx`: a sample annotation table describing the experimental group (type: T_KO, T_WT, N_KO, or N_WT), the tissue type (Tumor or Normal), the genotype (KO or WT), and the animal ID used to link paired samples from the same individual.

Samples are identified by IonCode barcodes (e.g., IonCode_0103), consistent with Ion Torrent multiplexed sequencing. Data were provided as pre-computed count tables and no upstream processing (e.g., read alignment or quantification) was performed within this pipeline.

After quality control and outlier removal, the final dataset used for the paired analyses contained 11 complete tumor-normal animal pairs: 6 KO pairs and 5 WT pairs.

### 1. Environment Setup

#### 1.1 Package manager

All packages are installed using **BiocManager**, the official package manager for Bioconductor, one of the two main R repositories for bioinformatics, specifically designed for genomic data analysis.

Using BiocManager instead of base `install.packages()` ensures that all packages are mutually compatible with the current versions of R and Bioconductor, avoiding version conflicts between CRAN and Bioconductor dependencies. For packages hosted on CRAN, BiocManager calls `install.packages()` internally, so there is no difference in the result, only in the version management.

This analysis uses **R version 4.6.0** and **Bioconductor version 3.23**.

#### 1.2 Packages

The packages installed for the analysis are:

| Package | Repository | Purpose |
|----------|------------|---------|
| **readxl** | CRAN | Reading Excel files (.xlsx) into R |
| **readr** | CRAN | Reading and writing CSV files |
| **dplyr** | CRAN | Data manipulation and filtering |
| **tidyr** | CRAN | Reshaping data frames for analysis and plotting |
| **stringr** | CRAN | String manipulation, used for formatting text in enrichment plots |
| **DESeq2** | Bioconductor | Main tool for differential gene expression analysis |
| **edgeR** | Bioconductor | Used strictly for filtering low-expression genes (`filterByExpr`) prior to modeling |
| **clusterProfiler** | Bioconductor | Gene Ontology (GO) and KEGG pathway enrichment analysis |
| **org.Mm.eg.db** | Bioconductor | Mouse gene annotation database; maps gene symbols to Entrez IDs |
| **AnnotationDbi** | Bioconductor | Interface for querying annotation databases such as `org.Mm.eg.db` |
| **enrichplot** | Bioconductor | Visualization methods for `clusterProfiler` enrichment results |
| **EnhancedVolcano** | Bioconductor | Volcano plots for differential expression results |
| **pheatmap** | CRAN | Clustered heatmaps for visualising expression patterns |
| **ggplot2** | CRAN | General-purpose data visualisation |
| **RColorBrewer** | CRAN | Colour palettes for figures |

&nbsp;

#### 1.3 Analysis settings and cut-offs

Low-expression genes were filtered before differential expression analysis using `edgeR::filterByExpr()`.

Differential expression analysis was performed with DESeq2. Genes were considered differentially expressed when they satisfied both `padj < 0.05` and `|log2FC| ≥ 1`.

For unpaired genotype comparisons, the model used `type` as the explanatory variable. For paired tumor-versus-normal comparisons, animal ID was included as a blocking factor to control for baseline differences between mice.

GO Biological Process enrichment was performed with `clusterProfiler::enrichGO()` using `ont = "BP"`, Benjamini-Hochberg correction, `pvalueCutoff = 0.05`, and `qvalueCutoff = 0.2`.

KEGG enrichment was performed with `clusterProfiler::enrichKEGG()` using the mouse organism code `mmu`, Benjamini-Hochberg correction, and `pvalueCutoff = 0.05`.

Comparisons with fewer than 10 mapped Entrez IDs were not interpreted for enrichment analysis.

### 2. Data Preprocessing and Quality Assessment

#### Data loading and alignment

Raw count data and sample metadata were imported from Excel using the `readxl package`. Metadata columns for tissue type, genotype, and animal ID were derived from the `type` column and stored as factors. The sample order in the count matrix was verified and, where necessary, the metadata was reordered to ensure row-by-column alignment between the two tables.

#### Low-expression gene filtering

Genes with very low or near-zero counts across all samples were removed prior to statistical modeling to reduce noise and improve multiple-testing correction power. Filtering was applied using the `filterByExpr()` function from the `edgeR` package, with genotype as the grouping variable. This function retains genes that have at least a minimum count threshold in a number of samples no smaller than the smallest experimental group.

#### Principal Component Analysis (PCA) and outlier identification

Sample-level quality was assessed using PCA. Count data were first transformed using the Variance Stabilizing Transformation (VST) implemented in DESeq2 (`vst()`, `blind = TRUE`), which stabilizes the variance across the range of counts independently of any experimental design. PCA was then computed on the transformed values using `plotPCA()`, and samples were plotted on the first two principal components, colored by tissue type and shaped by genotype.

Two samples (`IonCode_0103` and `IonCode_0105`) were identified as technical outliers. These samples had substantially lower sequencing depth (approximately 3–4 million reads) compared to the rest of the dataset (approximately 12–18 million reads per sample), and they clustered away from other samples of their respective experimental group on the PCA plot. Both samples were excluded from all downstream analyses.

#### Post-removal quality check

PCA was repeated on the cleaned dataset to confirm that, after outlier removal, samples clustered consistently by tissue type along the first principal component, as expected for a bulk RNA-seq dataset where the tissue of origin is the primary source of variation. The cleaned count matrix and metadata were saved as CSV files in `data/clean/` for use in all subsequent steps.

### 3. Differential expression analysis

#### Methods & DEGs Summary

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


## 4. Differential Expression Visualization

To visually evaluate the statistical outputs generated by the DESeq2 pipeline, high-throughput visualizations were implemented. These plots capture both the global landscape of transcriptional alterations and targeted expression patterns of biological pathway components across the experimental contrasts.

### 4.1 Differential Expression Landscape (Volcano Plots)
* **Comparison 1 (T_KO vs T_WT):** Stood as a clean baseline with 0 significant DEGs under the strict threshold ($padj < 0.05, |log2FC| \ge 1$), showing that removing the suppressor within already established tumor environments does not induce further global massive transcriptional shifts.
* **Comparison 2 (N_KO vs N_WT):** Identified strictly 2 significant DEGs (*Slc7a5* and *2340079G19Rik*), indicating that knocking out **KCASH2 (*Kctd21*)** in normal tissue acts as a quiet metabolic priming stage rather than an immediate, full-scale tumorigenic driver on its own.
* **Comparison 3 (T_KO vs N_KO):** Revealed a massive molecular explosion of 2,452 significant DEGs. High-intensity oncogenic and tissue-remodeling markers were automatically highlighted at the peak of statistical significance, prominently featuring genes such as *Tgm3*, *Adh1*, *Atp12a*, *Bmp3*, and *Hes6*, showcasing a profound systemic shift between the tumor and normal states in the knockout background.
* **Comparison 4 (T_WT vs N_WT):** Captured the classical colorectal cancer baseline with 3,223 significant DEGs. This landscape is visually marked by the extreme suppression of normal tissue homeostasis markers like *Aqp8* (water channels) and *Zg16* (protective mucosal barrier), against the sharp activation of tumor-associated markers like *Cxcl5* and *Stra6*.

### 4.2 High-Throughput Expression Profiles (Top 50 Clustered Heatmap)
The bidirectional hierarchical clustering of the top 50 most significant DEGs from Comparison 3 (Tumor KO vs Normal KO) successfully partitioned the samples into discrete, high-contrast expression blocks. Tumoral samples showed an intense, uniform shift in gene expression profiles compared to the normal backgrounds, confirming that the statistical sorting accurately captured the definitive transcriptional signature of the disease tissue.

### 4.3 Targeted Pathway Analysis (Hedgehog & Gli Signaling Heatmap)
To address the core project hypothesis, a dedicated, low-dimensional heatmap successfully isolated the key components of the Hedgehog pathway (*Gli1*, *Gli2*, *Gli3*, *Ptch1*, *Smo*, *Shh*, and **KCASH2/*Kctd21***). The visual profile demonstrated a subtle, coordinated regulatory fine-tuning across the cohorts. Crucially, the expression levels of the suppressor **KCASH2** visually partitioned across the sub-clusters alongside changes in the transcription factors *Gli1* and *Gli3*, confirming the visual presence of pathway modulation within specific sample groups.


#### 5. Functional Enrichment Analysis and Biological Interpretation

##### Input DEG lists

Functional enrichment analysis was performed using the DEG lists obtained from the official differential expression analysis. DEGs were defined using the same thresholds applied in the DESeq2 step: `padj < 0.05` and `|log2FC| ≥ 1`.

The direct genotype comparisons produced no or very few significant DEGs:

* `T_KO vs T_WT`: 0 DEGs
* `N_KO vs N_WT`: 2 DEGs

For this reason, the functional analysis focused mainly on the paired tumor-versus-normal comparisons:

* `T_KO vs N_KO`: 2452 DEGs
* `T_WT vs N_WT`: 3223 DEGs

This indicates that the dominant transcriptional signal in the dataset is associated with tumor development rather than with a broad genotype effect of KCASH2 loss.

##### Gene ID conversion and enrichment strategy

Mouse gene symbols were converted into Entrez IDs using the `org.Mm.eg.db` annotation database. These converted gene lists were then used for Gene Ontology Biological Process enrichment and KEGG pathway analysis with `clusterProfiler`.

GO enrichment was used to identify overrepresented biological processes, while KEGG enrichment was used to detect altered molecular pathways. In addition, a targeted screen was performed for pathways directly related to the project hypothesis, especially Hedgehog signaling, Wnt/β-catenin signaling, and inflammation-related genes.

##### GO Biological Process enrichment

In the `T_KO vs N_KO` comparison, the most enriched GO Biological Process terms were mainly related to:

* chemotaxis;
* taxis;
* cell chemotaxis;
* leukocyte migration;
* organic anion transport;
* ameboidal-type cell migration.

These results suggest that tumor development in the KO background is associated with immune-cell recruitment, cell migration, inflammatory signalling, and transport-related processes.

The `T_WT vs N_WT` comparison showed a similar functional pattern. The top enriched GO terms included:

* chemotaxis;
* taxis;
* ameboidal-type cell migration;
* regulation of epithelial cell proliferation;
* tissue migration;
* epithelial cell migration.

Overall, the GO results show that both paired tumor-versus-normal comparisons are enriched for processes related to migration, immune response, tissue remodelling, and epithelial behaviour.

##### KEGG pathway enrichment

KEGG enrichment analysis also supported the presence of cancer-relevant and microenvironment-related alterations.

In the `T_KO vs N_KO` comparison, enriched pathways included:

* cytokine-cytokine receptor interaction;
* viral protein interaction with cytokine and cytokine receptor;
* pancreatic secretion;
* ECM-receptor interaction;
* mineral absorption;
* rheumatoid arthritis.

These pathway names should not be interpreted literally as different diseases. Instead, they reflect pathway modules involving cytokines, extracellular matrix components, immune response, and tissue remodelling.

In the `T_WT vs N_WT` comparison, KEGG enrichment highlighted pathways such as:

* PI3K-Akt signalling;
* ECM-receptor interaction;
* cytokine-cytokine receptor interaction;
* integrin signalling;
* Ras signalling;
* cancer-related KEGG terms.

These pathways are relevant in a tumor context because they are associated with cell survival, proliferation, migration, extracellular matrix remodelling, cytokine signalling, and interaction with the tumor microenvironment.

##### Hedgehog and Wnt targeted analysis

Because the project focuses on the possible tumor suppressor role of KCASH2 and its relationship with Hedgehog signaling, Hedgehog-related genes were examined separately.

In the `T_KO vs N_KO` comparison, `Shh`, `Ptch2`, and `Sufu` were found among the DEGs and were upregulated in tumor tissue. In the `T_WT vs N_WT` comparison, `Shh` and `Ptch2` were also upregulated.

This shows that individual Hedgehog-associated genes change during the tumor-versus-normal transition in both genetic backgrounds. However, the targeted pathway screen did not identify Hedgehog signalling as significantly enriched after multiple-testing correction. Therefore, the results support changes in individual Hedgehog-related genes, but they do not provide strong evidence for global pathway-level enrichment of Hedgehog signalling.

The targeted pathway screen also detected several Wnt/β-catenin-related genes among the DEGs, especially in the tumor-versus-normal comparisons. These included genes such as `Axin2`, `Lef1`, `Wnt7b`, `Wnt5b`, `Fzd10`, `Tcf7`, `Dkk2`, `Dkk3`, `Wnt3`, `Wnt5a`, `Fzd8`, and `Fzd3`. However, after multiple-testing correction, Wnt/β-catenin was not significantly enriched either. Therefore, Wnt-related genes are present among the DEGs, but the pathway-level result should be interpreted cautiously.

##### Functional interpretation

Overall, the functional analysis suggests that the dominant biological signal in this dataset is the tumor-versus-normal transition. The direct KO-versus-WT comparisons showed no or very few DEGs, while both paired tumor-versus-normal comparisons produced large DEG sets and strong functional enrichment.

The enriched GO and KEGG results point mainly to:

* cell migration;
* chemotaxis;
* immune and cytokine signalling;
* epithelial proliferation;
* extracellular matrix interaction;
* integrin signalling;
* PI3K-Akt signalling;
* Ras-related pathways.

Taken together, these results do not support a strong global transcriptional effect of KCASH2 loss when comparing KO and WT samples directly. Instead, they show that tumor development is associated with broad transcriptional and functional changes in both genotypes. Hedgehog-related genes such as `Shh`, `Ptch2`, and `Sufu` were altered in tumor-versus-normal comparisons, but Hedgehog signalling was not significantly enriched as a complete pathway. Therefore, the role of Hedgehog in this dataset should be interpreted carefully.

### 6. Visualization

#### Quality control plots

PCA plots were generated during quality control to evaluate the overall structure of the RNA-seq samples. The initial PCA was used to identify technical outliers, while the PCA after outlier removal was used to confirm that the cleaned dataset showed a more consistent sample distribution.

This step was important because downstream differential expression analysis depends strongly on the quality and comparability of the samples.

#### Differential expression plots

Volcano plots were generated for each differential expression comparison. These plots summarize both the magnitude of expression change and the statistical significance of each gene.

The direct genotype comparisons, `T_KO vs T_WT` and `N_KO vs N_WT`, showed few or no significant genes, which is consistent with the DEG summary table. In contrast, the paired tumor-versus-normal comparisons showed a much stronger differential expression signal, especially `T_WT vs N_WT`.

Heatmaps were also generated for the most significant genes in the main tumor-versus-normal comparisons. These plots were used to visualize whether the selected DEGs showed coherent expression patterns across samples. The clustering structure helped confirm that the strongest separation was related to tumor versus normal tissue rather than to genotype alone.

#### Functional enrichment plots

GO and KEGG dotplots were generated from the enrichment results.

GO dotplots were used to display the most enriched Biological Process terms, including processes related to chemotaxis, leukocyte migration, epithelial migration, tissue migration, and organic anion transport.

KEGG dotplots were used to summarize pathway-level enrichment, including PI3K-Akt signalling, ECM-receptor interaction, cytokine-cytokine receptor interaction, integrin signalling, Ras signalling, and cancer-related pathway modules.

A `compareCluster` GO plot was also generated to compare enriched biological processes across DEG contrasts. This visualization helped show that the paired tumor-versus-normal comparisons carried the strongest functional signal.

#### Visualization summary

Overall, the visualization results support the main conclusion of the analysis: the dominant signal in the dataset is the tumor-versus-normal transition, while the direct KO-versus-WT comparisons show limited transcriptomic differences.

The figures also support the functional interpretation that tumor tissue is associated with changes in migration, immune and cytokine signalling, extracellular matrix interaction, epithelial behaviour, and cancer-related signalling pathways.
