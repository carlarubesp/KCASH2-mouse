# KCASH2 Mouse RNA-seq Analysis Pipeline

This repository contains the bioinformatics pipeline to analyze RNA sequencing (RNA-seq) data in a mouse model (*in vivo*), evaluating the role of the tumor suppressor **KCASH2** in colorectal cancer (CRC) and its impact on the Hedgehog (Hh) signaling pathway.

---

## Project Structure

For the scripts to run correctly, please ensure the following folder structure is maintained:

[TODO: folder structure]

&nbsp;

## Environment Requirements

The analysis was developed and validated under the following specifications:
- **R version**: 4.6.0 (or compatible)
- **Bioconductor version**: 3.23
- **Package manager**: BiocManager to ensure mutual compatibility between CRAN and Bioconductor libraries.

&nbsp;

## Execution Instructions
The scripts are designed to run sequentially. You must execute them in the following order (1 to 3):

### Step 1: Environment and library setup

Open an R terminal or your IDE (such as RStudio) and run the first script to install all necessary dependencies:

```R
source("scripts/01_setup.R")
```

*This script will automatically install key packages like DESeq2, edgeR, clusterProfiler, and graphical tools like pheatmap and EnhancedVolcano.*

### Step 2: Preprocessing and quality control

Once the libraries are installed, run the data cleaning script:

```R
source("scripts/02_preprocessing.R")
```

### Step 3: Differential expression analysis (DESeq2)

Next, run the main statistical analysis to find differentially expressed genes (DEGs):

```R
source("scripts/03_DESeq2.R")
```

### Step 4: Visualization

Finally, generate graphical representations of your results by running:

```R
source("scripts/04_visualization.R")
```

&nbsp;

## Generated Results
After completing the execution, you will find the following files ready for interpretation:

[TODO: finish this part]