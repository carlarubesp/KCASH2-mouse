# ==============================================================================
# SCRIPT 02: Data preprocessing and quality assessment
# ==============================================================================

library(readxl)
library(DESeq2)
library(edgeR)
library(ggplot2)

# Load data
counts_raw <- read_excel("data/raw_counts_mouse.xlsx")
metadata_raw <- read_excel("data/metadata_mouse.xlsx")

counts <- as.data.frame(counts_raw)
metadata <- as.data.frame(metadata_raw)

rownames(counts) <- counts[, 1]
counts <- counts[, -1]

# Inspect counts and metadata
cat("Counts matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
print(head(counts[, 1:5]))  # first 5 samples to keep it readable
print(metadata)

# Add factor columns
metadata$tissue <- as.factor(ifelse(grepl("^T", metadata$type), "tumor", "normal"))
metadata$genotype <- as.factor(ifelse(grepl("WT$", metadata$type), "WT", "KO"))
metadata$animal_id <- as.factor(metadata$Pairs)

print(head(colnames(counts)))
print(head(metadata$sample))

print(length(colnames(counts)))
print(nrow(metadata))
print(which(colnames(counts) != metadata$sample))

print(data.frame(
  counts = colnames(counts)[14:28],
  metadata = metadata$sample[14:28]
))

# Reorder metadata to match counts column order
metadata <- metadata[match(colnames(counts), metadata$sample), ]
print(all(colnames(counts) == metadata$sample))

# Filter low-expression genes
keep <- filterByExpr(counts, group = metadata$genotype)
counts_filtered <- counts[keep, ]
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# PCA on all samples
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = metadata,
  design    = ~ genotype + tissue
)
vsd_all <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd_all, intgroup = c("tissue", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot_all <- ggplot(pca_data, aes(PC1, PC2, color = tissue, shape = genotype)) +
  geom_point(size = 5, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - All Samples")

print(pca_plot_all)
ggsave("results/plots/01_PCA_all_samples.pdf", plot = pca_plot_all, width = 8, height = 6)

# Remove outliers (low sequencing depth: ~3-4M reads vs ~12-18M in other samples)
outliers <- c("IonCode_0103", "IonCode_0105")
counts_clean <- counts_filtered[, !colnames(counts_filtered) %in% outliers]
metadata_clean <- metadata[metadata$sample %in% colnames(counts_clean), ]

cat("Samples after outlier removal:", ncol(counts_clean), "\n")
print(table(metadata_clean$genotype, metadata_clean$tissue))

# Create clean data directory if it doesn't exist
dir.create("data/clean", showWarnings = FALSE)

# Save clean data
write.csv(counts_clean, "data/clean/filtered_counts_clean.csv")
write.csv(metadata_clean, "data/clean/metadata_clean.csv")

# PCA after outlier removal
dds_clean <- DESeqDataSetFromMatrix(
  countData = counts_clean,
  colData   = metadata_clean,
  design    = ~ genotype + tissue
)
vsd_clean <- vst(dds_clean, blind = TRUE)

pca_clean <- plotPCA(vsd_clean, intgroup = c("tissue", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pca_clean, "percentVar"))

pca_plot_clean <- ggplot(pca_clean, aes(PC1, PC2, color = tissue, shape = genotype)) +
  geom_point(size = 5, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - After Outlier Removal")

print(pca_plot_clean)
ggsave("results/plots/02_PCA_clean.pdf", plot = pca_plot_clean, width = 8, height = 6)

print("Script 02 completed successfully.")