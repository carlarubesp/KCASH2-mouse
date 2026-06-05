# ==============================================================================
# SCRIPT 04: Visualization
# ==============================================================================

library(EnhancedVolcano)
library(pheatmap)
library(dplyr)
library(ggplot2)

# Create the directory for saving plots automatically if it doesn't exist
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# LOAD DATA 

deg1 <- read.csv("results/deseq2/Comparison_1_T_KO_vs_T_WT_DEGs.csv")
deg2 <- read.csv("results/deseq2/Comparison_2_N_KO_vs_N_WT_DEGs.csv")
deg3 <- read.csv("results/deseq2/Comparison_3_T_KO_vs_N_KO_DEGs.csv")
deg4 <- read.csv("results/deseq2/Comparison_4_T_WT_vs_N_WT_DEGs.csv")

# VOLCANO PLOTS

# Check if the dataset has no significant genes to plot
make_volcano <- function(df, title, file) {
  
  if (nrow(df) == 0) {
    p <- ggplot() + 
      scale_x_continuous(limits = c(-4, 4), name = expression(Log[2] ~ fold ~ change)) +
      scale_y_continuous(limits = c(0, 6), name = expression(-Log[10] ~ P)) +
      annotate("text", x = 0, y = 3, 
               label = "No Significant DEGs Found\n(padj < 0.05, |log2FC| >= 1)", 
               size = 5, color = "red", fontface = "bold") + 
      theme_classic(base_size = 15) + 
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray95")
      ) +
      labs(title = title, subtitle = "Total Significant DEGs: 0")
    
    ggsave(file, plot = p, width = 8, height = 6, device = "pdf")
    cat("Generated empty-axes volcano for", title, "- No genes to plot.\n")
    return(NULL)
  }
  
  # Generate the full-scale global volcano plot
  p <- EnhancedVolcano(
    df,
    lab = df$gene,
    x = "log2FoldChange",
    y = "padj", 
    pCutoff = 0.05,
    FCcutoff = 1,
    title = title,
    subtitle = paste("Total Significant DEGs:", nrow(df)),
    legendPosition = "right",
    pointSize = 2.0,
    labSize = 4.0
  )
  ggsave(file, plot = p, width = 8, height = 6, device = "pdf")
  
  
  cat("Successfully generated Volcano plot for", title, "- Saved as PDF.\n")
}

# Execute the volcano plotting function for the four comparisons
make_volcano(deg1, "Comparison 1: T_KO vs T_WT", "results/plots/Comparison_1_T_KO_vs_T_WT_DEGs.pdf")
make_volcano(deg2, "Comparison 2: N_KO vs N_WT", "results/plots/Comparison_2_N_KO_vs_N_WT_DEGs.pdf")
make_volcano(deg3, "Comparison 3: T_KO vs N_KO", "results/plots/Comparison_3_T_KO_vs_N_KO_DEGs.pdf")
make_volcano(deg4, "Comparison 4: T_WT vs N_WT", "results/plots/Comparison_4_T_WT_vs_N_WT_DEGs.pdf")


# 3. HEATMAP FOR TOP 50 GENES

counts <- read.csv("data/clean/filtered_counts_clean.csv", row.names = 1)
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

top50_genes <- deg3 %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(50) %>%
  pull(gene)

common_genes <- intersect(top50_genes, rownames(counts))

if (length(common_genes) > 0) {
  heatmap_matrix <- counts[common_genes, , drop = FALSE]
  heatmap_matrix <- log2(heatmap_matrix + 1) 
  
  pheatmap::pheatmap(
    heatmap_matrix,
    scale = "row", 
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = "Top 50 Most Significant DEGs (Tumor KO vs Normal KO)",
    fontsize_row = 8,
    fontsize_col = 9,
    filename = "results/plots/Heatmap_Comparison_3_T_KO_vs_N_KO_DEGs.pdf", 
    width = 10,  
    height = 12
  )
  
  cat("Heatmap generated and saved successfully as PDF.\n")
} else {
  cat("Warning: No common genes found for Heatmap.\n")
}


# 4. TARGETED PATHWAY ANALYSIS: Hedgehog & Gli Signaling

# Define the core Hedgehog pathway components
# Note: Kctd21 is the official mouse gene symbol for the tumor suppressor KCASH2
hedgehog_target_genes <- c("Gli1", "Gli2", "Gli3", "Ptch1", "Smo", "Shh", "Kctd21")

# Intersect target genes to verify their presence within the clean count matrix rows
hh_common_genes <- intersect(hedgehog_target_genes, rownames(counts))

if (length(hh_common_genes) > 0) {
  hh_matrix <- counts[hh_common_genes, , drop = FALSE]
  hh_matrix <- log2(hh_matrix + 1)
  
  # Generate the customized targeted pathway expression profile
  pheatmap::pheatmap(
    hh_matrix,
    scale = "row", 
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = "Expression Profile of Key Hedgehog Pathway Genes (Targeted Analysis)",
    fontsize_row = 10,
    fontsize_col = 9,
    filename = "results/plots/Targeted_Hedgehog_Pathway_Heatmap.pdf", 
    width = 10,  
    height = 5
  )
  
  cat("Targeted Hedgehog Pathway Heatmap generated and saved successfully as PDF.\n")
} else {
  cat("Warning: Target Hedgehog pathway genes were not found in the clean counts matrix.\n")
}

print("Script 04 completed successfully.")