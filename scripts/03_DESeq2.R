library(DESeq2)
library(dplyr)

counts_clean <- read.csv("data/clean/filtered_counts_clean.csv", row.names = 1)
metadata_clean <- read.csv("data/clean/metadata_clean.csv", row.names = 1)
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

filter_data <- function(
  test_type,
  reference_type,
  paired
) {
  # filter the data according to the types we want to compare
  metadata_sub <- metadata_clean |>
    filter(type %in% c(test_type, reference_type))

  # convert type into factor ("category")
  metadata_sub$type <- factor(metadata_sub$type)

  # make reference_type the reference group
  metadata_sub$type <- relevel(metadata_sub$type, ref = reference_type)

  if (paired) {
    # keep only animals that have both sample types
    complete_animals <- metadata_sub |>
      group_by(animal_id) |>
      filter(
        n() == 2,
        n_distinct(type) == 2
      ) |>
      pull(animal_id) |>
      unique()

    metadata_sub <- metadata_sub |>
      filter(animal_id %in% complete_animals)

    # blocked by animal_id if paired
    metadata_sub$animal_id <- factor(metadata_sub$animal_id)
    design_formula <- ~ animal_id + type
  } else {
    design_formula <- ~type
  }

  counts <- counts_clean[, metadata_sub$sample]

  # failsafe if for any reason metadata and counts dont align
  stopifnot(all(colnames(counts) == metadata_sub$sample))

  rownames(metadata_sub) <- metadata_sub$sample

  list(
    filtered_samples = metadata_sub,
    counts = counts,
    design_formula = design_formula,
    contrast = c("type", test_type, reference_type)
  )
}

run_deseq <- function(filtered_data) {
  dds <- DESeqDataSetFromMatrix(
    countData = filtered_data$counts,
    colData = filtered_data$filtered_samples,
    design = filtered_data$design_formula,
  )

  # run DESeq2
  dds <- DESeq(dds)

  # extract results for test_type vs reference_type
  res <- results(
    dds,
    contrast = filtered_data$contrast
  )

  # convert to data frame and add gene names
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  res_df <- res_df |>
    select(gene, everything())

  # filter DEGs
  degs <- res_df |>
    filter(!is.na(padj)) |>
    filter(padj < 0.05, abs(log2FoldChange) >= 1) |>
    mutate(
      direction = ifelse(log2FoldChange > 0, "up", "down")
    ) |>
    select(gene, everything()) |>
    arrange(padj)

  # create summary of up/down regulated genes
  summary <- data.frame(
    comparison = paste(
      filtered_data$contrast[2],
      "vs",
      filtered_data$contrast[3]
    ),
    upregulated = sum(degs$direction == "up"),
    downregulated = sum(degs$direction == "down"),
    total_degs = nrow(degs)
  )

  list(
    degs = degs,
    summary = summary
  )
}

# run the 4 comparisons
comparison_1_data <- filter_data("T_KO", "T_WT", paired = FALSE)
comparison_1 <- run_deseq(comparison_1_data)

comparison_2_data <- filter_data("N_KO", "N_WT", paired = FALSE)
comparison_2 <- run_deseq(comparison_2_data)

comparison_3_data <- filter_data("T_KO", "N_KO", paired = TRUE)
comparison_3 <- run_deseq(comparison_3_data)

comparison_4_data <- filter_data("T_WT", "N_WT", paired = TRUE)
comparison_4 <- run_deseq(comparison_4_data)

# write the deg results
write.csv(
  comparison_1$degs,
  "results/deseq2/Comparison_1_T_KO_vs_T_WT_DEGs.csv",
  row.names = FALSE
)

write.csv(
  comparison_2$degs,
  "results/deseq2/Comparison_2_N_KO_vs_N_WT_DEGs.csv",
  row.names = FALSE
)

write.csv(
  comparison_3$degs,
  "results/deseq2/Comparison_3_T_KO_vs_N_KO_DEGs.csv",
  row.names = FALSE
)

write.csv(
  comparison_4$degs,
  "results/deseq2/Comparison_4_T_WT_vs_N_WT_DEGs.csv",
  row.names = FALSE
)

# up/down regulation summary across different comparisons
up_down_summary <- bind_rows(
  comparison_1$summary,
  comparison_2$summary,
  comparison_3$summary,
  comparison_4$summary
)

write.csv(
  up_down_summary,
  "results/deseq2/Up_Down_DEGs_Summary_Table.csv",
  row.names = FALSE
)

cat("\nSample counts per comparison:\n")

cat("\nComparison 1: T_KO vs T_WT\n")
print(table(comparison_1_data$filtered_samples$type))

cat("\nComparison 2: N_KO vs N_WT\n")
print(table(comparison_2_data$filtered_samples$type))

cat("\nComparison 3: T_KO vs N_KO, paired\n")
print(table(comparison_3_data$filtered_samples$type))
cat("Number of paired animals:", length(unique(comparison_3_data$filtered_samples$animal_id)), "\n")

cat("\nComparison 4: T_WT vs N_WT, paired\n")
print(table(comparison_4_data$filtered_samples$type))
cat("Number of paired animals:", length(unique(comparison_4_data$filtered_samples$animal_id)), "\n")

cat("\nDesign formulas:\n")
cat("Comparison 1: "); print(comparison_1_data$design_formula)
cat("Comparison 2: "); print(comparison_2_data$design_formula)
cat("Comparison 3: "); print(comparison_3_data$design_formula)
cat("Comparison 4: "); print(comparison_4_data$design_formula)

cat("\nContrasts:\n")
print(comparison_1_data$contrast)
print(comparison_2_data$contrast)
print(comparison_3_data$contrast)
print(comparison_4_data$contrast)

cat("\nUp/Down DEG summary:\n")
print(up_down_summary)
