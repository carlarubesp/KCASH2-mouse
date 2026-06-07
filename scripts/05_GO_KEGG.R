

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

results_dir <- "results"
de_dir <- file.path(results_dir, "DE_tables")
func_dir <- file.path(results_dir, "functional")
plot_dir <- file.path(results_dir, "plots")
dir.create(func_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- c("T_KO_vs_T_WT", "N_KO_vs_N_WT", "T_KO_vs_N_KO_paired", "T_WT_vs_N_WT_paired")

load_deg <- function(comp) {
  candidates <- c(file.path(de_dir, paste0(comp, "_DEGs.csv")), file.path(de_dir, paste0(comp, "_DESeq2_DEGs.csv")), file.path(de_dir, paste0(comp, "_DESeq2_results_DEGs.csv")))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) { warning("No DEG table found for: ", comp); return(tibble()) }
  x <- read_csv(f, show_col_types = FALSE)
  if (!"gene" %in% colnames(x)) {
    gene_col <- intersect(c("Gene", "gene_id", "symbol", "SYMBOL"), colnames(x))[1]
    if (!is.na(gene_col)) x <- rename(x, gene = all_of(gene_col))
  }
  x$comparison <- comp
  x
}

deg_list <- setNames(lapply(comparisons, load_deg), comparisons)

convert_to_entrez <- function(deg_df, comp) {
  if (nrow(deg_df) == 0) return(tibble())
  symbols <- unique(na.omit(deg_df$gene))
  converted <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  write_csv(converted, file.path(func_dir, paste0(comp, "_symbol_to_entrez.csv")))
  converted
}

entrez_list <- mapply(convert_to_entrez, deg_list, names(deg_list), SIMPLIFY = FALSE)

run_go <- function(entrez_df, comp) {
  if (nrow(entrez_df) < 10) { warning("Too few mapped DEGs for GO enrichment in ", comp); return(NULL) }
  ego <- enrichGO(gene = unique(entrez_df$ENTREZID), OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  write_csv(as.data.frame(ego), file.path(func_dir, paste0(comp, "_GO_BP.csv")))
  if (nrow(as.data.frame(ego)) > 0) {
    p <- dotplot(ego, showCategory = 20) + ggtitle(paste("GO Biological Process:", comp))
    ggsave(file.path(plot_dir, paste0(comp, "_GO_BP_dotplot.pdf")), p, width = 9, height = 7)
    ggsave(file.path(plot_dir, paste0(comp, "_GO_BP_dotplot.png")), p, width = 9, height = 7, dpi = 300)
  }
  ego
}

go_list <- mapply(run_go, entrez_list, names(entrez_list), SIMPLIFY = FALSE)

run_kegg <- function(entrez_df, comp) {
  if (nrow(entrez_df) < 10) { warning("Too few mapped DEGs for KEGG enrichment in ", comp); return(NULL) }
  ekegg <- enrichKEGG(gene = unique(entrez_df$ENTREZID), organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  ekegg_df <- as.data.frame(ekegg)
  write_csv(ekegg_df, file.path(func_dir, paste0(comp, "_KEGG.csv")))
  targeted <- ekegg_df %>% filter(grepl("Hedgehog|Wnt|colorectal|cell cycle|TGF|MAPK|PI3K", Description, ignore.case = TRUE))
  write_csv(targeted, file.path(func_dir, paste0(comp, "_KEGG_targeted_pathways.csv")))
  if (nrow(ekegg_df) > 0) {
    p <- dotplot(ekegg, showCategory = 10) +
      ggtitle(paste("KEGG pathways:", comp)) +
      scale_y_discrete(labels = function(x) {
        x <- gsub(" - Mus musculus \\(house mouse\\)", "", x)
        stringr::str_wrap(x, width = 35)
      }) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    ggsave(file.path(plot_dir, paste0(comp, "_KEGG_dotplot.pdf")), p, width = 12, height = 8)
    ggsave(file.path(plot_dir, paste0(comp, "_KEGG_dotplot.png")), p, width = 12, height = 8, dpi = 300)
  }
  ekegg
}

kegg_list <- mapply(run_kegg, entrez_list, names(entrez_list), SIMPLIFY = FALSE)

entrez_clusters <- lapply(entrez_list, function(x) unique(x$ENTREZID))
entrez_clusters <- entrez_clusters[lengths(entrez_clusters) >= 10]
if (length(entrez_clusters) >= 2) {
  cc_go <- compareCluster(geneCluster = entrez_clusters, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  write_csv(as.data.frame(cc_go), file.path(func_dir, "compareCluster_GO_BP.csv"))
  if (nrow(as.data.frame(cc_go)) > 0) {
    p <- dotplot(cc_go, showCategory = 10) + ggtitle("Functional comparison across DEG contrasts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(plot_dir, "compareCluster_GO_BP.pdf"), p, width = 11, height = 7)
    ggsave(file.path(plot_dir, "compareCluster_GO_BP.png"), p, width = 11, height = 7, dpi = 300)
  }
}

hh_genes <- c("Shh", "Ihh", "Dhh", "Ptch1", "Ptch2", "Smo", "Sufu", "Gli1", "Gli2", "Gli3", "Hhip", "Kctd21")
hh_table <- bind_rows(lapply(names(deg_list), function(comp) {
  x <- deg_list[[comp]]
  if (nrow(x) == 0) return(tibble())
  x %>% filter(gene %in% hh_genes) %>% mutate(comparison = comp)
}))
write_csv(hh_table, file.path(func_dir, "Hedgehog_DEG_hits.csv"))

sink(file.path(func_dir, "sessionInfo_block4.txt"))
print(sessionInfo())
sink()
