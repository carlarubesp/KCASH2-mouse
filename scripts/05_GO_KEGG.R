
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
})

# ----------------------------------------------------------
# 1. Define directories
# ----------------------------------------------------------

results_dir <- "results"

# IMPORTANT:
# This must match the output directory used by 03_DESeq2.R
de_dir <- file.path(results_dir, "deseq2")

func_dir <- file.path(results_dir, "functional")
plot_dir <- file.path(results_dir, "plots")

dir.create(func_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 2. Define official DEG files from 03_DESeq2.R
# ----------------------------------------------------------

deg_files <- c(
  "T_KO_vs_T_WT" = "Comparison_1_T_KO_vs_T_WT_DEGs.csv",
  "N_KO_vs_N_WT" = "Comparison_2_N_KO_vs_N_WT_DEGs.csv",
  "T_KO_vs_N_KO" = "Comparison_3_T_KO_vs_N_KO_DEGs.csv",
  "T_WT_vs_N_WT" = "Comparison_4_T_WT_vs_N_WT_DEGs.csv"
)

comparisons <- names(deg_files)

# ----------------------------------------------------------
# 3. Load DEG tables
# ----------------------------------------------------------

load_deg <- function(comp) {
  
  f <- file.path(de_dir, deg_files[[comp]])
  
  if (!file.exists(f)) {
    warning("No official DEG table found for: ", comp, " at ", f)
    return(tibble())
  }
  
  x <- read_csv(f, show_col_types = FALSE)
  
  # Standardise gene column name
  if (!"gene" %in% colnames(x)) {
    gene_col <- intersect(
      c("gene", "Gene", "gene_id", "symbol", "SYMBOL", "external_gene_name"),
      colnames(x)
    )[1]
    
    if (is.na(gene_col)) {
      warning("No gene column found in: ", f)
      return(tibble())
    }
    
    x <- rename(x, gene = all_of(gene_col))
  }
  
  x$comparison <- comp
  
  return(x)
}

deg_list <- setNames(lapply(comparisons, load_deg), comparisons)

# Save DEG counts used by this script
deg_counts <- tibble(
  comparison = names(deg_list),
  n_DEGs_used_for_enrichment = sapply(deg_list, nrow)
)

write_csv(
  deg_counts,
  file.path(func_dir, "DEG_counts_used_for_functional_analysis.csv")
)

print(deg_counts)

# ----------------------------------------------------------
# 4. Convert mouse gene symbols to Entrez IDs
# ----------------------------------------------------------

convert_to_entrez <- function(deg_df, comp) {
  
  if (nrow(deg_df) == 0) {
    warning("No DEGs available for Entrez conversion in ", comp)
    return(tibble())
  }
  
  symbols <- unique(na.omit(deg_df$gene))
  
  if (length(symbols) == 0) {
    warning("No valid gene symbols found for ", comp)
    return(tibble())
  }
  
  converted <- bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  converted <- converted %>%
    distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
  
  write_csv(
    converted,
    file.path(func_dir, paste0(comp, "_symbol_to_entrez.csv"))
  )
  
  return(converted)
}

entrez_list <- mapply(
  convert_to_entrez,
  deg_list,
  names(deg_list),
  SIMPLIFY = FALSE
)

# Save Entrez mapping counts
entrez_counts <- tibble(
  comparison = names(entrez_list),
  mapped_entrez_ids = sapply(entrez_list, nrow)
)

write_csv(
  entrez_counts,
  file.path(func_dir, "Entrez_mapping_counts.csv")
)

print(entrez_counts)

# ----------------------------------------------------------
# 5. GO Biological Process enrichment
# ----------------------------------------------------------

run_go <- function(entrez_df, comp) {
  
  if (nrow(entrez_df) < 10) {
    warning("Too few mapped DEGs for GO enrichment in ", comp)
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = unique(entrez_df$ENTREZID),
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  
  write_csv(
    ego_df,
    file.path(func_dir, paste0(comp, "_GO_BP.csv"))
  )
  
  if (nrow(ego_df) > 0) {
    
    p <- dotplot(ego, showCategory = 15) +
      ggtitle(paste("GO Biological Process:", comp)) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    ggsave(
      file.path(plot_dir, paste0(comp, "_GO_BP_dotplot.pdf")),
      p,
      width = 10,
      height = 8
    )
    
    ggsave(
      file.path(plot_dir, paste0(comp, "_GO_BP_dotplot.png")),
      p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  return(ego)
}

go_list <- mapply(
  run_go,
  entrez_list,
  names(entrez_list),
  SIMPLIFY = FALSE
)

# ----------------------------------------------------------
# 6. KEGG pathway enrichment
# ----------------------------------------------------------

run_kegg <- function(entrez_df, comp) {
  
  if (nrow(entrez_df) < 10) {
    warning("Too few mapped DEGs for KEGG enrichment in ", comp)
    return(NULL)
  }
  
  ekegg <- enrichKEGG(
    gene = unique(entrez_df$ENTREZID),
    organism = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  ekegg <- setReadable(
    ekegg,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID"
  )
  
  ekegg_df <- as.data.frame(ekegg)
  
  write_csv(
    ekegg_df,
    file.path(func_dir, paste0(comp, "_KEGG.csv"))
  )
  
  targeted <- ekegg_df %>%
    filter(
      grepl(
        "Hedgehog|Wnt|colorectal|cell cycle|TGF|MAPK|PI3K|cancer|pluripotency",
        Description,
        ignore.case = TRUE
      )
    )
  
  write_csv(
    targeted,
    file.path(func_dir, paste0(comp, "_KEGG_targeted_pathways.csv"))
  )
  
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
    
    ggsave(
      file.path(plot_dir, paste0(comp, "_KEGG_dotplot.pdf")),
      p,
      width = 12,
      height = 8
    )
    
    ggsave(
      file.path(plot_dir, paste0(comp, "_KEGG_dotplot.png")),
      p,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  return(ekegg)
}

kegg_list <- mapply(
  run_kegg,
  entrez_list,
  names(entrez_list),
  SIMPLIFY = FALSE
)

# ----------------------------------------------------------
# 7. compareCluster GO analysis
# ----------------------------------------------------------

entrez_clusters <- lapply(entrez_list, function(x) {
  if (!"ENTREZID" %in% colnames(x)) return(character(0))
  unique(x$ENTREZID)
})

entrez_clusters <- entrez_clusters[lengths(entrez_clusters) >= 10]

if (length(entrez_clusters) >= 2) {
  
  cc_go <- compareCluster(
    geneCluster = entrez_clusters,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  cc_go_df <- as.data.frame(cc_go)
  
  write_csv(
    cc_go_df,
    file.path(func_dir, "compareCluster_GO_BP.csv")
  )
  
  if (nrow(cc_go_df) > 0) {
    
    p <- dotplot(cc_go, showCategory = 10) +
      ggtitle("Functional comparison across DEG contrasts") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    ggsave(
      file.path(plot_dir, "compareCluster_GO_BP.pdf"),
      p,
      width = 11,
      height = 7
    )
    
    ggsave(
      file.path(plot_dir, "compareCluster_GO_BP.png"),
      p,
      width = 11,
      height = 7,
      dpi = 300
    )
  }
  
} else {
  warning("Not enough DEG clusters with >=10 mapped genes for compareCluster.")
}

# ----------------------------------------------------------
# 8. Specific Hedgehog gene check
# ----------------------------------------------------------

hh_genes <- c(
  "Shh", "Ihh", "Dhh",
  "Ptch1", "Ptch2",
  "Smo", "Sufu",
  "Gli1", "Gli2", "Gli3",
  "Hhip",
  "Kctd21"
)

hh_table <- bind_rows(lapply(names(deg_list), function(comp) {
  
  x <- deg_list[[comp]]
  
  if (nrow(x) == 0) return(tibble())
  
  x %>%
    filter(gene %in% hh_genes) %>%
    mutate(comparison = comp)
}))

write_csv(
  hh_table,
  file.path(func_dir, "Hedgehog_DEG_hits.csv")
)

# ----------------------------------------------------------
# 9. Targeted pathway screen from DEG gene lists
# ----------------------------------------------------------

targeted_gene_sets <- list(
  Hedgehog_signaling = c(
    "Shh", "Ihh", "Dhh", "Ptch1", "Ptch2", "Smo", "Sufu",
    "Gli1", "Gli2", "Gli3", "Hhip", "Kctd21"
  ),
  Wnt_beta_catenin = c(
    "Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt3a", "Wnt4", "Wnt5a", "Wnt5b",
    "Wnt6", "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b",
    "Wnt10a", "Wnt10b", "Wnt11", "Wnt16",
    "Fzd1", "Fzd2", "Fzd3", "Fzd4", "Fzd5", "Fzd6", "Fzd7", "Fzd8", "Fzd9", "Fzd10",
    "Ctnnb1", "Apc", "Axin1", "Axin2", "Dkk1", "Dkk2", "Dkk3", "Dkk4",
    "Tcf7", "Tcf7l1", "Tcf7l2", "Lef1"
  ),
  Inflammation_immune = c(
    "Tnf", "Il1b", "Il6", "Il10", "Ifng", "Cxcl1", "Cxcl2", "Cxcl10",
    "Ccl2", "Ccl5", "Nfkb1", "Nfkb2", "RelA", "Stat1", "Stat3",
    "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd19", "Cd68", "Adgre1"
  )
)

all_detected_genes <- unique(unlist(lapply(deg_list, function(x) x$gene)))

targeted_screen <- bind_rows(lapply(names(deg_list), function(comp) {
  
  deg_symbols <- unique(deg_list[[comp]]$gene)
  
  bind_rows(lapply(names(targeted_gene_sets), function(pathway_name) {
    
    pathway_genes <- targeted_gene_sets[[pathway_name]]
    
    pathway_detected <- intersect(pathway_genes, all_detected_genes)
    degs_in_pathway <- intersect(deg_symbols, pathway_genes)
    
    # Fisher exact test
    a <- length(degs_in_pathway)
    b <- length(deg_symbols) - a
    c <- length(pathway_detected) - a
    d <- length(all_detected_genes) - a - b - c
    
    if (any(c(a, b, c, d) < 0)) {
      fisher_p <- NA_real_
      odds <- NA_real_
    } else {
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      fisher_p <- ft$p.value
      odds <- unname(ft$estimate)
    }
    
    tibble(
      comparison = comp,
      pathway_screen = pathway_name,
      pathway_genes_detected_in_dataset = length(pathway_detected),
      DEGs_in_pathway = length(degs_in_pathway),
      DEG_symbols = paste(degs_in_pathway, collapse = ";"),
      fisher_pvalue = fisher_p,
      odds_ratio = odds
    )
  }))
}))

targeted_screen <- targeted_screen %>%
  group_by(comparison) %>%
  mutate(padj_within_screen = p.adjust(fisher_pvalue, method = "BH")) %>%
  ungroup()

write_csv(
  targeted_screen,
  file.path(func_dir, "targeted_pathway_screen.csv")
)

targeted_hits <- targeted_screen %>%
  filter(!is.na(padj_within_screen), padj_within_screen < 0.05)

write_csv(
  targeted_hits,
  file.path(func_dir, "targeted_pathway_hits_only.csv")
)

# ----------------------------------------------------------
# 10. Save session information for reproducibility
# ----------------------------------------------------------

sink(file.path(func_dir, "sessionInfo_block4.txt"))
print(sessionInfo())
sink()

message("Functional enrichment analysis completed successfully.")