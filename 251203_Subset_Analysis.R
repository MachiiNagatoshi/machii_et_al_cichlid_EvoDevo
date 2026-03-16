# =========================================================================
# Subset Analysis: Fibroblast and Epithelial cells
# =========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)

# =========================================================================
# Paths
# =========================================================================
output_dir_epi <- "/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/Detail_epithelial"
output_dir_fib <- "/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/Detail_fibroblast"

dir.create(output_dir_epi, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_fib, showWarnings = FALSE, recursive = TRUE)

# =========================================================================
# Candidate genes (defined once)
# =========================================================================
known_markers <- list(
  vcanb = c("vcanb"),
  POSTN = c("POSTN"),
  col1a1a = c("ENSMZEG00005018722"),
  col1a1b = c("col1a1b"),
  col1a2 = c("col1a2")
)

# =========================================================================
# Load data
# =========================================================================
seurat_obj <- readRDS("/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/Integrated_Analysis/integrated_with_markers.rds")

# =========================================================================
# Function for subset analysis
# =========================================================================
analyze_subset <- function(obj, cell_type, output_dir) {
  cat("\n=== Analyzing", cell_type, "===\n")
  
  subset_obj <- subset(obj, idents = cell_type)
  cat("Cells extracted:", ncol(subset_obj), "\n")
  
  subset_obj <- NormalizeData(subset_obj)
  subset_obj <- FindVariableFeatures(subset_obj)
  subset_obj <- ScaleData(subset_obj)
  subset_obj <- RunPCA(subset_obj, npcs = 30)
  subset_obj <- FindNeighbors(subset_obj, dims = 1:20)
  subset_obj <- FindClusters(subset_obj, resolution = 0.5)
  subset_obj <- RunUMAP(subset_obj, dims = 1:20)
  
  # UMAP plots
  p1 <- DimPlot(subset_obj, reduction = "umap", label = TRUE) + 
    ggtitle(paste(cell_type, "- Clusters"))
  ggsave(paste0(output_dir, "/umap_clusters.svg"), p1, width = 8, height = 6)
  
  p2 <- DimPlot(subset_obj, reduction = "umap", group.by = "sample") + 
    ggtitle(paste(cell_type, "- by Sample"))
  ggsave(paste0(output_dir, "/umap_sample.svg"), p2, width = 8, height = 6)
  
  # DEGs
  markers <- FindAllMarkers(subset_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, paste0(output_dir, "/cluster_markers.csv"), row.names = FALSE)
  
  # Sample proportion per cluster
  prop_table <- table(subset_obj$seurat_clusters, subset_obj$sample)
  prop_pct <- prop.table(prop_table, margin = 1) * 100
  
  write.csv(prop_table, paste0(output_dir, "/cluster_sample_counts.csv"))
  write.csv(prop_pct, paste0(output_dir, "/cluster_sample_proportions.csv"))
  
  cat("\nSample proportions per cluster (%):\n")
  print(round(prop_pct, 2))
  
  saveRDS(subset_obj, paste0(output_dir, "/subset_object.rds"))
  cat("Results saved to:", output_dir, "\n")
  
  return(subset_obj)
}

# =========================================================================
# Run subset analysis
# =========================================================================
epi_obj <- analyze_subset(seurat_obj, "Epithelial cell", output_dir_epi)
fib_obj <- analyze_subset(seurat_obj, "Fibroblast", output_dir_fib)

cat("\n=== Subset Analysis Complete ===\n")

# =========================================================================
# Visualize candidate genes (FeaturePlot)
# =========================================================================
for(cell_type in names(known_markers)) {
  cat("\nProcessing", cell_type, "markers...\n")
  
  available <- intersect(known_markers[[cell_type]], rownames(seurat_obj))
  
  if(length(available) > 0) {
    p_fib <- FeaturePlot(fib_obj, features = available)
    p_epi <- FeaturePlot(epi_obj, features = available)
    
    ggsave(paste0(output_dir_fib, "/candidate_genes_", cell_type, "_featureplot.png"), 
           p_fib, width = 6, height = 5)
    ggsave(paste0(output_dir_epi, "/candidate_genes_", cell_type, "_featureplot.png"), 
           p_epi, width = 6, height = 5)
  } else {
    cat("  No markers found for", cell_type, "\n")
  }
}

# =========================================================================
# Correlation analysis
# =========================================================================
analyze_gene_correlation <- function(seurat_obj, target_gene, output_dir, prefix) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  if(!target_gene %in% rownames(seurat_obj)) {
    cat("ERROR:", target_gene, "not found\n")
    return(NULL)
  }
  
  target_expr <- FetchData(seurat_obj, vars = target_gene)[,1]
  expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  
  results <- apply(expr_matrix, 1, function(gene_expr) {
    test <- cor.test(target_expr, gene_expr, method = "spearman", exact = FALSE)
    c(cor = test$estimate, p = test$p.value)
  })
  
  cor_df <- data.frame(
    gene = colnames(results),
    correlation = results["cor.rho", ],
    p.value = results["p", ]
  ) %>%
    filter(gene != target_gene) %>%
    mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    arrange(desc(correlation))
  
  write.csv(cor_df, 
            paste0(output_dir, "/", prefix, "_", target_gene, "_correlation.csv"), 
            row.names = FALSE)
  
  top_genes <- c(target_gene, head(cor_df %>% filter(FDR < 0.05), 5)$gene)
  
  if(length(top_genes) > 1) {
    p <- FeaturePlot(seurat_obj, features = top_genes, ncol = 3)
    ggsave(paste0(output_dir, "/", prefix, "_", target_gene, "_featureplot.png"), 
           p, width = 12, height = 8)
  }
  
  cat(target_gene, ":", sum(cor_df$FDR < 0.05), "significant genes (FDR < 0.05)\n")
  
  return(cor_df)
}

# --- Run correlation analysis ---
genes_of_interest <- c("c1galt1lb")

cat("\n### Epithelial ###\n")
for(gene in genes_of_interest) {
  analyze_gene_correlation(epi_obj, gene, output_dir_epi, "epi")
}

cat("\nCorrelation Analysis Complete.\n")

# =========================================================================
# Blend plot (Fibroblast)
# =========================================================================
gene_pairs <- list(
  c("vcanb", "POSTN"),
  c("vcanb", "TCF4"),
  c("vcanb", "tnfaip6"),
  c("vcanb", "PTX3")
)

cat("\n### Fibroblast Blend Plots ###\n")
DefaultAssay(fib_obj) <- "RNA"

for(pair in gene_pairs) {
  gene1 <- pair[1]
  gene2 <- pair[2]
  
  if(!gene1 %in% rownames(fib_obj) | !gene2 %in% rownames(fib_obj)) {
    cat("Skipping:", gene1, "vs", gene2, "(not found)\n")
    next
  }
  
  p <- FeaturePlot(fib_obj, 
                   features = c(gene1, gene2), 
                   blend = TRUE, 
                   blend.threshold = 0.5,
                   cols = c("#E0E0DF", "#B22235", "#2983BB"))
  
  filename <- paste0(output_dir_fib, "/fib_", gene1, "_vs_", gene2, "_blend.png")
  ggsave(filename, p, width = 15, height = 4)
  cat("Saved:", filename, "\n")
}

cat("\n=== All Done ===\n")