# =========================================================================
# Integration Analysis of Hsa and Hch samples (Seurat v5)
# =========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)

# =========================================================================
# Paths
# =========================================================================
OUTPUT_DIR <- "/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/Integrated_Analysis"
MARKER_DIR <- paste0(OUTPUT_DIR, "/Markers")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MARKER_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================================================================
# 1. Load singlet data
# =========================================================================
cat("Loading Hsa and Hch singlet data...\n")
hsa <- readRDS("/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/DoubletFinder_Hsa/seurat_singlets_only.rds")
hch <- readRDS("/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/DoubletFinder_Hch/seurat_singlets_only.rds")

cat("Hsa cells:", ncol(hsa), "\n")
cat("Hch cells:", ncol(hch), "\n")

# Add sample identity
hsa$sample <- "Hsa"
hch$sample <- "Hch"

# =========================================================================
# 2. CCA Integration
# =========================================================================
cat("\n=== CCA Integration ===\n")

merged <- merge(hsa, y = hch, add.cell.ids = c("Hsa", "Hch"), project = "Cichlid_Lip")

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

merged <- IntegrateLayers(
  object = merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.8)
merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:30)

p1 <- DimPlot(merged, reduction = "umap", group.by = "sample") + 
  ggtitle("CCA Integration by Sample")
p2 <- DimPlot(merged, reduction = "umap", label = TRUE) + 
  ggtitle("CCA Integration by Cluster")
ggsave(paste0(OUTPUT_DIR, "/UMAP_CCA_integration.svg"), p1 + p2, width = 14, height = 6)

saveRDS(merged, paste0(OUTPUT_DIR, "/integrated_CCA.rds"))

# =========================================================================
# 3. Prepare for marker detection
# =========================================================================
seurat_obj <- JoinLayers(merged)

# =========================================================================
# 4. Find markers for all clusters
# =========================================================================
cat("\nFinding markers for all clusters...\n")

all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(all_markers, paste0(MARKER_DIR, "/all_cluster_markers.csv"), row.names = FALSE)

# =========================================================================
# 5. Heatmap of top markers
# =========================================================================
cat("\nCreating heatmap of top markers...\n")

top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

p_heatmap <- DoHeatmap(seurat_obj, features = top10$gene, size = 3) + 
  NoLegend()
ggsave(paste0(MARKER_DIR, "/top10_markers_heatmap.png"), 
       p_heatmap, width = 12, height = 16)

# =========================================================================
# 6. Known marker genes - FeaturePlot per cell type
# =========================================================================
cat("\nChecking known marker genes...\n")

known_markers <- list(
  Lymphoid_1 = c("ENSMZEG00005006170", "CXCR4", "ccr9a"),
  Epithelial = c("epcam.1", "TGM1", "krt4", "DSP", "DSP.1", "dsc2l", "ppl", "evpla", "pkp3a", "klf4", "tp63", "itgb4", "itga6b", "LAMC2"),
  Fibroblasts = c("sparc", "dcn", "pdgfra", "prrx1b"),
  Myeloid = c("mpeg1.1", "irf8", "csf1ra", "ENSMZEG00005018567", "cyba", "cybb", "ncf1", "AIF1", "mpx", "mmp9"),
  Proliferating = c("pcna", "cdk1", "ccna2", "ccnb1", "smc2", "smc4"),
  Lymphoid_2 = c("cd79a", "ENSMZEG00005006303"),
  Chondrocyte = c("acanb", "COL2A1"),
  Erythroid = c("HBE1.2", "HBE1.1", "hbae4"),
  Osteoblast = c("runx2b", "sp7"),
  Goblet = c("agr2", "SPDEF"),
  Endothelial = c("pecam1a", "plvapa", "egfl7")
)

for(cell_type in names(known_markers)) {
  cat("\nProcessing", cell_type, "markers...\n")
  
  available <- intersect(known_markers[[cell_type]], rownames(seurat_obj))
  
  if(length(available) > 0) {
    cat("  Found", length(available), "markers:", paste(available, collapse=", "), "\n")
    
    n_markers <- length(available)
    ncol_plot <- min(4, n_markers)
    nrow_plot <- ceiling(n_markers / ncol_plot)
    plot_height <- max(4, nrow_plot * 3)
    plot_width <- ncol_plot * 3
    
    p <- FeaturePlot(seurat_obj, features = available, ncol = ncol_plot)
    ggsave(paste0(MARKER_DIR, "/known_markers_", cell_type, "_featureplot.png"),
           p, width = plot_width, height = plot_height)
  } else {
    cat("  No markers found for", cell_type, "\n")
  }
}

cat("\nAll known marker plots saved to:", MARKER_DIR, "\n")

# DotPlot
cat("\nCreating dot plot...\n")
p_dot <- DotPlot(seurat_obj, features = known_markers) + 
  RotatedAxis() +
  ggtitle("known markers dotplot")
ggsave(paste0(MARKER_DIR, "/known_markers_dotplot_test.svg"), 
       p_dot, width = 14, height = 8)

# =========================================================================
# 7. Cluster annotation
# =========================================================================
Idents(seurat_obj) <- "seurat_clusters"

new.cluster.ids <- c(
  "Lymphoid cell 1 (Tcell?)",
  "Lymphoid cell 1 (Tcell?)",
  "Epithelial cell",
  "Fibroblast",
  "Epithelial cell",
  "Epithelial cell",
  "Fibroblast",
  "Epithelial cell",
  "Fibroblast",
  "Fibroblast",
  "Epithelial cell",
  "Myeloid cell 1",
  "Proliferating cell(epithelial?)",
  "Lymphoid cell 2",
  "Myeloid cell 2",
  "Epithelial cell",
  "Chondrocyte",
  "Fibroblast",
  "Lymphoid cell 1 (Tcell?)",
  "Myeloid cell 3 (Neutrophil?)",
  "Erythroid cell",
  "Osteoblast",
  "Epithelial cell",
  "Goblet cell",
  "Myeloid cell 4",
  "Endothelial cell"
)

names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$cell_type <- Idents(seurat_obj)

p_annot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave(paste0(MARKER_DIR, "/umap_annotated_labels.svg"), p_annot, width = 12, height = 10)

p_split <- DimPlot(seurat_obj, reduction = "umap", split.by = "sample")
ggsave(paste0(MARKER_DIR, "/umap_annotated_split_sample.svg"), p_split, width = 20, height = 8)

# =========================================================================
# 8. Save final object
# =========================================================================
saveRDS(seurat_obj, paste0(OUTPUT_DIR, "/integrated_with_markers.rds"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", MARKER_DIR, "\n")
cat("Key outputs:\n")
cat("  - all_cluster_markers.csv: All cluster markers\n")
cat("  - top10_markers_heatmap.png: Top 10 markers heatmap\n")
cat("  - Various visualization plots\n")