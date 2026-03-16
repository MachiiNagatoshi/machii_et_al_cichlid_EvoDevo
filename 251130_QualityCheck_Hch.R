# # =========================================================================
# # INSTALL
# # =========================================================================
# # ライブラリを削除
# remove.packages("Seurat")
# remove.packages("leidenbase")
# remove.packages("ggplot2")
# remove.packages("dplyr")
# remove.packages("patchwork")
# remove.packages("Matrix")
# remove.packages("SeuratObject")
# 
# 
# # 依存パッケージも含めて完全にクリーンアップしたい場合は、
# # Rを再起動してから以下を実行することをお勧めします
# 
# # ライブラリを再インストール
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("patchwork")
# install.packages("Matrix")
# 
# # Seuratのインストール　バイナリ版をインストール（コンパイル不要）
# install.packages("Seurat", type = "binary")
# 
# # 必要な依存パッケージ
# BiocManager::install("SingleCellExperiment")
# 
# # BiocManagerが必要
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # scDblFinderをインストール
# BiocManager::install("scDblFinder")



# =========================================================================
# scDblFinder Analysis (Seurat v5 compatible)
# =========================================================================
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

# Parameters
SAMPLE_NAME <- "Hch"  # または "Hsa"
CELLRANGER_DIR <- paste0("/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/CellRanger/", 
                         SAMPLE_NAME, "_outs/filtered_feature_bc_matrix")
OUTPUT_DIR <- paste0("/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/scRNAseq/DoubletFinder_", 
                     SAMPLE_NAME)
EXPECTED_DOUBLET_RATE <- 0.10

# Create output directorys
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
cts <- Read10X(data.dir = CELLRANGER_DIR)
seurat_obj_Hch <- CreateSeuratObject(counts = cts, project = SAMPLE_NAME, 
                                     min.cells = 3, min.features = 200)

# QC
seurat_obj_Hch[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_Hch, pattern = "^APK84-")


p4 <- FeatureScatter(seurat_obj_Hch, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(seurat_obj_Hch, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p4 + p5

ggsave(paste0(OUTPUT_DIR, "/FeatureScatter.png"), p4 + p5, width = 8, height = 6)

# percent.mtの外れ値を除去
seurat_obj_Hch <- subset(seurat_obj_Hch, subset = percent.mt < 10) 


# -------------------------------------------------------------------------
# Doublet detection with scDblFinder
# -------------------------------------------------------------------------
cat("Detecting doublets with scDblFinder...\n")

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj_Hch)

# Run scDblFinder
sce <- scDblFinder(sce, dbr = EXPECTED_DOUBLET_RATE)

# Add results back to Seurat object
seurat_obj_Hch$scDblFinder.class <- sce$scDblFinder.class
seurat_obj_Hch$scDblFinder.score <- sce$scDblFinder.score

# Summary
cat("\nDoublet detection summary:\n")
table(seurat_obj_Hch$scDblFinder.class)
cat("\nDoublet rate:", 
    sum(seurat_obj_Hch$scDblFinder.class == "doublet") / ncol(seurat_obj_Hch), "\n")


# -------------------------------------------------------------------------
# Filter doublets and save
# -------------------------------------------------------------------------
cat("Step 6: Filtering doublets...\n")

# Keep only singlets
seurat_singlets <- subset(seurat_obj_Hch, subset = scDblFinder.class == "singlet")

cat("Cells before filtering:", ncol(seurat_obj_Hch), "\n")
cat("Cells after filtering:", ncol(seurat_singlets), "\n")
cat("Removed:", ncol(seurat_obj_Hch) - ncol(seurat_singlets), "doublets\n")

# Save results
saveRDS(seurat_obj_Hch, paste0(OUTPUT_DIR, "/seurat_with_doublets.rds"))
saveRDS(seurat_singlets, paste0(OUTPUT_DIR, "/seurat_singlets_only.rds"))

cat("\nAnalysis complete! Results saved to:", OUTPUT_DIR, "\n")


# -------------------------------------------------------------------------
# Step 7: Check the result
# -------------------------------------------------------------------------


# Standard workflow
cat("Preprocessing...\n")
seurat_obj_Hch <- NormalizeData(seurat_obj_Hch)
seurat_obj_Hch <- FindVariableFeatures(seurat_obj_Hch)
seurat_obj_Hch <- ScaleData(seurat_obj_Hch)
seurat_obj_Hch <- RunPCA(seurat_obj_Hch, npcs = 50, verbose = FALSE)

# Clustering
cat("Clustering...\n")
seurat_obj_Hch <- FindNeighbors(seurat_obj_Hch, dims = 1:20)
seurat_obj_Hch <- FindClusters(seurat_obj_Hch, resolution = 0.8)
seurat_obj_Hch <- RunUMAP(seurat_obj_Hch, dims = 1:20)


# UMAP colored by doublet classification
p1 <- DimPlot(seurat_obj_Hch, reduction = "umap", group.by = "scDblFinder.class",
              cols = c("singlet" = "gray80", "doublet" = "red")) +
  ggtitle("Doublet Detection Results")
p1

# UMAP colored by doublet score
p2 <- FeaturePlot(seurat_obj_Hch, features = "scDblFinder.score", pt.size = 0.3) +
  scale_color_viridis_c() +
  ggtitle("Doublet Score")
p2

# Doublet score by cluster
p3 <- VlnPlot(seurat_obj_Hch, features = "scDblFinder.score", 
              group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("Doublet Score by Cluster")
p3

# Save plots
ggsave(paste0(OUTPUT_DIR, "/doublet_UMAP.png"), p1, width = 8, height = 6)
ggsave(paste0(OUTPUT_DIR, "/doublet_score_UMAP.png"), p2, width = 8, height = 6)
ggsave(paste0(OUTPUT_DIR, "/doublet_score_violin.png"), p3, width = 10, height = 6)


# Doublet detection summary
cat("\n=== Doublet Detection Summary ===\n")
cat("Total cells:", ncol(seurat_obj_Hch), "\n")
cat("Singlets:", sum(seurat_obj_Hch$scDblFinder.class == "singlet"), "\n")
cat("Doublets:", sum(seurat_obj_Hch$scDblFinder.class == "doublet"), "\n")

doublet_rate <- sum(seurat_obj_Hch$scDblFinder.class == "doublet") / ncol(seurat_obj_Hch) * 100
cat("Doublet rate:", round(doublet_rate, 2), "%\n")
cat("Expected rate:", EXPECTED_DOUBLET_RATE * 100, "%\n")